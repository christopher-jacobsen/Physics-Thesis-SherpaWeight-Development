////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaWeight.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeight.h"
#include "SherpaMECalculator.h"

// Sherpa includes
#include <SHERPA/Main/Sherpa.H>
#include <SHERPA/Initialization/Initialization_Handler.H>
#include <ATOOLS/Org/Data_Reader.H>
#include <ATOOLS/Math/Vector.H>
//#include <ATOOLS/Phys/Cluster_Amplitude.H>

// Root includes
#include <TMatrixD.h>
#include <TVectorD.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::SherpaWeight()
    : m_upSherpa( new SHERPA::Sherpa )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::~SherpaWeight() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::Initialize( const std::string & eventFileName, const std::vector<const char *> & argv )
{
    m_eventFileName = eventFileName;
    
    //m_argv.push_back( "INIT_ONLY=2" ); // prevent Sherpa from starting the cross section integration

    // TODO: check if already initialized
    
    // initialize sherpa

    if (!m_upSherpa->InitializeTheRun( static_cast<int>(argv.size()), const_cast<char **>(argv.data()) ))
    {
        LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        throw int(-2);  // TODO
    }
    
    // get the various configuration paths and file names
    
    SHERPA::Initialization_Handler * pInitHandler = m_upSherpa->GetInitHandler();
    if (!pInitHandler)
        throw int(-2); // TODO
    
    std::string runPath  = pInitHandler->Path();
    std::string runFile  = pInitHandler->File();
    std::string baseFile = runFile.substr( 0, runFile.find("|") );  // strip off section declaration following '|'

    // read the run file/section for optional SHERPA_WEIGHT_FILE definition
    std::string weightFile, paramFile;
    {
        ATOOLS::Data_Reader reader(" ",";","!","=");
        reader.AddComment("#");
        reader.AddWordSeparator("\t");
        reader.SetInputPath( runPath );
        reader.SetInputFile( runFile );
        
        // Note: By default, Data_Reader appends the (run) section and command line parameters to the section of the file it is reading.
      
        weightFile = reader.GetValue<std::string>( "SHERPA_WEIGHT_FILE", baseFile + "|(SherpaWeight){|}(SherpaWeight)" );
        paramFile  = reader.GetValue<std::string>( "FR_PARAMCARD",       std::string("param_card.dat") );
    }
    
    // read list of reweight parameters

    {
        ATOOLS::Data_Reader reader(" ",";","!","=");
        reader.AddComment("#");
        reader.AddWordSeparator("\t");
        reader.SetInputPath( runPath );
        reader.SetInputFile( weightFile );
        
        // Note: By default, Data_Reader appends the (run) section and command line parameters to the section of the file it is reading.
        reader.SetAddCommandLine(false);  // read the raw file input

        std::vector< std::vector<std::string> > stringMatrix;
        if (reader.MatrixFromFile( stringMatrix ))
        {
            m_parameters.reserve( stringMatrix.size() );    // reserve space for all entries
            
            for ( const auto & row : stringMatrix )
            {
                if (!row.size())    // skip empty rows (probably not necessary, but just for caution sake)
                    continue;
            
                ReweightParameter param;

                param.name = row[0];
                
                if (row.size() > 1)
                {
                    try
                    {
                        // could use ATOOLS::ToType<double>(), but that doesn't check if the conversion failed
                        // use std::stringstream instead
                    
                        std::stringstream stream;
                        stream.precision(12);
                    
                        stream << row[1];
                        stream >> param.scale;
                        if (stream.fail() || !stream.eof())
                            throw -2;  // TODO
                        
                        if (row.size() > 2)
                        {
                            stream.clear();
                            stream << row[2];
                            stream >> param.offset;
                            if (stream.fail() || !stream.eof())
                                throw -2;  // TODO
                        }
                    }
                    catch (const std::exception & error)
                    {
                        LogMsgError( "Failed to read scale/offset for reweight parameter %hs", FMT_HS(param.name.c_str()) );
                        throw;
                    }
                }
                
                m_parameters.push_back( std::move(param) );
            }
        }
        
        // TODO:: ensure all names are unique
    }
    
    if (NParameters())
    {
        GetBilinearMatrices( m_parameters, m_evalMatrix, m_invCoefMatrix, m_coefNames );
        
        // create list of param_card files
        
        std::string srcFilePath = runPath + paramFile;
        std::string dstPath     = runPath + "SherpaWeight/";

        // create destination param_card directory
        {
            std::string command = "mkdir \"" + dstPath + "\"";
            system( command.c_str() );
        }

        for (size_t i = 0; i < m_evalMatrix.size(); ++i)
        {
            // construct destination file path
            char dstFile[40];
            sprintf( dstFile, "param_card_%03u.dat", FMT_U(i) );  // %i is max 10 chars
            
            std::string dstFilePath( dstPath + dstFile );

            // create param_card with new parameter values
            CreateFeynRulesParamCard( srcFilePath, dstFilePath, m_parameters, m_evalMatrix[i] );

            m_paramCards.push_back( dstFilePath );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::EvaluateEvents()
{
    //std::string                 argParamCard( "FR_PARAMCARD=" + m_paramCards[run] );
    //std::vector<const char *>   runArgv(m_argv);

    // TODO: implement this method
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::GetBilinearMatrices( const ParameterVector & parameters,
                                        DoubleMatrix & evalMatrix,
                                        DoubleMatrix & invCoefMatrix,
                                        StringVector & coefNames )  // static
{
    // clear outputs
    evalMatrix    .clear();
    invCoefMatrix .clear();
    coefNames     .clear();

    size_t nParam = parameters.size();
    if (nParam < 1)
        return;

    size_t nCoefs = (nParam + 1) * (nParam + 2) / 2;  // nCoefs >= 3 [nCoefs = 3,6,10,15,...]

    // eval matrix  rows: nCoefs   columns: nParam
    // coef matrix  rows: nCoefs   columns: nCoefs
    double eval[nCoefs][nParam];
    double coef[nCoefs][nCoefs];
    
    memset( eval, 0, sizeof(eval) );    // zero all elements
    memset( coef, 0, sizeof(coef) );    // zero all elements

    //// fill in eval matrix
    {
        // eval[0] = [0,...,0]
        
        // eval[1..nParam] = diagonal of 1's
        for (size_t i = 0; i < nParam; ++i)
        {
            eval[i+1][i] = 1.0;
        }

        // eval[nParam+1..2*nParam] = diagonal of -1's
        for (size_t i = 0; i < nParam; ++i)
        {
            eval[i+1+nParam][i] = -1.0;
        }
        
        // remainder with rotations of [1,-1,0,...], [1,0,-1,0,...], ...
        for (size_t r = 2 * nParam + 1, g = 1; r < nCoefs; )
        {
            double row[nParam];
            
            memset( row, 0, sizeof(row) );  // zero all elements
            
            row[0]   =  1.0;
            row[g++] = -1.0;
            
            for (size_t i = 0; (i < nParam) && (r < nCoefs); ++i)
            {
                memcpy( eval[r++], row, sizeof(row) );
                std::rotate( row, row + nParam - 1, row + nParam );
            }
        }
        
        // multiply with scale and add offset
        
        for (size_t r = 0; r < nCoefs; ++r)
        {
            for (size_t c = 0; c < nParam; ++c)
            {
                const ReweightParameter & param = parameters[c];
                double &                  entry = eval[r][c];

                entry = entry * param.scale + param.offset;
            }
        }
    }

    // fill in coef matrix
    {
        // filling of Mvector
        size_t k;
        for (size_t r = 0; r < nCoefs; ++r)
        {
            coef[r][0] = 1.0;
            k = 1 + nParam;

            for (size_t c = 0; c < nParam; ++c)
            {
                double value   = eval[r][c];
                double sqValue = value * value;

                coef[r][c + 1] = value;     // single terms
                coef[r][k]     = sqValue;   // square terms
                
                k += nParam - c;
            }

            // interference terms
            k = 1 + nParam;
            for (size_t c = 0; c < nParam - 1; ++c)
            {
                for (size_t j = c + 1; j < nParam; ++j)
                {
                    double value = eval[r][c] * eval[r][j];
                
                    coef[r][k + j - c] = value;
                }

                k += nParam - c;
            }
        }
    }
    
    // invert coefficient matrix
    {
        TMatrixD matrix( (Int_t)nCoefs, (Int_t)nCoefs, static_cast<const Double_t *>( coef[0] ) );
        matrix.Invert();
        
        /*
        // get rid of numerical precision residuals
        for (size_t i = 0; i < nCoefs; ++i)
        {
            for (size_t j = 0; j < nCoefs; ++j)
            {
                if (abs(matrix[i][j]) < 1.0E-10)
                    matrix[i][j] = 0.0;
            }
        }
        */
        
        memcpy( coef[0], matrix.GetMatrixArray(), sizeof(coef) );
    }

    // fill in output matrices
    {
        evalMatrix   .resize( nCoefs );
        invCoefMatrix.resize( nCoefs );

        for (size_t r = 0; r < nCoefs; ++r)
        {
            evalMatrix   [r].assign( eval[r], eval[r] + nParam );
            invCoefMatrix[r].assign( coef[r], coef[r] + nCoefs );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::CreateFeynRulesParamCard( const std::string & srcFilePath, const std::string & dstFilePath,
                                             const ParameterVector & parameters, const DoubleVector & paramValues )    // static
{
    struct Local  // local object for automated cleanup
    {
        FILE * fpSrc = nullptr;
        FILE * fpDst = nullptr;

        ~Local()
        {
            if (fpSrc) { fclose(fpSrc); fpSrc = nullptr; }
            if (fpDst) { fclose(fpDst); fpDst = nullptr; }
        }
    }
    local;


    if (parameters.size() != paramValues.size())
        throw -4; // TODO
    
    local.fpSrc = fopen( srcFilePath.c_str(), "rt" );
    if (!local.fpSrc)
        throw -5;   // TODO
    
    local.fpDst = fopen( dstFilePath.c_str(), "wt" );
    if (!local.fpDst)
        throw -5;   // TODO
    
    const size_t maxLine = 1024;
    char srcBuffer[maxLine];
    char dstBuffer[maxLine];
    
    while (fgets( srcBuffer, (int)maxLine, local.fpSrc ) != nullptr)
    {
        strcpy( dstBuffer, srcBuffer );  // default dst is same as src
        
        char * pSrc = srcBuffer;
        char * pDst = dstBuffer;

        while (*pSrc && isblank(*pSrc)) pSrc++;                 // skip white-space

        if (*pSrc == '#') goto PUT_LINE;                        // keep comment lines
        
        if (strncmp( pSrc, "Block", 5 ) == 0) goto PUT_LINE;    // keep Block lines

        if (strncmp( pSrc, "DECAY", 5 ) == 0)                   // keep DECAY prefix on DECAY lines
        {
            pSrc += 5;
            pDst += 5;
        }
        
        {
            int    id           = 0;
            double value        = 0.0;
            char   name[121]    = {};       // 1 extra for null character
        
            if (sscanf( pSrc, "%i %lf # %120s", static_cast<int *>(&id), static_cast<double *>(&value), static_cast<char *>(name) ) == 3)
            {
                for (size_t i = 0; i < parameters.size(); ++i)
                {
                    if (parameters[i].name == name)
                    {
                        // replace value
                        value = paramValues[i];
                        sprintf( pDst, "  %i  %.13E  # %s\n", FMT_I(id), FMT_F(value), FMT_HS(name) );
                        break;
                    }
                }
            }
        }

        PUT_LINE :
            fputs( dstBuffer, local.fpDst );
    }
    
    if (!feof(local.fpSrc))
    {
        if (ferror(local.fpSrc))
        {
            throw -6;
        }
        else
        {
            throw -7;
        }
    }
}
