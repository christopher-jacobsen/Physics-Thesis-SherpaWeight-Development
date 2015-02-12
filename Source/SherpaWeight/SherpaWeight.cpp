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
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::~SherpaWeight() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::Initialize( const std::vector<const char *> & argv )
{
    m_argv = argv;
    m_argv.push_back( "INIT_ONLY=2" ); // prevent Sherpa from starting the cross section integration

    //////

#if 0
    {
        SHERPA::Sherpa bob;
        
        if (!bob.InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
        {
            LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
            throw int(-2);  // TODO
        }

        if (!bob.InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
        {
            LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
            throw int(-2);  // TODO
        }
    }
#endif

    SHERPA::Sherpa sherpa;
    
    if (!sherpa.InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
    {
        LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        throw int(-2);  // TODO
    }
    
    // get the various configuration paths and file names
    
    SHERPA::Initialization_Handler * pInitHandler = sherpa.GetInitHandler();
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
    
    Int_t nParam = (Int_t)m_parameters.size();
    if (m_parameters.size())
    {
        TMatrixD evalMatrix, invMatrix;
        GetBilinearMatrices( m_parameters, evalMatrix, invMatrix );
        
        if (evalMatrix.GetNcols() != nParam)
            throw -3; // TODO
        
        m_nEvaluations = evalMatrix.GetNrows();
        
        // create list of param_card files
        
        std::string srcFilePath = runPath + paramFile;
        std::string dstPath     = runPath + "SherpaWeight/";

        // create destination param_card directory
        {
            std::string command = "mkdir \"" + dstPath + "\"";
            system( command.c_str() );
        }

        for (Int_t i = 0; i < m_nEvaluations; ++i)
        {
            // copy param values from TMatrixTRow<Double_t> to std::vector<double>
            std::vector<double> paramValues( nParam );
            std::copy_n( evalMatrix[i].GetPtr(), nParam, paramValues.begin() );

            // construct destination file path
            char dstFile[40];
            sprintf( dstFile, "param_card_%03i.dat", FMT_I(i) );  // %i is max 10 chars
            
            std::string dstFilePath( dstPath + dstFile );

            // create param_card with new parameter values
            CreateParamCard( srcFilePath, dstFilePath, m_parameters, paramValues );

            m_paramCards.push_back( dstFilePath );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::ProcessEvent( int32_t eventId, size_t nInParticles, const std::vector<int> & particleCodes, const ATOOLS::Vec4D_Vector & particleMomenta )
{
    if (!nInParticles || (nInParticles > particleCodes.size()))
    {
        throw int(-2);  // TODO
    }
    
    if (particleCodes.size() != particleMomenta.size())
    {
        throw int(-2);  // TODO
    }
    
    if (m_paramCards.size() != m_nEvaluations)
    {
        throw int(-2);  // TODO
    }
    
    for (size_t run = 0; run < m_nEvaluations; ++run)
    {
        //std::string                 argParamCard( "FR_PARAMCARD=" + m_paramCards[run] );
        //std::vector<const char *>   runArgv(m_argv);
        
        //runArgv.push_back(argParamCard.c_str());

        SHERPA::Sherpa sherpa;
        
        if (!sherpa.InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
        {
            LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
            throw int(-2);  // TODO
        }
        
        SherpaMECalculator meCalc( &sherpa );

        for (size_t i = 0; i < particleCodes.size(); ++i)
        {
            if (i < nInParticles)
                meCalc.AddInFlav( particleCodes[i] );
            else
                meCalc.AddOutFlav( particleCodes[i] );
        }

        try
        {
            meCalc.Initialize();
        }
        catch (const ATOOLS::Exception & error)
        {
            LogMsgError( "Sherpa exception caught: \"%hs\" in %hs::%hs",
                FMT_HS(error.Info().c_str()), FMT_HS(error.Class().c_str()), FMT_HS(error.Method().c_str()) );
            
            LogMsgInfo( "Sherpa process name: \"%hs\"", FMT_HS(meCalc.Name().c_str()) );
            
            const ATOOLS::Cluster_Amplitude * pAmp = meCalc.GetAmp();
            if (!pAmp)
                LogMsgInfo("Sherpa process cluster amplitude: null");
            else
            {
                //const ATOOLS::ClusterLeg_Vector & legs = pAmp->Legs();
                //int bob = 5;
            }

            throw;
        }

        meCalc.SetMomenta( particleMomenta );
        
        double me = meCalc.MatrixElement();
        LogMsgInfo( "Event %i: ME=%E", FMT_I(eventId), FMT_F(me) );

    #if 0
        {
            MODEL::ScalarConstantsMap * pConstants = meCalc.GetModelScalarConstants();
            if (pConstants)
            {
                /**/
                for (const auto & entry : *pConstants)
                {
                    std::cout << entry.first << ": " << entry.second << std::endl;
                }
                /**/
                
                double aEWM1 = 227.9;
                
                MODEL::ScalarConstantsMap::iterator itrFind = pConstants->find( "aEWM1" );
                if (itrFind != pConstants->end())
                {
                    MODEL::ScalarConstantsMap::value_type & constantPair = *itrFind;
                    constantPair.second = aEWM1;
                
                    double sw = pConstants->at("sw");
                    double cw = pConstants->at("cw");

                    double aEW = pow(aEWM1,-1.);
                    double ee = 2.*sqrt(aEW)*sqrt(M_PI);
                    double gw = ee*pow(sw,-1.);
                    double g1 = ee*pow(cw,-1.);

                    pConstants->at("aEW") = aEW;
                    pConstants->at("ee")  = ee;
                    pConstants->at("gw")  = gw;
                    pConstants->at("g1")  = g1;
                    
                }

                /**/
                std::cout << "------------" << std::endl;
                for (const auto & entry : *pConstants)
                {
                    std::cout << entry.first << ": " << entry.second << std::endl;
                }
                /**/
            }
        }
        
        double me2 = meCalc.MatrixElement();
        LogMsgInfo( "Event %i: ME2=%E", FMT_I(eventId), FMT_F(me2) );
    #endif
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::GetBilinearMatrices( const std::vector<ReweightParameter> & parameters, TMatrixD & evalMatrix, TMatrixD & invCoefMatrix )  // static
{
    // clear outputs
    evalMatrix    .Clear();
    invCoefMatrix .Clear();
    
    Int_t nParam = (Int_t)parameters.size();
    if (nParam < 1)
        return;

    Int_t nCoefs = (nParam + 1) * (nParam + 2) / 2;  // nCoefs >= 3 [nCoefs = 3,6,10,15,...]

    // evalMatrix  rows: nCoefs   columns: nParam
    // invmatrix   rows: nCoefs   columns: nCoefs

    //// fill in evalMatrix
    {
        evalMatrix.ResizeTo( nCoefs, nParam );     // new elements are zeroed

        // zero all elements (just in case)
        evalMatrix = 0.0;
        
        // evalMatrix[0] = [0,...,0]
        
        // evalMatrix[1..nParam] = diagonal of 1's
        for (Int_t i = 0; i < nParam; ++i)
        {
            evalMatrix[i+1][i] = 1.0;
        }

        // evalMatrix[nParam+1..2*nParam] = diagonal of -1's
        for (Int_t i = 0; i < nParam; ++i)
        {
            evalMatrix[i+1+nParam][i] = -1.0;
        }
        
        // remainder with rotations of [1,-1,0,...], [1,0,-1,0,...], ...
        for (Int_t r = 2 * nParam + 1, g = 1; r < nCoefs; )
        {
            Double_t row[nParam];

            memset( row, 0, sizeof(row) );
            
            row[0]   =  1.0;
            row[g++] = -1.0;
            
            for (Int_t i = 0; (i < nParam) && (r < nCoefs); ++i)
            {
                evalMatrix[r++] = TVectorD( nParam, row );

                std::rotate( row, row + nParam - 1, row + nParam );
            }
        }
        
        // multiply with scale and add offset
        
        for (Int_t r = 0; r < nCoefs; ++r)
        {
            for (Int_t c = 0; c < nParam; ++c)
            {
                const ReweightParameter & param = parameters[c];
                Double_t &                entry = evalMatrix[r][c];

                entry = entry * param.scale + param.offset;
            }
        }
    }

    {
        Double_t Mvector[nCoefs][nCoefs];

        // filling of Mvector
        Int_t k;
        for (Int_t r = 0; r < nCoefs; ++r)
        {
            Mvector[r][0] = 1.0;
            k = 1 + nParam;

            for (Int_t c = 0; c < nParam; ++c)
            {
                Double_t value   = evalMatrix[r][c];
                Double_t sqValue = value * value;

                Mvector[r][c + 1] = value;          // single terms
                Mvector[r][k]     = sqValue;        // square terms
                
                k += nParam - c;
            }

            // interference terms
            k = 1 + nParam;
            for (Int_t c = 0; c < nParam - 1; ++c)
            {
                for (Int_t j = c + 1; j < nParam; ++j)
                {
                    Double_t value = evalMatrix[r][c] * evalMatrix[r][j];
                
                    Mvector[r][k + j - c] = value;
                }

                k += nParam - c;
            }
        }
        
        invCoefMatrix.ResizeTo( nCoefs, nCoefs  );     // new elements are zeroed
        invCoefMatrix.SetMatrixArray( Mvector[0] );
        invCoefMatrix.Invert();
        
        // get rid of numerical precision residuals
        for (Int_t i = 0; i < nCoefs; ++i)
        {
            for (Int_t j = 0; j < nCoefs; ++j)
            {
                if (abs(invCoefMatrix[i][j]) < 1.0E-10)
                    invCoefMatrix[i][j] = 0.0;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::CreateParamCard( const std::string & srcFilePath, const std::string & dstFilePath,
                                    const std::vector<ReweightParameter> & parameters,
                                    const std::vector<double> & paramValues )    // static
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
