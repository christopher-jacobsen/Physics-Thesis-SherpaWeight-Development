////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaWeight.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeight.h"
#include "MERootEvent.h"

#include "common.h"

#include <limits>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Sherpa and Root include files

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"

// Sherpa includes
#include <SHERPA/Main/Sherpa.H>
#include <SHERPA/Initialization/Initialization_Handler.H>
#include <MODEL/Main/Model_Base.h>
#include <ATOOLS/Org/Run_Parameter.H>
#include <ATOOLS/Org/Data_Reader.H>

// Root includes
#include <TFile.h>
#include <TTree.h>
#include <TMatrixD.h>

#pragma clang diagnostic pop

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
    if (m_bOwnModel)
        delete m_pModel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::Initialize( const std::string & eventFileName, const std::vector<const char *> & argv,
                               bool bReadParameters /*= true*/ )
{
    // TODO: check if already initialized
    
    m_eventFileName = eventFileName;
    m_argv          = argv;
    
    {
        std::string application( argv[0] );
        m_appRunPath = application.substr( 0, application.rfind("/") ) + "/";
    }
    
    // initialize sherpa
    {
        //std::vector<const char *> runArgv(m_argv);
        //runArgv.push_back( "INIT_ONLY=2" ); // prevent Sherpa from starting the cross section integration
        
        if (!m_upSherpa->InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
            ThrowError( "Failed to initialize Sherpa framework. Check Run.dat file." );
    }
    
    SHERPA::Initialization_Handler * pInitHandler = m_upSherpa->GetInitHandler();
    if (!pInitHandler)
        ThrowError( "Failed to get Sherpa initialization handler. Is Sherpa initialized?" );

    m_sherpaRunPath = pInitHandler->Path();

    // determine the model interface

    if (!m_pModel)
    {
        MODEL::Model_Base * pModel = pInitHandler->GetModel();
        if (!pModel)
            ThrowError( "Model is undefined. Set MODEL in (model) section of Run.dat file." );

        std::string modelName = pModel->Name();

        if (modelName == "FeynRules")
            m_pModel = new FeynRulesModel;

        if (!m_pModel)
            ThrowError( "Model " + modelName + " is not supported." );
    }

    // read the run file/section for run parameters
    {
        std::string runFileSection  = pInitHandler->File();
        std::string runFileBase     = runFileSection.substr( 0, runFileSection.find("|") );  // strip off section declaration following '|'

        ATOOLS::Data_Reader reader(" ",";","!","=");
        reader.AddComment("#");
        reader.AddWordSeparator("\t");
        reader.SetInputPath( SherpaRunPath() );
        reader.SetInputFile( runFileSection );
        
        // Note: By default, Data_Reader appends the (run) section and command line parameters to the section of the file it is reading.

        m_sherpaWeightFileSection = reader.GetValue<std::string>( "SHERPA_WEIGHT_FILE", runFileBase + "|(SherpaWeight){|}(SherpaWeight)" );
    }

    // determine and create temporary work directory
    {
        m_tmpPath = SherpaRunPath() + "SherpaWeight.tmp/";
        std::string command = "mkdir \"" + m_tmpPath + "\"";
        system( command.c_str() );
    }

    // read the model file/section and initialize the model interface
    {
        std::string modelFile = ATOOLS::rpa->gen.Variable("MODEL_DATA_FILE");
        if (modelFile.empty())
            ThrowError( "No model file/section defined. Set MODEL_DATA_FILE." );

        ATOOLS::Data_Reader reader(" ",";","!","=");
        reader.AddComment("#");
        reader.AddWordSeparator("\t");
        reader.SetInputPath( SherpaRunPath() );
        reader.SetInputFile( modelFile );

        m_pModel->Initialize( *this, reader );
    }

    // optionally read reweight parameters from file
    
    if (bReadParameters)
        ReadParametersFromFile();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::ReadParametersFromFile( const char * filePath /*= nullptr*/ )
{
    ParameterVector params;
    
    ATOOLS::Data_Reader reader(" ",";","!","=");
    reader.AddComment("#");
    reader.AddWordSeparator("\t");
    reader.SetInputPath( SherpaRunPath() );
    reader.SetInputFile( m_sherpaWeightFileSection );
    
    // Note: By default, Data_Reader appends the (run) section and command line parameters to the section of the file it is reading.
    reader.SetAddCommandLine(false);  // read the raw file input
    
    std::vector< std::vector<std::string> > stringMatrix;
    if (reader.MatrixFromFile( stringMatrix ))
    {
        params.reserve( stringMatrix.size() );    // reserve space for all entries
        
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
                    
                    // TODO: try using std::stof
                    
                    std::stringstream stream;
                    stream.precision(12);
                    
                    stream << row[1];
                    stream >> param.scale;
                    if (stream.fail() || !stream.eof())
                        ThrowError( "Could not convert string to floating point value." );  // caught below
                    
                    if (row.size() > 2)
                    {
                        stream.clear();
                        stream << row[2];
                        stream >> param.offset;
                        if (stream.fail() || !stream.eof())
                            ThrowError( "Could not convert string to floating point value." );  // caught below
                    }
                }
                catch (const std::exception & error)
                {
                    ThrowError( "Failed to read scale/offset for reweight parameter " + param.name );
                }
            }
            
            params.push_back( std::move(param) );
        }
    }
    
    SetParameters( params );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::SetParameters( const ParameterVector & params )
{
    // ensure all names are unique
    {
        std::set<std::string> duplicateNames;

        for (auto itr1 = params.cbegin(), end = params.cend(); itr1 != end; ++itr1)
        {
            for (auto itr2 = itr1 + 1; itr2 != end; ++itr2)
            {
                if (itr1->name == itr2->name)
                {
                    duplicateNames.insert( itr1->name );
                    break;
                }
            }
        }

        if (!duplicateNames.empty())
        {
            std::string strDup;
            for (const std::string & s : duplicateNames)
            {
                if (!strDup.empty()) strDup += ", ";
                strDup += s;
            }
            
            ThrowError( "Duplicate reweight parameters: " + strDup );
        }
    }

    // validate parameters with model
    m_pModel->ValidateParameters(params);

    // set the parameters
    m_parameters = params;
    
    // calculate bilinear matrices
    
    GetBilinearMatrices( m_parameters, m_evalMatrix, m_invCoefMatrix, m_coefNames );

    m_matrixElements.clear();  // cleanup any previous run
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::EvaluateEvents()
{
    size_t nEvaluations = NEvaluations();
    
    m_matrixElements.clear();  // cleanup any previous run

    if (nEvaluations == 0)
        return;

    // run evaluations
    
    std::string outputFile( TemporaryPath() + "SherpaME-tmp.root" );
    std::string baseCommand( ApplicationRunPath() + "SherpaME" );
    
    baseCommand += " " + m_eventFileName;
    baseCommand += " " + outputFile;

    for (size_t i = 1; i < m_argv.size(); ++i)
        baseCommand += std::string(" ") + m_argv[i];
    
    for (size_t run = 0; run < nEvaluations; ++run)
    {
        std::string modelArgs = m_pModel->CommandLineArgs( m_parameters, m_evalMatrix[run] );
        std::string command   = baseCommand + " " + modelArgs;

        LogMsgInfo( "" );
        LogMsgInfo( "+----------------------------------------------------------+");
        LogMsgInfo( "|  Evaluation Run %u                                       |", FMT_U(run + 1) );
        LogMsgInfo( "+----------------------------------------------------------+");
        LogMsgInfo( "" );

        LogMsgInfo( "Running command:" );
        LogMsgInfo( "%hs\n", FMT_HS(command.c_str()) );

        int result = system( command.c_str() );
        if (result != 0)
            ThrowError( "Evaluation failed." );
        
        AddMatrixElementsFromFile( outputFile.c_str() );
    }
    
    // validate matrix elements
    {
        EventMatrixElementMap::iterator itr = m_matrixElements.begin();
        EventMatrixElementMap::iterator end = m_matrixElements.end();
        while (itr != end)
        {
            if (itr->second.size() != nEvaluations)
            {
                LogMsgWarning( "Discarding event %i. Evaluations: %u, require: %u.", FMT_I(itr->first), FMT_U(itr->second.size()), FMT_U(nEvaluations) );
                itr = m_matrixElements.erase(itr);  // returned itr is next itr in map
                continue;
            }
            ++itr;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
const SherpaWeight::DoubleVector & SherpaWeight::MatrixElements( int32_t eventId ) const
{
    EventMatrixElementMap::const_iterator itrFind = m_matrixElements.find(eventId);

    if (itrFind == m_matrixElements.end())
    {
        static const DoubleVector empty;
        return empty;
    }
    
    return itrFind->second;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
const SherpaWeight::DoubleVector SherpaWeight::CoefficientValues( int32_t eventId ) const
{
    DoubleVector coefs;
    
    const DoubleVector & matrixElements = MatrixElements(eventId);

    if (matrixElements.size() == m_invCoefMatrix.size())
    {
        size_t nCoefs = m_invCoefMatrix.size();
        
        coefs.resize(nCoefs);
        
        for (size_t i = 0; i < nCoefs; ++i)
        {
            double & coef = coefs[i];

            for (size_t j = 0; j < nCoefs; ++j)
            {
                coef += m_invCoefMatrix[i][j] * matrixElements[j];
            }
        }
    }

    return coefs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::AddMatrixElementsFromFile( const char * filePath )
{
    // open input file
    
    LogMsgInfo( "Adding ME from root file: %hs", FMT_HS(filePath) );

    std::unique_ptr<TFile> upInputFile( new TFile(filePath) );

    if (upInputFile->IsZombie() || !upInputFile->IsOpen())      // IsZombie is true if constructor failed
        ThrowError( "Failed to open ME root file (" + std::string(filePath) + ")" );
    
    // get and setup input tree
    
    TTree * pInputTree = nullptr;
    upInputFile->GetObject( "SherpaME", pInputTree );
    if (!pInputTree)
        ThrowError( "Failed to load input tree." );
    
    MERootEvent inputEvent;
    inputEvent.SetInputTree( pInputTree );
    
    // loop through and process each input event
    
    const Long64_t nEntries = pInputTree->GetEntries();
    
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        if (pInputTree->LoadTree(iEntry) < 0)
            ThrowError( "LoadTree failed on entry " + std::to_string(iEntry) );
        
        if (pInputTree->GetEntry(iEntry) < 0)
            ThrowError( "GetEntry failed on entry " + std::to_string(iEntry) );

        AddMatrixElement( inputEvent.id, inputEvent.me );
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::AddMatrixElement( int32_t eventId, double me )
{
    DoubleVector & vector = m_matrixElements[eventId];
    vector.push_back( me );
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
        coefNames.resize( nCoefs );

        for (size_t r = 0; r < nCoefs; ++r)
        {
            size_t c = 0;
            
            // 0th order term
            {
                coef[r][c] = 1.0;
                if (r == 0) coefNames[c] = "F00";
                ++c;
            }

            // 1st order terms
            for (size_t i = 0; i < nParam; ++i)
            {
                coef[r][c] = eval[r][i];
                if (r == 0) coefNames[c] = "F0" + std::to_string(i+1) + "_" + parameters[i].name;
                ++c;
            }

            // square and cross terms
            for (size_t i = 0; i < nParam; ++i)
            {
                for (size_t j = i; j < nParam; ++j)
                {
                    coef[r][c] = eval[r][i] * eval[r][j];
                    if (r == 0)
                    {
                        coefNames[c] = "F" + std::to_string(i+1) + std::to_string(j+1) + "_" + parameters[i].name;
                        if (i != j) coefNames[c] += "_" + parameters[j].name;
                    }
                    ++c;
                }
            }
        }
    }
    
    DoubleMatrix coefMatrix;  // for debugging
    {
        coefMatrix.resize( nCoefs );
        for (size_t r = 0; r < nCoefs; ++r)
        {
            coefMatrix[r].assign( coef[r], coef[r] + nCoefs );
        }
    }
    
    // invert coefficient matrix
    {
        TMatrixD matrix( (Int_t)nCoefs, (Int_t)nCoefs, static_cast<const Double_t *>( coef[0] ) );  // copy coef to TMatrixD

        Double_t determinant = 0;
        matrix.Invert( &determinant );
        if (fabs(determinant) < std::numeric_limits<Double_t>::epsilon())
            ThrowError( std::logic_error( "Bilinear coefficient matrix is singular and cannot be inverted." ) );

        matrix.GetMatrix2Array( static_cast<double *>( coef[0] ) );     // copy back to coef
    }

    /*
    // get rid of numerical precision residuals
    {
        const double epsilon = std::numeric_limits<float>::epsilon();
        
        for (size_t r = 0; r < nCoefs; ++r)
        {
            for (size_t c = 0; c < nCoefs; ++c)
            {
                double & entry = coef[r][c];

                if (fabs(entry) < epsilon)
                    entry = 0.0;
            }
        }
    }
    */

    // copy to output matrices
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
// class SherpaWeight::FeynRulesModel
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::FeynRulesModel::FeynRulesModel()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::FeynRulesModel::~FeynRulesModel() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::FeynRulesModel::Initialize( SherpaWeight & parent, ATOOLS::Data_Reader & modelFileSectionReader )
{
    LogMsgInfo( "------------------------------------------------------------" );
    LogMsgInfo( "  FeynRules Model Interface" );
    LogMsgInfo( "------------------------------------------------------------" );

    m_pParent = &parent;

    std::string paramCard = modelFileSectionReader.GetValue<std::string>( "FR_PARAMCARD", std::string("param_card.dat") );
    if (paramCard.find("|") != std::string::npos)
        ThrowError( "FR_PARAMCARD must define only a file and not a section within a file." );

    m_sourceParamCardFile = paramCard;

    LogMsgInfo( ("  FR_PARAMCARD=" + m_sourceParamCardFile).c_str() );
    LogMsgInfo( "------------------------------------------------------------" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::FeynRulesModel::ValidateParameters( const ParameterVector & /*params*/ )
{
    // do nothing: validate during creation of FeynRules param card
}

////////////////////////////////////////////////////////////////////////////////////////////////////
std::string SherpaWeight::FeynRulesModel::CommandLineArgs( const ParameterVector & params, const DoubleVector & paramValues )
{
    std::string srcFilePath = m_pParent->SherpaRunPath() + m_sourceParamCardFile;
    std::string dstFilePath = m_pParent->TemporaryPath() + "param_card_me.dat";
        
    // create param_card with new parameter values
    CreateFeynRulesParamCard( srcFilePath, dstFilePath, params, paramValues );

    std::string cmdArgs = "\"FR_PARAMCARD=" + dstFilePath + "\"";
    return cmdArgs;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::FeynRulesModel::CreateFeynRulesParamCard( const std::string & srcFilePath, const std::string & dstFilePath,
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
        ThrowError( std::invalid_argument( "CreateFeynRulesParamCard: mismatch in size of parameter and value vectors." ) );

    // setup remainingNames to track outstanding parameter substitutions
    std::set<std::string> remainingNames;

    for (const ReweightParameter & p : parameters)
        remainingNames.insert( p.name );

    // open files

    local.fpSrc = fopen( srcFilePath.c_str(), "rt" );
    if (!local.fpSrc)
        ThrowError( std::system_error( errno, std::generic_category(), "Failed to open FeynRules param card (" + srcFilePath + ")" ) );

    local.fpDst = fopen( dstFilePath.c_str(), "wt" );
    if (!local.fpDst)
        ThrowError( std::system_error( errno, std::generic_category(), "Failed to create FeynRules param card (" + srcFilePath + ")" ) );

    // copy all lines from source to destination, substituting parameter values

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
                        
                        remainingNames.erase( parameters[i].name );
                        
                        break;
                    }
                }
            }
        }

        PUT_LINE :
            fputs( dstBuffer, local.fpDst );
    }
    
    if (!feof(local.fpSrc))
        ThrowError( "Failed to read FeynRules param card (" + srcFilePath + ")" );

    if (!remainingNames.empty())
    {
        std::string strList;
        for (const std::string & s : remainingNames)
        {
            if (!strList.empty()) strList += ", ";
            strList += s;
        }
        
        ThrowError( "FeynRules param card (" + srcFilePath + ") missing reweight parameters: " + strList );
    }
}
