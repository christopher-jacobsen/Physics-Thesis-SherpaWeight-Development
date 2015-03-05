////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeightProgram.h"
#include "SherpaWeight.h"
#include "SherpaRootEvent.h"

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Sherpa and Root include files

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"

// Sherpa includes
#include <ATOOLS/Org/Exception.H>

// Root includes
#include <TFile.h>
#include <TTree.h>

#pragma clang diagnostic pop

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeightProgram::SherpaWeightProgram()
    : m_upSherpaWeight( new SherpaWeight )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeightProgram::~SherpaWeightProgram() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SherpaWeightProgram::ParseCommandLine( int argc, const char * argv[], RunParameters & param )
{
    param = RunParameters();  // clear all values in case of error

    if (argc < 3)
        goto USAGE;
    
    param.argv.push_back( argv[0] );

    param.inputRootFileName  = argv[1];
    param.outputRootFileName = argv[2];
 
    if (param.inputRootFileName.empty() || param.outputRootFileName.empty())
        goto USAGE;
    
    if (param.inputRootFileName == param.outputRootFileName)
    {
        LogMsgError( "Output root file cannot be the same as the input root file." );
        return -1;
    }
    
    for (int a = 3; a < argc; ++a)
        param.argv.push_back( argv[a] );

    return 0;

 USAGE:
    LogMsgInfo("Usage: SherpaWeight input_root_file output_root_file <sherpa_arguments ...>");
    return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SherpaWeightProgram::Run( const RunParameters & param )
{
    try
    {
        time_t timeStartRun = time(nullptr);

        // initialize sherpa weight
        m_upSherpaWeight->Initialize( param.inputRootFileName, param.argv );

        if (!m_upSherpaWeight->NCoefficients())
        {
            LogMsgInfo( "No reweighting parameters defined. Nothing to do." );
            return EXIT_SUCCESS;
        }
        
        // evaluate the events
        m_upSherpaWeight->EvaluateEvents();
        
        SaveCoefficients( param );

        time_t timeStopRun = time(nullptr);

        LogMsgInfo( "\nDone. (%u seconds)", FMT_U(timeStopRun - timeStartRun) );

        return EXIT_SUCCESS;
    }
    catch (const ATOOLS::Exception & error)
    {
        LogMsgError( "Sherpa Exception: %hs\n\t[Source %hs::%hs]",
                     FMT_HS(error.Info().c_str()), FMT_HS(error.Class().c_str()), FMT_HS(error.Method().c_str()) );
    }

    return EXIT_FAILURE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeightProgram::SaveCoefficients( const RunParameters & param )
{
    LogMsgInfo( "\n+----------------------------------------------------------+"   );
    LogMsgInfo(   "|  Calculating and Saving Coefficients                     |"   );
    LogMsgInfo(   "+----------------------------------------------------------+\n" );

    // open input file

    LogMsgInfo( "Input file : %hs", FMT_HS(param.inputRootFileName.c_str()) );
    std::unique_ptr<TFile> upInputFile( new TFile( param.inputRootFileName.c_str() ) );
    if (upInputFile->IsZombie() || !upInputFile->IsOpen())      // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to open input file (%hs).", FMT_HS(param.inputRootFileName.c_str()) );
        ThrowError( std::invalid_argument( param.inputRootFileName ) );
    }

    // create output file

    LogMsgInfo( "Output file: %hs", FMT_HS(param.outputRootFileName.c_str()) );
    std::unique_ptr<TFile> upOutputFile( new TFile( param.outputRootFileName.c_str(), "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(param.outputRootFileName.c_str()) );
        ThrowError( std::invalid_argument( param.outputRootFileName ) );
    }

    // get and setup input tree

    TTree * pInputTree = nullptr;
    upInputFile->GetObject( "t3", pInputTree );
    if (!pInputTree)
        ThrowError( "Failed to load input tree." );

    // create output tree

    TTree * pOutputTree( new TTree( "t3", "SherpaWeight" ) );   // owned by current directory
    if (pOutputTree->IsZombie())
        ThrowError( "Failed to construct output tree." );

    pOutputTree->SetDirectory( upOutputFile.get() );   // attach to output file, output file now owns tree and will call delete
    pOutputTree->SetAutoSave(0);                       // disable autosave
    
    // create event and connect/declare variables

    SherpaRootEvent currentEvent;

    currentEvent.SetInputTree(  pInputTree  );
    currentEvent.SetOutputTree( pOutputTree );
    
    // add coefficient output variables
    
    const SherpaWeight::StringVector & coefNames = m_upSherpaWeight->CoefficientNames();

    SherpaWeight::DoubleVector eventCoefs( m_upSherpaWeight->NCoefficients() );  // storage for event reweight coefficients

    if (eventCoefs.size())
    {
        if (coefNames.size() != eventCoefs.size())
        {
            ThrowError( std::logic_error( "Number of coefficient names (" + std::to_string(coefNames.size()) +
                                          ") does not match number of coefficients (" + std::to_string(eventCoefs.size())  + ")." ) );
        }

        // add entire coefficient vector as a branch
        pOutputTree->Branch( "Fij", &eventCoefs );

        // add each coefficient as an individual branch
        size_t index = 0;
        for (const std::string & name : coefNames)
        {
            size_t pos = name.find( "_", name.find( "_", name.find("_") + 1 ) + 1 );  // find 3rd underscore
            std::string shortName = name.substr(0, pos );
            TBranch * pBranch = pOutputTree->Branch( shortName.c_str(), &eventCoefs[index] );
            pBranch->SetTitle( (name + "/D").c_str() );
            ++index;
        }
    }

    // loop through entries
    
    const Long64_t nEntries     = pInputTree->GetEntries();
    const Long64_t logFrequency = std::max( nEntries / 10, Long64_t(1) );

    LogMsgInfo( "\nCalculating coefficients for %lli events ...", FMT_LLI(nEntries) );

    time_t timeStartProcess = time(nullptr);

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        if (pInputTree->LoadTree(iEntry) < 0)
            ThrowError( "LoadTree failed on entry " + std::to_string(iEntry+1) );
        
        if (pInputTree->GetEntry(iEntry) < 0)
            ThrowError( "GetEntry failed on entry " + std::to_string(iEntry+1) );
        
        {
            SherpaWeight::DoubleVector coefs = m_upSherpaWeight->CoefficientValues( currentEvent.id );
            if (coefs.empty())
            {
                LogMsgWarning( "No coefficients for event %%li (id %i).", FMT_LLI(iEntry+1), FMT_I(currentEvent.id) );
            }
            else if (iEntry % logFrequency == 0)
            {
                LogMsgInfo( "Event %lli (id %i):", FMT_LLI(iEntry+1), FMT_I(currentEvent.id) );

                size_t index = 0;
                for (double value : coefs)
                {
                    const char * pName = (index < coefNames.size()) ? coefNames[index].c_str() : "Unknown";
                    LogMsgInfo( "%-20hs:\t%.15E", FMT_HS(pName), FMT_F(value) );
                    ++index;
                }
                
                LogMsgInfo( "" );
            }
            
            if (coefs.size() == eventCoefs.size())
                std::copy( coefs.begin(), coefs.end(), eventCoefs.begin() );
            else
                std::fill( eventCoefs.begin(), eventCoefs.end(), 0.0 );
        }
        
        
        if (pOutputTree->Fill() < 0)
            ThrowError( "Fill failed on entry " + std::to_string(iEntry+1) );
    }
    
    upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();

    time_t timeStopProcess = time(nullptr);

    LogMsgInfo( "%lli events completed. (%u seconds)", FMT_LLI(nEntries), FMT_U(timeStopProcess - timeStartProcess) );
}
