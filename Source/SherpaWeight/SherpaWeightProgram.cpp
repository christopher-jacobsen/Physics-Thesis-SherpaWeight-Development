////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeightProgram.h"
#include "SherpaWeight.h"
#include "SherpaRootEvent.h"

#include "common.h"

// TODO: pragma warnings away?

// Sherpa includes
#include <ATOOLS/Org/Exception.H>

// Root includes
#include <TFile.h>
#include <TTree.h>

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
        // initialize sherpa weight
        m_upSherpaWeight->Initialize( param.inputRootFileName, param.argv );

        // evaluate the events
        m_upSherpaWeight->EvaluateEvents();
        
        SaveCoefficients( param );

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
        ThrowError( "Failed to load input file tree." );

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

    // loop through entries
    
    const Long64_t nEntries = pInputTree->GetEntries();
    
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        if (pInputTree->LoadTree(iEntry) < 0)
        {
            LogMsgError( "LoadTree failed on entry %lli", FMT_LLI(iEntry) );
            ThrowError( std::errc::io_error, "Save of coefficients is incomplete." );
        }
        
        if (pInputTree->GetEntry(iEntry) < 0)
        {
            LogMsgError( "GetEntry failed on entry %lli", FMT_LLI(iEntry) );
            ThrowError( std::errc::io_error, "Save of coefficients is incomplete." );
        }
        
        {
            SherpaWeight::DoubleVector coefs = m_upSherpaWeight->CoefficientValues( currentEvent.id );
            if (coefs.empty())
                LogMsgWarning( "No coefficients for entry %%li, event id %i", FMT_LLI(iEntry), FMT_I(currentEvent.id) );
            else
            {
                LogMsgInfo( "Event %i coefficients:", FMT_I(currentEvent.id) );

                const SherpaWeight::StringVector & names = m_upSherpaWeight->CoefficientNames();
                size_t index = 0;
                for (double value : coefs)
                {
                    const char * pName = (index < names.size()) ? names[index].c_str() : "Unknown";
                    LogMsgInfo( "%-20hs:\t%.15E", FMT_HS(pName), FMT_F(value) );
                    ++index;
                }
                
                LogMsgInfo( "" );
            }
        }
        
        
        if (pOutputTree->Fill() < 0)
        {
            LogMsgError( "Fill failed on entry %lli", FMT_LLI(iEntry) );
            ThrowError( std::errc::io_error, "Save of coefficients is incomplete." );
        }
    }
    
    upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();
}
