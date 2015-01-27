////////////////////////////////////////////////////////////////////////////////////////////////////
//  SMEX.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SMEX.h"

#include "SherpaEvent.h"

#include "common.h"

// TODO: pragma warnings away?
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SMEX
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SMEX::SMEX()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SMEX::~SMEX() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SMEX::ParseCommandLine( int argc, const char * argv[], RunParameters & param )
{
    param = RunParameters();  // clear all values in case of error

    if (argc != 3)
        goto USAGE;
    
    param.inputRootFileName  = argv[1];
    param.outputRootFileName = argv[2];
 
    if (param.inputRootFileName.empty() || param.outputRootFileName.empty())
        goto USAGE;
    
    if (param.inputRootFileName == param.outputRootFileName)
    {
        LogMsgError( "Output root file cannot be the same as the input root file." );
        return -1;
    }
    
    return 0;

 USAGE:
    LogMsgInfo("Usage: SMEX input_root_file output_root_file");
    return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SMEX::Run( const RunParameters & param )
{
    /*
    struct Local  // local object for automated cleanup
    {
        ~Local()
        {
        }
    }
    local;
    */
    
    //gSystem->Load("libTree");
    //gROOT->ProcessLine( "#include <vector>" );

    // open input file
    LogMsgInfo( "Input file : %hs", FMT_HS(param.inputRootFileName.c_str()) );
    std::unique_ptr<TFile> upInputFile( new TFile( param.inputRootFileName.c_str() ) );
    if (upInputFile->IsZombie() || !upInputFile->IsOpen())      // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to open input file (%hs).", FMT_HS(param.inputRootFileName.c_str()) );
        return -2;
    }

    // create output file
    LogMsgInfo( "Output file: %hs", FMT_HS(param.outputRootFileName.c_str()) );
    std::unique_ptr<TFile> upOutputFile( new TFile( param.outputRootFileName.c_str(), "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(param.outputRootFileName.c_str()) );
        return -2;
    }

    // get and setup input tree
    TTree * pInputTree = nullptr;
    upInputFile->GetObject( "t3", pInputTree );
    if (!pInputTree)
    {
        LogMsgError( "Failed to load input file tree." );
        return -2;
    }

    // create output tree
    TTree * pOutputTree( new TTree( "t3", "SMEX" ) );   // owned by current directory
    if (pOutputTree->IsZombie())
    {
        LogMsgError("Failed to construct output tree.");
        return -2;
    }
    pOutputTree->SetDirectory( upOutputFile.get() );   // attach to output file, output file now owns tree and will call delete
    pOutputTree->SetAutoSave(0);                       // disable autosave

    // create event and connect/declare variables
    SherpaEvent currentEvent;
    currentEvent.SetInputTree( pInputTree );
    currentEvent.SetOutputTree( pOutputTree );
    
    // loop through entries
    
    const Long64_t nEntries = pInputTree->GetEntries();
    
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
    {
        if (pInputTree->LoadTree(iEntry) < 0)
        {
            LogMsgError( "LoadTree failed on entry %lli", FMT_LLI(iEntry) );
            break;
        }
        
        if (pInputTree->GetEntry(iEntry) < 0)
        {
            LogMsgError( "GetEntry failed on entry %lli", FMT_LLI(iEntry) );
            break;
        }
        
        ProcessEvent( currentEvent );
        
        if (pOutputTree->Fill() < 0)
        {
            LogMsgError( "Fill failed on entry %lli", FMT_LLI(iEntry) );
            break;
        }
    }
    
    upOutputFile->Write( 0, TFile::kOverwrite );
    upOutputFile->Close();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SMEX::ProcessEvent( SherpaEvent & event )
{
    


}
