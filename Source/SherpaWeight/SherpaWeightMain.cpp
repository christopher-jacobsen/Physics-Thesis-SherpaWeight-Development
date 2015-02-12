////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightMain.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeightMain.h"
#include "SherpaWeight.h"
#include "SherpaRootEvent.h"

#include "common.h"

#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Math/Vector.H>

// TODO: pragma warnings away?
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeightMain::SherpaWeightMain()
    : m_upSherpaWeight( new SherpaWeight )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeightMain::~SherpaWeightMain() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SherpaWeightMain::ParseCommandLine( int argc, const char * argv[], RunParameters & param )
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
int SherpaWeightMain::Run( const RunParameters & param )
{
    try
    {
        m_upSherpaWeight->Initialize( param.argv );
        
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
        TTree * pOutputTree( new TTree( "t3", "SherpaWeight" ) );   // owned by current directory
        if (pOutputTree->IsZombie())
        {
            LogMsgError("Failed to construct output tree.");
            return -2;
        }
        pOutputTree->SetDirectory( upOutputFile.get() );   // attach to output file, output file now owns tree and will call delete
        pOutputTree->SetAutoSave(0);                       // disable autosave

        // create event and connect/declare variables
        SherpaRootEvent currentEvent;
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
    catch (const ATOOLS::Exception & error)
    {
        LogMsgError( "Sherpa exception caught: \"%hs\" in %hs::%hs",
            FMT_HS(error.Info().c_str()), FMT_HS(error.Class().c_str()), FMT_HS(error.Method().c_str()) );
        return -3;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeightMain::ProcessEvent( SherpaRootEvent & event )
{
    // validate event
    
    if (event.nparticle > SherpaRootEvent::max_nparticle )
    {
        LogMsgError( "Number of outgoing particles in event exceeds maximum. nparticles=%i (max %u).", FMT_I(event.nparticle), FMT_U(SherpaRootEvent::max_nparticle) );
        return;
    }

    if (event.nparticle <= 0)
    {
        LogMsgWarning( "Skipping event %i. No outgoing particles.", FMT_I(event.id) );
        return;
    }

    // calculate incoming 4-momentum
    ATOOLS::Vec4D P_in1, P_in2;
    {
        ATOOLS::Vec4D P_out;
        for (Int_t i = 0; i < event.nparticle; ++i)
        {
            ATOOLS::Vec4D P_part( event.E[i], event.px[i], event.py[i], event.pz[i] );
            P_out += P_part;
        }
        
        double mass = P_out.Mass();
        //double tau  = event.x1 * event.x2;          // x1 * x2 = mass^2 / s = tau
        //double s    = mass * mass / tau;
        //double P    = sqrt(s/4-1);                  // P = proton momentum, s = 4(M_p^2 + P^2)
        //double pz1  =  event.x1 * P;
        //double pz2  = -event.x2 * P;
        
        // mass^2 = (e1+e2)^2-(pz1+pz2)^2 = (x1+x2)^2 * P^2 - x^2 * P^2

        // assume massless incoming partons
        
        double E1 = (P_out[0] + P_out[3]) / 2;
        double E2 = (P_out[0] - P_out[3]) / 2;

        P_in1[0] = E1; //pz1;   // fabs(pz1)
        P_in1[3] = E1; //pz1;

        P_in2[0] =  E2; //-pz2;  // fabs(pz2)
        P_in2[3] = -E2; // pz2;
        /*
        P_in1[0] = pz1;   // fabs(pz1)
        P_in1[3] = pz1;

        P_in2[0] = -pz2;  // fabs(pz2)
        P_in2[3] =  pz2;
        */
        {
            // validate assumptions and calculation
    
            ATOOLS::Vec4D P_in    = P_in1 + P_in2;
            ATOOLS::Vec4D P_delta = P_out - P_in;
        
            double mass_in = P_in.Mass();
            double mass_delta = mass_in - mass;

            if ((fabs(P_delta[0]) > 0.001) ||
                (fabs(P_delta[1]) > 0.001) ||
                (fabs(P_delta[2]) > 0.001) ||
                (fabs(P_delta[3]) > 0.001) ||
                (fabs(mass_delta) > 0.001))
            {
                LogMsgWarning( "Massless approximation does not hold in event %i", FMT_I(event.id) );
            }
        }
    }
    
    // define the flavors and momenta
    std::vector<int>     particles;
    ATOOLS::Vec4D_Vector momenta;
    {
        particles.push_back( event.id1 );
        particles.push_back( event.id2 );
        
        momenta.push_back( P_in1 );
        momenta.push_back( P_in2 );

        for (Int_t i = 0; i < event.nparticle; ++i)
        {
            Int_t code = event.kf[i];
            particles.push_back( code );
            
            ATOOLS::Vec4D P_part( event.E[i], event.px[i], event.py[i], event.pz[i] );
            momenta.push_back( P_part );
        }
    }

    // process the event
    m_upSherpaWeight->ProcessEvent( event.id, 2, particles, momenta );    
}
