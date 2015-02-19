////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaMEProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaMEProgram.h"
#include "SherpaRootEvent.h"
#include "MERootEvent.h"
#include "SherpaMECalculator.h"

#include "common.h"

// Sherpa includes
#include <SHERPA/Main/Sherpa.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Math/Vector.H>

// TODO: pragma warnings away?
// Root includes
#include <TFile.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaMEProgram
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaMEProgram::SherpaMEProgram()
    : m_upSherpa( new SHERPA::Sherpa )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaMEProgram::~SherpaMEProgram() throw()
{
    try
    {
        m_upSherpa.reset();
    }
    catch (...)
    {
        LogMsgError( "Unexpected exception while terminating Sherpa framework." );
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SherpaMEProgram::ParseCommandLine( int argc, const char * argv[], RunParameters & param )
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
    LogMsgInfo("Usage: SherpaME input_root_file output_root_file <sherpa_arguments ...>");
    return -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int SherpaMEProgram::Run( const RunParameters & param )
{
    try
    {
        // initialize sherpa
        
        if (!m_upSherpa->InitializeTheRun( static_cast<int>(param.argv.size()), const_cast<char **>(param.argv.data()) ))
            ThrowError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        
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

        TTree * pOutputTree( new TTree( "SherpaME", "SherpaME" ) );   // owned by current directory
        if (pOutputTree->IsZombie())
            ThrowError("Failed to construct output tree.");

        pOutputTree->SetDirectory( upOutputFile.get() );   // attach to output file, output file now owns tree and will call delete
        pOutputTree->SetAutoSave(0);                       // disable autosave

        // create event containers and connect/declare variables

        SherpaRootEvent inputEvent;
        MERootEvent     outputEvent;
        
        inputEvent.SetInputTree(   pInputTree  );
        outputEvent.SetOutputTree( pOutputTree );
        
        // loop through and process each input event
        
        const Long64_t nEntries = pInputTree->GetEntries();
        
        for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)
        {
            if (pInputTree->LoadTree(iEntry) < 0)
                ThrowError( "LoadTree failed on entry " + std::to_string(iEntry) );
            
            if (pInputTree->GetEntry(iEntry) < 0)
                ThrowError( "GetEntry failed on entry " + std::to_string(iEntry) );
            
            if (!ProcessEvent( inputEvent, outputEvent ))
                continue;
            
            if (pOutputTree->Fill() < 0)
                ThrowError( "Fill failed on entry " + std::to_string(iEntry) );
        }
        
        // write and close the output file (not really necessary as would be done in destructor)
        
        upOutputFile->Write( 0, TFile::kOverwrite );
        upOutputFile->Close();

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
bool SherpaMEProgram::ProcessEvent( const SherpaRootEvent & inputEvent, MERootEvent & outputEvent )
{
    outputEvent = MERootEvent();  // clear values
    
    // validate event
    
    if (inputEvent.nparticle > SherpaRootEvent::max_nparticle )
    {
        LogMsgError( "Number of outgoing particles in event exceeds maximum. nparticles=%i (max %u).",
                     FMT_I(inputEvent.nparticle), FMT_U(SherpaRootEvent::max_nparticle) );
        return false;
    }

    if (inputEvent.nparticle <= 0)
    {
        LogMsgWarning( "Skipping event %i. No outgoing particles.", FMT_I(inputEvent.id) );
        return false;
    }

    // calculate incoming 4-momentum
    ATOOLS::Vec4D P_in1, P_in2;
    {
        ATOOLS::Vec4D P_out;
        for (Int_t i = 0; i < inputEvent.nparticle; ++i)
        {
            ATOOLS::Vec4D P_part( inputEvent.E[i], inputEvent.px[i], inputEvent.py[i], inputEvent.pz[i] );
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
                LogMsgWarning( "Massless approximation does not hold in event %i", FMT_I(inputEvent.id) );
            }
        }
    }
    
    // define the flavors and momenta
    std::vector<int>     particles;
    ATOOLS::Vec4D_Vector momenta;
    {
        particles.push_back( inputEvent.id1 );
        particles.push_back( inputEvent.id2 );
        
        momenta.push_back( P_in1 );
        momenta.push_back( P_in2 );

        for (Int_t i = 0; i < inputEvent.nparticle; ++i)
        {
            Int_t code = inputEvent.kf[i];
            particles.push_back( code );
            
            ATOOLS::Vec4D P_part( inputEvent.E[i], inputEvent.px[i], inputEvent.py[i], inputEvent.pz[i] );
            momenta.push_back( P_part );
        }
    }

    // get the matrix element for the event

    double me = GetEventME( 2, particles, momenta );

    LogMsgInfo( "Event %i: ME=%E", FMT_I(inputEvent.id), FMT_F(me) );
    
    outputEvent.id = inputEvent.id;
    outputEvent.me = me;
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
double SherpaMEProgram::GetEventME( size_t nInParticles, const std::vector<int> & particleCodes, const ATOOLS::Vec4D_Vector & particleMomenta )
{
    if (!nInParticles || (nInParticles > particleCodes.size()))
        ThrowError( std::invalid_argument( "GetEventME: invalid number of input particles " + std::to_string(nInParticles) ) );
    
    if (particleCodes.size() != particleMomenta.size())
        ThrowError( std::invalid_argument( "GetEventME: mismatch in number of particle codes and momenta" ) );
    
    SherpaMECalculator meCalc( m_upSherpa.get() );

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
        // TODO: remove this catch, refactor Initialize to return a bool
    
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
    
    return meCalc.MatrixElement();
}
