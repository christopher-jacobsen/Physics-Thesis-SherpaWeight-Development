////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaMEProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaMEProgram.h"

#include "SherpaRootEventFile.h"
#include "HepMCEventFile.h"
#include "MERootEvent.h"

#include "SherpaMECalculator.h"

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Sherpa and Root include files

// Sherpa includes
#include <SHERPA/Main/Sherpa.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Math/Vector.H>

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
        time_t timeStartRun = time(nullptr);

        // initialize sherpa
        {
            std::vector<const char *> runArgv(param.argv);

            runArgv.push_back( "INIT_ONLY=2"   );   // prevent Sherpa from starting the cross section integration
            runArgv.push_back( "EVENT_OUTPUT=" );   // prevent Sherpa from overwriting any event file specified in the dat file

            if (!m_upSherpa->InitializeTheRun( static_cast<int>(runArgv.size()), const_cast<char **>(runArgv.data()) ))
                ThrowError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        }

        // open input file

        LogMsgInfo( "Input file : %hs", FMT_HS(param.inputRootFileName.c_str()) );
        //SherpaRootEventFile inputFile;
        HepMCEventFile inputFile;
        inputFile.Open( param.inputRootFileName, EventFileInterface::OpenMode::Read );

        // create output file

        LogMsgInfo( "Output file: %hs", FMT_HS(param.outputRootFileName.c_str()) );
        std::unique_ptr<TFile> upOutputFile( new TFile( param.outputRootFileName.c_str(), "RECREATE" ) );
        if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
        {
            LogMsgError( "Failed to create output file (%hs).", FMT_HS(param.outputRootFileName.c_str()) );
            ThrowError( std::invalid_argument( param.outputRootFileName ) );
        }

        // create output tree

        TTree * pOutputTree( new TTree( "SherpaME", "SherpaME" ) );   // owned by current directory
        if (pOutputTree->IsZombie())
            ThrowError("Failed to construct output tree.");

        pOutputTree->SetDirectory( upOutputFile.get() );   // attach to output file, output file now owns tree and will call delete
        pOutputTree->SetAutoSave(0);                       // disable autosave

        // create event containers and connect/declare variables

        EventFileEvent  inputEvent;
        MERootEvent     outputEvent;

        outputEvent.SetOutputTree( pOutputTree );
        
        // loop through and process each input event

        uint64_t    iEvent          = 0;
        uint64_t    nEvents         = inputFile.Count();
        uint64_t    logFrequency    = 1;
        uint32_t    logCount        = 0;

        if (nEvents)
            LogMsgInfo( "\nGetting matrix elements for %llu events ...", FMT_LLU(nEvents) );
        else
            LogMsgInfo( "\nGetting matrix elements for events ..." );

        time_t timeStartProcess = time(nullptr);
        
        for ( ; inputFile.ReadEvent( inputEvent ); ++iEvent)
        {
            if (!ProcessEvent( inputEvent, outputEvent ))
                continue;

            if ((iEvent + 1) % logFrequency == 0)
            {
                LogMsgInfo( "Event %llu (id %i): ME = %E", FMT_LLU(iEvent+1), FMT_I(inputEvent.eventId), FMT_F(outputEvent.me) );
                if (++logCount == 10)
                {
                    logFrequency *= 10;
                    logCount      = 1;
                }
            }

            if (pOutputTree->Fill() < 0)
                ThrowError( "Fill failed on event " + std::to_string(iEvent+1) );
        }

        // write and close the output file (not really necessary as would be done in destructor)
        
        upOutputFile->Write( 0, TFile::kOverwrite );
        upOutputFile->Close();

        time_t timeStopProcess = time(nullptr);

        LogMsgInfo( "%llu events completed. (%u seconds)", FMT_LLU(iEvent), FMT_U(timeStopProcess - timeStartProcess) );
        LogMsgInfo( "Done. (%u seconds)", FMT_U(timeStopProcess - timeStartRun) );

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
bool SherpaMEProgram::ProcessEvent( const EventFileEvent & inputEvent, MERootEvent & outputEvent )
{
    outputEvent = MERootEvent();  // clear values

    // define the flavors and momenta
    std::vector<int>     particles;
    ATOOLS::Vec4D_Vector momenta;
    {
        for (const EventFileEvent::Particle & part : inputEvent.input)
        {
            ATOOLS::Vec4D P_part( part.E, part.px, part.py, part.pz );

            particles.push_back( part.pdg );
            momenta  .push_back( P_part   );
        }

        for (const EventFileEvent::Particle & part : inputEvent.output)
        {
            ATOOLS::Vec4D P_part( part.E, part.px, part.py, part.pz );

            particles.push_back( part.pdg );
            momenta  .push_back( P_part   );
        }
    }

    // get the matrix element for the event

    double me = GetEventME( inputEvent.input.size(), particles, momenta );

    outputEvent.id = inputEvent.eventId;
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
