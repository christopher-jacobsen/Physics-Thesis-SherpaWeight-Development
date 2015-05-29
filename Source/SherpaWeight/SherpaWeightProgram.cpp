////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeightProgram.h"
#include "SherpaWeight.h"

#include "SherpaRootEventFile.h"
#include "HepMCEventFile.h"

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// Sherpa and Root include files

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
    //SherpaRootEventFile inputFile;
    HepMCEventFile inputFile;
    inputFile.Open( param.inputRootFileName, EventFileInterface::OpenMode::Read );

    // open output file

    LogMsgInfo( "Output file: %hs", FMT_HS(param.outputRootFileName.c_str()) );
    //SherpaRootEventFile outputFile;
    HepMCEventFile outputFile;
    outputFile.Open( param.outputRootFileName, EventFileInterface::OpenMode::Write );

    // add coefficient output variables
    
    const SherpaWeight::StringVector &  coefNames = m_upSherpaWeight->CoefficientNames();
    const size_t                        nCoefs    = m_upSherpaWeight->NCoefficients();

    outputFile.SetCoefficientNames( coefNames );

    // create event containers

    EventFileEvent::UniquePtr   upCurrentEvent  = inputFile.AllocateEvent();
    EventFileEvent &            currentEvent    = *upCurrentEvent;

    // loop through and process each input event

    uint64_t    iEvent          = 1;
    uint64_t    nEvents         = inputFile.Count();
    uint64_t    logFrequency    = 1;
    uint32_t    logCount        = 0;

    if (nEvents)
        LogMsgInfo( "\nCalculating coefficients for %llu events ...", FMT_LLU(nEvents) );
    else
        LogMsgInfo( "\nCalculating coefficients for events ..." );

    time_t timeStartProcess = time(nullptr);

    for ( ; inputFile.ReadEvent( currentEvent ); ++iEvent)
    {
        {
            SherpaWeight::DoubleVector coefs = m_upSherpaWeight->CoefficientValues( currentEvent.eventId );
            if (coefs.empty())
            {
                LogMsgWarning( "No coefficients for event %llu (id %i).", FMT_LLU(iEvent), FMT_I(currentEvent.eventId) );
            }

            if (coefs.size() != nCoefs)
            {
                LogMsgWarning( "Missing coefficients for event %llu (id %i). Setting all to zero.", FMT_LLU(iEvent), FMT_I(currentEvent.eventId) );
                coefs.resize(nCoefs);
                std::fill( coefs.begin(), coefs.end(), 0.0 );
            }

            if (iEvent % logFrequency == 0)
            {
                if (++logCount == 10)
                {
                    logFrequency *= 10;
                    logCount      = 1;
                }

                LogMsgInfo( "Event %llu (id %i):", FMT_LLU(iEvent), FMT_I(currentEvent.eventId) );

                size_t index = 0;
                for (double value : coefs)
                {
                    const char * pName = (index < coefNames.size()) ? coefNames[index].c_str() : "Unknown";
                    LogMsgInfo( "%-20hs:\t%.15E", FMT_HS(pName), FMT_F(value) );
                    ++index;
                }
                
                LogMsgInfo( "" );
            }

            currentEvent.SetCoefficients( coefs );
        }

        outputFile.WriteEvent( currentEvent );
    }

    outputFile.Close(); // Close flushes events to disk

    time_t timeStopProcess = time(nullptr);

    LogMsgInfo( "%llu events completed. (%u seconds)", FMT_LLU(iEvent-1), FMT_U(timeStopProcess - timeStartProcess) );
}
