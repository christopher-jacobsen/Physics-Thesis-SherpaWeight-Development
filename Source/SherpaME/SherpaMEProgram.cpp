////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaMEProgram.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaMEProgram.h"
#include "SherpaRootEvent.h"
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

struct EventFileEvent
{
    struct Particle
    {
        int32_t pdg;

        double  E;
        double  px;
        double  py;
        double  pz;
    };

    int32_t                 eventId;
    std::vector<Particle>   input;
    std::vector<Particle>   output;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

struct EventFileInterface
{
    enum class OpenMode
    {
        Read,
        Write
    };

    virtual ~EventFileInterface() throw()   = default;

    virtual void Open( const std::string & fileName, OpenMode mode )   = 0;
    virtual void Close() throw()        = 0;

    virtual uint64_t Count() const      = 0;

    virtual bool ReadEvent( EventFileEvent & event ) = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaRootEventFile : public EventFileInterface
{
public:
    static bool IsSupported( const std::string & fileName ) throw();

    SherpaRootEventFile() = default;
    virtual ~SherpaRootEventFile() throw() override;

    virtual void Open( const std::string & fileName, OpenMode mode ) override;
    virtual void Close() throw() override;

    virtual uint64_t Count() const override;

    virtual bool ReadEvent( EventFileEvent & event ) override;

private:
    std::string             m_fileName;

    std::unique_ptr<TFile>  m_upFile;
    TTree *                 m_pTree     = nullptr;

    SherpaRootEvent         m_event;
    Long64_t                m_nEntries  = 0;
    Long64_t                m_iEntry    = 0;
};

SherpaRootEventFile::~SherpaRootEventFile() throw()
{
    Close();    // [noexcept]
}

void SherpaRootEventFile::Open( const std::string & fileName, OpenMode mode )
{
    Close();

    m_fileName = fileName;

    try
    {
        m_upFile.reset( new TFile( fileName.c_str() ) );
    }
    catch (...)
    {
        LogMsgError( "Failed to construct root object for file (%hs).", FMT_HS(m_fileName.c_str()) );
        throw;
    }

    if (m_upFile->IsZombie() || !m_upFile->IsOpen())      // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to open root file (%hs).", FMT_HS(m_fileName.c_str()) );
        ThrowError( std::invalid_argument( m_fileName ) );
    }

    if (mode == OpenMode::Read)
    {
        m_upFile->GetObject( "t3", m_pTree );
        if (!m_pTree)
        {
            LogMsgError( "Failed to load tree (t3) in root file (%hs).", FMT_HS(m_fileName.c_str()) );
            ThrowError( std::invalid_argument( m_fileName ) );
        }

        m_event.SetInputTree( m_pTree );

        m_nEntries = m_pTree->GetEntries();
    }
}

void SherpaRootEventFile::Close() throw()
{
    m_nEntries  = 0;
    m_iEntry    = 0;
    m_pTree     = nullptr;

    try
    {
        if (m_upFile)
            m_upFile->Close();
    }
    catch (...)
    {
        LogMsgError( "Unexpected exception while closing root file (%hs).", FMT_HS(m_fileName.c_str()) );
    }

    try
    {
        m_upFile.reset();
    }
    catch (...)
    {
        LogMsgError( "Unexpected exception while destructing root file object." );
    }

    m_fileName.clear();     // [noexcept]
}

uint64_t SherpaRootEventFile::Count() const
{
    return static_cast<uint64_t>( std::max(m_nEntries, Long64_t(0)) );
}

bool SherpaRootEventFile::ReadEvent( EventFileEvent & event )
{
    if (!m_pTree)
        ThrowError( "Next() called on closed root file." );

    if (m_iEntry >= m_nEntries)
        return false;

    if (m_pTree->LoadTree(m_iEntry) < 0)
        ThrowError( "LoadTree failed on entry " + std::to_string(m_iEntry+1) );

    if (m_pTree->GetEntry(m_iEntry) < 0)
        ThrowError( "GetEntry failed on entry " + std::to_string(m_iEntry+1) );

    ++m_iEntry;

    // fill in event

    if (m_event.nparticle <= 0)
        ThrowError( "No outgoing particles in event id " + std::to_string(m_event.id) );

    if ((size_t)m_event.nparticle > SherpaRootEvent::max_nparticle )
    {
        LogMsgError( "Number of outgoing particles in event id %i exceeds maximum. nparticles=%i (max %u).",
                     FMT_I(m_event.id), FMT_I(m_event.nparticle), FMT_U(SherpaRootEvent::max_nparticle) );
        ThrowError( "Number of outgoing particles in event exceeds maximum." );
    }

    event.eventId = m_event.id;

    size_t nOutput = static_cast<size_t>(m_event.nparticle);
    event.output.resize( nOutput );

    double E_out  = 0;
    double pz_out = 0;

    for (size_t i = 0; i < nOutput; ++i)
    {
        EventFileEvent::Particle & particle = event.output[i];

        particle.pdg = m_event.kf[i];
        particle.E   = m_event.E [i];
        particle.px  = m_event.px[i];
        particle.py  = m_event.py[i];
        particle.pz  = m_event.pz[i];

        E_out  += particle.E;
        pz_out += particle.pz;
    }

    {
        event.input.resize(2);

        EventFileEvent::Particle & in1 = event.input[0];
        EventFileEvent::Particle & in2 = event.input[1];

        in1.pdg = m_event.id1;
        in2.pdg = m_event.id2;

        // assume massless incoming partons
        double E1 = (E_out + pz_out) / 2;
        double E2 = (E_out - pz_out) / 2;

        in1.E  =  E1;
        in2.E  =  E2;

        in1.pz =  E1;
        in2.pz = -E2;

        in1.px = 0;
        in2.px = 0;

        in1.py = 0;
        in2.py = 0;

        /* TODO
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
        */
    }

    return true;
}

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

            runArgv.push_back( "INIT_ONLY=2" ); // prevent Sherpa from starting the cross section integration

            if (!m_upSherpa->InitializeTheRun( static_cast<int>(runArgv.size()), const_cast<char **>(runArgv.data()) ))
                ThrowError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        }

        // open input file

        LogMsgInfo( "Input file : %hs", FMT_HS(param.inputRootFileName.c_str()) );
        SherpaRootEventFile inputFile;
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

        const uint64_t nEntries     = inputFile.Count();
        const uint64_t logFrequency = std::max( nEntries / 10, uint64_t(1) );

        LogMsgInfo( "\nGetting matrix elements for %llu events ...", FMT_LLU(nEntries) );

        time_t timeStartProcess = time(nullptr);
        
        for (uint64_t iEntry = 0; inputFile.ReadEvent( inputEvent ); ++iEntry)
        {
            if (!ProcessEvent( inputEvent, outputEvent ))
                continue;

            if (iEntry % logFrequency == 0)
                LogMsgInfo( "Event %llu (id %i): ME = %E", FMT_LLU(iEntry+1), FMT_I(inputEvent.eventId), FMT_F(outputEvent.me) );

            if (pOutputTree->Fill() < 0)
                ThrowError( "Fill failed on entry " + std::to_string(iEntry+1) );
        }

        // write and close the output file (not really necessary as would be done in destructor)
        
        upOutputFile->Write( 0, TFile::kOverwrite );
        upOutputFile->Close();

        time_t timeStopProcess = time(nullptr);

        LogMsgInfo( "%llu events completed. (%u seconds)", FMT_LLU(nEntries), FMT_U(timeStopProcess - timeStartProcess) );
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
