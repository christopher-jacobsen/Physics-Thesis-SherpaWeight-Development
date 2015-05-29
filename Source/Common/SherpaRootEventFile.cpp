////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaRootEventFile.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaRootEventFile.h"

// Root includes
#include <TFile.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaEventFileEvent
////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaRootEventFileEvent : public EventFileEvent
{
public:
    SherpaRootEventFileEvent();

    virtual void Clear() override;

    virtual void GetSignalVertex( EventFileVertex & vertex ) const override;

    virtual void SetCoefficients( const DoubleVector & coefs ) override;

private:
    SherpaRootEvent m_event;
    DoubleVector    m_coefs;

    friend SherpaRootEventFile;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaRootEventFile
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaRootEventFile::~SherpaRootEventFile() throw()
{
    Close();    // [noexcept]
}

////////////////////////////////////////////////////////////////////////////////////////////////////
EventFileEvent::UniquePtr SherpaRootEventFile::AllocateEvent() const
{
    return EventFileEvent::UniquePtr( new SherpaRootEventFileEvent );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFile::Open( const std::string & fileName, OpenMode mode )
{
    Close();

    m_fileName = fileName;
    m_mode     = mode;

    try
    {
        Option_t * pOption = (mode == OpenMode::Write) ? "RECREATE" : "";

        m_upFile.reset( new TFile( fileName.c_str(), pOption ) );
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
        // get existing tree

        m_upFile->GetObject( "t3", m_pTree );
        if (!m_pTree)
        {
            LogMsgError( "Failed to load tree (t3) in root file (%hs).", FMT_HS(m_fileName.c_str()) );
            ThrowError( std::invalid_argument( m_fileName ) );
        }

        m_event.SetInputTree( m_pTree );

        m_nEntries = m_pTree->GetEntries();
    }

    if (mode == OpenMode::Write)
    {
        // create output tree

        m_pTree = new TTree( "t3", "SherpaWeight" );   // owned by current directory
        if (m_pTree->IsZombie())
            ThrowError( "Failed to construct output tree." );

        m_pTree->SetDirectory( m_upFile.get() );   // attach to output file, output file now owns tree and will call delete
        m_pTree->SetAutoSave(0);                       // disable autosave

        m_event.SetOutputTree( m_pTree );

        // add entire coefficient vector as a branch
        m_pTree->Branch( "Fij", &m_coefs );
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFile::Close() throw()
{
    m_nEntries  = 0;
    m_iEntry    = 0;
    m_pTree     = nullptr;

    if (m_mode == OpenMode::Write)
    {
        try
        {
            if (m_upFile)
                m_upFile->Write( 0, TFile::kOverwrite );
        }
        catch (...)
        {
            LogMsgError( "Unexpected exception while writing root file (%hs).", FMT_HS(m_fileName.c_str()) );
        }
    }

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

    m_fileName.clear(); // [noexcept]
    m_coefs.clear();    // [noexcept]
}

////////////////////////////////////////////////////////////////////////////////////////////////////
uint64_t SherpaRootEventFile::Count() const
{
    return static_cast<uint64_t>( std::max(m_nEntries, Long64_t(0)) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool SherpaRootEventFile::ReadEvent( EventFileEvent & vEvent )
{
    SherpaRootEventFileEvent & event = static_cast<SherpaRootEventFileEvent &>(vEvent);

    try
    {
        if (!m_pTree)
            ThrowError( "ReadEvent() called on closed file." );

        if (m_iEntry >= m_nEntries)
        {
            event.Clear();  // clear event
            return false;
        }

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
        event.m_event = m_event;

        return true;
    }
    catch (...)
    {
        event.Clear();  // clear event
        throw;          // rethrow
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFile::SetCoefficientNames( const StringVector & coefNames )
{
    if (coefNames.empty())
        ThrowError( "Called SetCoefficientNames() with empty string vector." );

    if (!m_coefs.empty())
        ThrowError( "SetCoefficientNames() must only be called once and before WriteEvent()." );

    m_coefs.resize( coefNames.size() );

    // add each coefficient as an individual branch
    size_t index = 0;
    for (const std::string & name : coefNames)
    {
        size_t pos = name.find( "_", name.find( "_", name.find("_") + 1 ) + 1 );  // find 3rd underscore
        std::string shortName = name.substr(0, pos );
        TBranch * pBranch = m_pTree->Branch( shortName.c_str(), &m_coefs[index] );
        pBranch->SetTitle( (name + "/D").c_str() );
        ++index;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFile::WriteEvent( const EventFileEvent & vEvent )
{
    const SherpaRootEventFileEvent & event = static_cast<const SherpaRootEventFileEvent &>(vEvent);

    if (!m_pTree)
        ThrowError( "WriteEvent() called on closed file." );

    if (m_coefs.empty())
        m_coefs.resize( event.m_coefs.size() );

    if (event.m_coefs.size() != m_coefs.size())
        ThrowError( "Event has " + std::to_string(event.m_coefs.size()) + " coefficients. Expected " + std::to_string(m_coefs.size()) + "." );

    std::copy( event.m_coefs.begin(), event.m_coefs.end(), m_coefs.begin() );  // use copy to ensure buffer not reallocated

    m_event = event.m_event;

    if (m_pTree->Fill() < 0)
        ThrowError( "Fill failed for event id " + std::to_string(m_event.id) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaEventFileEvent
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaRootEventFileEvent::SherpaRootEventFileEvent()
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFileEvent::Clear()
{
    eventId = 0;
    m_event = SherpaRootEvent();
    m_coefs.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFileEvent::GetSignalVertex( EventFileVertex & vertex ) const
{
    vertex = EventFileVertex();  // clear vertex

    size_t nOutput = static_cast<size_t>(m_event.nparticle);
    vertex.output.resize( nOutput );

    double E_out  = 0;
    double pz_out = 0;

    for (size_t i = 0; i < nOutput; ++i)
    {
        EventFileVertex::Particle & particle = vertex.output[i];

        particle.pdg = m_event.kf[i];
        particle.E   = m_event.E [i];
        particle.px  = m_event.px[i];
        particle.py  = m_event.py[i];
        particle.pz  = m_event.pz[i];

        E_out  += particle.E;
        pz_out += particle.pz;
    }

    {
        vertex.input.resize(2);

        EventFileVertex::Particle & in1 = vertex.input[0];
        EventFileVertex::Particle & in2 = vertex.input[1];

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
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaRootEventFileEvent::SetCoefficients( const DoubleVector & coefs )
{
    m_coefs = coefs;
}

