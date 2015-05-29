////////////////////////////////////////////////////////////////////////////////////////////////////
//  HepMCEventFile.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "HepMCEventFile.h"

#include "common.h"

// HepMC includes
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class HepMCEventFileEvent
////////////////////////////////////////////////////////////////////////////////////////////////////

class HepMCEventFileEvent : public EventFileEvent
{
public:
    HepMCEventFileEvent();

    virtual void Clear() override;

    virtual void GetSignalVertex( EventFileVertex & vertex ) const override;

    virtual void SetCoefficients( const DoubleVector & coefs ) override;

private:
    static EventFileVertex::Particle ConvertParticle( const HepMC::GenParticle & part );

private:
    std::unique_ptr<HepMC::GenEvent>    m_upGenEvent;

    friend HepMCEventFile;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// class HepMCEventFile
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
HepMCEventFile::HepMCEventFile()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
HepMCEventFile::~HepMCEventFile() throw()
{
    Close();    // [noexcept]
}

////////////////////////////////////////////////////////////////////////////////////////////////////
EventFileEvent::UniquePtr HepMCEventFile::AllocateEvent() const
{
    return EventFileEvent::UniquePtr( new HepMCEventFileEvent );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFile::Open( const std::string & fileName, OpenMode mode )
{
    Close();

    m_fileName = fileName;

    try
    {
        std::ios::openmode ioMode = (mode == OpenMode::Write) ? std::ios::out : std::ios::in;
        m_upIO.reset( new HepMC::IO_GenEvent( fileName, ioMode ) );
    }
    catch (...)
    {
        LogMsgError( "Failed to construct HepMC IO object for file (%hs).", FMT_HS(m_fileName.c_str()) );
        throw;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFile::Close() throw()
{
    try
    {
        m_upIO.reset();
    }
    catch (...)
    {
        LogMsgError( "Unexpected exception while destructing HepMC IO object." );
    }

    m_fileName.clear();     // [noexcept]
}

////////////////////////////////////////////////////////////////////////////////////////////////////
uint64_t HepMCEventFile::Count() const
{
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool HepMCEventFile::ReadEvent( EventFileEvent & vEvent )
{
    HepMCEventFileEvent & event = static_cast<HepMCEventFileEvent &>(vEvent);

    event.Clear();  // clear event

    if (!m_upIO)
        ThrowError( "ReadEvent() called on closed file." );

    HepMC::GenEvent & genEvent = *event.m_upGenEvent;

    if (!m_upIO->fill_next_event( &genEvent ))
        return false;  // no more events

    event.eventId = genEvent.event_number();

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFile::SetCoefficientNames( const StringVector & coefNames )
{
    if (coefNames.empty())
        ThrowError( "Called SetCoefficientNames() with empty string vector." );

//    if (!m_coefs.empty())
//        ThrowError( "SetCoefficientNames() must only be called once and before WriteEvent()." );

//    m_coefs.resize( coefNames.size() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFile::WriteEvent( const EventFileEvent & vEvent )
{
    const HepMCEventFileEvent & event = static_cast<const HepMCEventFileEvent &>(vEvent);

    if (!event.m_upGenEvent)
        ThrowError( "WriteEvent() called on uninitialized event." );

    if (!m_upIO)
        ThrowError( "WriteEvent() called on closed file." );

    m_upIO->write_event( event.m_upGenEvent.get() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class HepMCEventFileEvent
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
HepMCEventFileEvent::HepMCEventFileEvent()
{
    Clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFileEvent::Clear()
{
    eventId = 0;

    if (!m_upGenEvent)
        m_upGenEvent.reset( new HepMC::GenEvent );

    m_upGenEvent->clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFileEvent::GetSignalVertex( EventFileVertex & vertex ) const
{
    vertex = EventFileVertex();  // clear vertex

    if (!m_upGenEvent)
        ThrowError( "GetSignalVertex() called on uninitialized event." );

    // get signal process vertex
    HepMC::GenVertex * pSignal = m_upGenEvent->signal_process_vertex();
    if (!pSignal)
        ThrowError( "Missing signal vertex for event." );

    // input events
    {
        vertex.input.reserve( (size_t) std::max(pSignal->particles_in_size(), 0) );

        auto itrGenPart = pSignal->particles_in_const_begin();
        auto endGenPart = pSignal->particles_in_const_end();
        for ( ; itrGenPart != endGenPart; ++itrGenPart)
        {
            const HepMC::GenParticle * pGenPart = *itrGenPart;
            vertex.input.push_back( ConvertParticle( *pGenPart ) );
        }
    }

    // output events
    {
        vertex.output.reserve( (size_t) std::max(pSignal->particles_out_size(), 0) );

        auto itrGenPart = pSignal->particles_out_const_begin();
        auto endGenPart = pSignal->particles_out_const_end();
        for ( ; itrGenPart != endGenPart; ++itrGenPart)
        {
            const HepMC::GenParticle * pGenPart = *itrGenPart;
            vertex.output.push_back( ConvertParticle( *pGenPart ) );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void HepMCEventFileEvent::SetCoefficients( const DoubleVector & coefs )
{
    HepMC::WeightContainer & weights = m_upGenEvent->weights();

    for (double c : coefs)
        weights.push_back(c);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
EventFileVertex::Particle HepMCEventFileEvent::ConvertParticle( const HepMC::GenParticle & part )
{
    EventFileVertex::Particle result;

	HepMC::FourVector mom = part.momentum();

    result.pdg = part.pdg_id();
    result.E   = mom.e();
    result.px  = mom.px();
    result.py  = mom.py();
    result.pz  = mom.pz();

    return result;
}
