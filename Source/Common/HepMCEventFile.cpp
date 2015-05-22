////////////////////////////////////////////////////////////////////////////////////////////////////
//  HepMCEventFile.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "HepMCEventFile.h"

#include "common.h"

// HepMC includes
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>

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
    return 0; // TODO
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool HepMCEventFile::ReadEvent( EventFileEvent & event )
{
    event = EventFileEvent();  // clear event

    if (!m_upIO)
        ThrowError( "ReadEvent() called on closed file." );

    HepMC::GenEvent genEvent;

    if (!m_upIO->fill_next_event( &genEvent ))
        return false;  // no more events

    event.eventId = genEvent.event_number();

    // get signal process vertex
    HepMC::GenVertex * pSignal = genEvent.signal_process_vertex();
    if (!pSignal)
        ThrowError( "Missing signal vertex for event" );

    // input events
    {
        event.input.reserve( (size_t) std::max(pSignal->particles_in_size(), 0) );

        auto itrGenPart = pSignal->particles_in_const_begin();
        auto endGenPart = pSignal->particles_in_const_end();
        for ( ; itrGenPart != endGenPart; ++itrGenPart)
        {
            const HepMC::GenParticle * pGenPart = *itrGenPart;
            event.input.push_back( ConvertParticle( *pGenPart ) );
        }
    }

    // output events
    {
        event.output.reserve( (size_t) std::max(pSignal->particles_out_size(), 0) );

        auto itrGenPart = pSignal->particles_out_const_begin();
        auto endGenPart = pSignal->particles_out_const_end();
        for ( ; itrGenPart != endGenPart; ++itrGenPart)
        {
            const HepMC::GenParticle * pGenPart = *itrGenPart;
            event.output.push_back( ConvertParticle( *pGenPart ) );
        }
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
EventFileEvent::Particle HepMCEventFile::ConvertParticle( const HepMC::GenParticle & part )
{
    EventFileEvent::Particle result;

	HepMC::FourVector mom = part.momentum();

    result.pdg = part.pdg_id();
    result.E   = mom.e();
    result.px  = mom.px();
    result.py  = mom.py();
    result.pz  = mom.pz();

    return result;
}

