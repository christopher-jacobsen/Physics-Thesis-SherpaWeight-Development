////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaRootEventFile.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaRootEventFile.h"

// Root includes
#include <TFile.h>
#include <TTree.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaRootEventFile
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaRootEventFile::~SherpaRootEventFile() throw()
{
    Close();    // [noexcept]
}

////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////////////////////////
uint64_t SherpaRootEventFile::Count() const
{
    return static_cast<uint64_t>( std::max(m_nEntries, Long64_t(0)) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool SherpaRootEventFile::ReadEvent( EventFileEvent & event )
{
    if (!m_pTree)
        ThrowError( "ReadEvent() called on closed file." );

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
