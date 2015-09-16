////////////////////////////////////////////////////////////////////////////////////////////////////
//  HepMCEventFile.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HEPMC_EVENT_FILE_H
#define HEPMC_EVENT_FILE_H

#include "EventFile.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

namespace HepMC
{
class IO_GenEvent;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class HepMCEventFile
////////////////////////////////////////////////////////////////////////////////////////////////////

class HepMCEventFile : public EventFileInterface
{
public:
    static bool IsSupported( const std::string & fileName ) throw();

    HepMCEventFile();
    virtual ~HepMCEventFile() throw() override;

    virtual EventFileEvent::UniquePtr AllocateEvent() const override;

    virtual void Open( const std::string & fileName, OpenMode mode ) override;
    virtual void Close() throw() override;

    virtual uint64_t Count() const override;

    virtual bool ReadEvent( EventFileEvent & event ) override;

    virtual void SetCoefficientNames( const StringVector & coefNames )  override;

    virtual void WriteEvent( const EventFileEvent & event ) override;

private:
    std::string                             m_fileName;
    std::unique_ptr<std::istream>           m_upIStream;
    std::unique_ptr<std::ostream>           m_upOStream;
    std::unique_ptr<HepMC::IO_GenEvent>     m_upIO;
    StringVector                            m_coefNames;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // HEPMC_EVENT_FILE_H
