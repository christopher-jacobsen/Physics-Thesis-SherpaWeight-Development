////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaRootEventFile.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPA_ROOT_EVENT_FILE_H
#define SHERPA_ROOT_EVENT_FILE_H

#include "EventFile.h"

#include "common.h"
#include "SherpaRootEvent.h"

#include <Rtypes.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

// Root classes
class TFile;
class TTree;

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

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPA_ROOT_EVENT_FILE_H
