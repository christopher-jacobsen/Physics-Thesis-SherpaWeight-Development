////////////////////////////////////////////////////////////////////////////////////////////////////
//  EventFile.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_FILE_H
#define EVENT_FILE_H

#include "common.h"

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

    virtual void Open( const std::string & fileName, OpenMode mode )    = 0;
    virtual void Close() throw()                                        = 0;

    virtual uint64_t Count() const                                      = 0;

    virtual bool ReadEvent( EventFileEvent & event )                    = 0;
};

#endif // EVENT_FILE_H
