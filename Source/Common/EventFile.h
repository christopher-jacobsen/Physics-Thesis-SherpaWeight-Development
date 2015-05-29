////////////////////////////////////////////////////////////////////////////////////////////////////
//  EventFile.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_FILE_H
#define EVENT_FILE_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

struct EventFileVertex
{
    struct Particle
    {
        int32_t pdg     = 0;

        double  E       = 0;
        double  px      = 0;
        double  py      = 0;
        double  pz      = 0;
    };

    std::vector<Particle>   input;
    std::vector<Particle>   output;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

struct EventFileEvent
{
    int32_t eventId = 0;

public:
    typedef std::vector<double>                 DoubleVector;
    typedef std::unique_ptr<EventFileEvent>     UniquePtr;

    virtual ~EventFileEvent() throw()                                   = default;

    virtual void Clear()                                                = 0;

    virtual void GetSignalVertex( EventFileVertex & vertex ) const      = 0;

    virtual void SetCoefficients( const DoubleVector & coefs )          = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

struct EventFileInterface
{
    enum class OpenMode
    {
        Read,
        Write
    };

    typedef std::vector<std::string> StringVector;

public:
    virtual ~EventFileInterface() throw()                               = default;

    virtual EventFileEvent::UniquePtr AllocateEvent() const             = 0;

    virtual void Open( const std::string & fileName, OpenMode mode )    = 0;
    virtual void Close() throw()                                        = 0;

    // reading
    virtual uint64_t Count() const                                      = 0;

    virtual bool ReadEvent( EventFileEvent & event )                    = 0;  // returns false if no more events

    // writing
    virtual void SetCoefficientNames( const StringVector & coefNames )  = 0;

    virtual void WriteEvent( const EventFileEvent & event )             = 0;
};

#endif // EVENT_FILE_H
