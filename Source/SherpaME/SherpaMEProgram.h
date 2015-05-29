////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaMEProgram.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPA_ME_PROGRAM_H
#define SHERPA_ME_PROGRAM_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

namespace SHERPA
{
class Sherpa;
}

namespace ATOOLS
{
template<typename>
class Vec4;

typedef Vec4<double>            Vec4D;
typedef std::vector<Vec4D>      Vec4D_Vector;
}

struct EventFileVertex;
struct MERootEvent;

////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaMEProgram
{
public:
    struct RunParameters
    {
        std::string     inputRootFileName;
        std::string     outputRootFileName;

        std::vector<const char *> argv;
    };

public:
    SherpaMEProgram();
    ~SherpaMEProgram() throw();
    
    int ParseCommandLine( int argc, const char * argv[], RunParameters & param );
    
    int Run( const RunParameters & param );
    
private:
    bool ProcessEvent( int32_t eventId, const EventFileVertex & inputEvent, MERootEvent & outputEvent );     // returns true if valid outputEvent generated

    double GetEventME( size_t nInParticles, const std::vector<int> & particleCodes, const ATOOLS::Vec4D_Vector & particleMomenta );
    
private:
    std::unique_ptr<SHERPA::Sherpa> m_upSherpa;

private:
    SherpaMEProgram(const SherpaMEProgram &) = delete;
    SherpaMEProgram & operator=(const SherpaMEProgram &) = delete;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPA_ME_PROGRAM_H
