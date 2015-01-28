////////////////////////////////////////////////////////////////////////////////////////////////////
//  SMEX.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SMEX_H
#define SMEX_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations
struct SherpaEvent;

namespace SHERPA
{
class Sherpa;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

class SMEX
{
public:
    struct RunParameters
    {
        std::string     inputRootFileName;
        std::string     outputRootFileName;

        std::vector<const char *> argv;
    };

public:
    SMEX();
    ~SMEX() throw();
    
    int ParseCommandLine( int argc, const char * argv[], RunParameters & param );
    
    int Run( const RunParameters & param );
    
private:
    void ProcessEvent( SherpaEvent & event );
    
private:
    std::unique_ptr<SHERPA::Sherpa> m_upSherpa;

private:
    SMEX(const SMEX &) = delete;
    SMEX & operator=(const SMEX &) = delete;
};

#endif // SMEX_H
