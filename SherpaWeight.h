////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeight.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAWEIGHT_H
#define SHERPAWEIGHT_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations
struct SherpaEvent;

namespace SHERPA
{
class Sherpa;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeight
{
public:
    struct RunParameters
    {
        std::string     inputRootFileName;
        std::string     outputRootFileName;

        std::vector<const char *> argv;
    };

public:
    SherpaWeight();
    ~SherpaWeight() throw();
    
    int ParseCommandLine( int argc, const char * argv[], RunParameters & param );
    
    int Run( const RunParameters & param );
    
private:
    void ProcessEvent( SherpaEvent & event );
    
private:
    std::unique_ptr<SHERPA::Sherpa> m_upSherpa;

private:
    SherpaWeight(const SherpaWeight &) = delete;
    SherpaWeight & operator=(const SherpaWeight &) = delete;
};

#endif // SHERPAWEIGHT_H
