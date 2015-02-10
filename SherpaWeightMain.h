////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightMain.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAWEIGHTMAIN_H
#define SHERPAWEIGHTMAIN_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations
struct SherpaRootEvent;

namespace SHERPA
{
class Sherpa;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeightMain
{
public:
    struct RunParameters
    {
        std::string     inputRootFileName;
        std::string     outputRootFileName;

        std::vector<const char *> argv;
    };

public:
    SherpaWeightMain();
    ~SherpaWeightMain() throw();
    
    int ParseCommandLine( int argc, const char * argv[], RunParameters & param );
    
    int Run( const RunParameters & param );
    
private:
    void ProcessEvent( SherpaRootEvent & event );
    
private:
    std::unique_ptr<SHERPA::Sherpa> m_upSherpa;

private:
    SherpaWeightMain(const SherpaWeightMain &) = delete;
    SherpaWeightMain & operator=(const SherpaWeightMain &) = delete;
};

#endif // SHERPAWEIGHTMAIN_H
