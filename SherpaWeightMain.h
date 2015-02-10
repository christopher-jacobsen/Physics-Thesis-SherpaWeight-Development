////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightMain.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAWEIGHTMAIN_H
#define SHERPAWEIGHTMAIN_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

struct SherpaRootEvent;
class  SherpaWeight;

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
    std::unique_ptr<SherpaWeight> m_upSherpaWeight;

private:
    SherpaWeightMain(const SherpaWeightMain &) = delete;
    SherpaWeightMain & operator=(const SherpaWeightMain &) = delete;
};

#endif // SHERPAWEIGHTMAIN_H
