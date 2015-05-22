////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeightProgram.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPA_WEIGHT_PROGRAM_H
#define SHERPA_WEIGHT_PROGRAM_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

class SherpaWeight;

////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeightProgram
{
public:
    struct RunParameters
    {
        std::string     inputRootFileName;
        std::string     outputRootFileName;

        std::vector<const char *> argv;
    };

public:
    SherpaWeightProgram();
    ~SherpaWeightProgram() throw();
    
    int ParseCommandLine( int argc, const char * argv[], RunParameters & param );
    
    int Run( const RunParameters & param );

private:
    void SaveCoefficients( const RunParameters & param );
    
private:
    std::unique_ptr<SherpaWeight> m_upSherpaWeight;

private:
    SherpaWeightProgram(const SherpaWeightProgram &) = delete;
    SherpaWeightProgram & operator=(const SherpaWeightProgram &) = delete;
};

#endif // SHERPA_WEIGHT_PROGRAM_H
