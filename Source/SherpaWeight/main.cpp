////////////////////////////////////////////////////////////////////////////////////////////////////
//  main.cpp
//  SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeightProgram.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    try
    {
        int result = 0;

        SherpaWeightProgram SherpaWeight;
        
        SherpaWeightProgram::RunParameters param;
        result = SherpaWeight.ParseCommandLine( argc, argv, param );
        if (result != 0)
            return result;
        
        return SherpaWeight.Run( param );
    }
    catch (const std::exception & error)
    {
        LogMsgError( "Exception: %hs", FMT_HS(error.what()) );
    }
    catch (...)
    {
        LogMsgError( "Unknown Exception!" );
    }

    return EXIT_FAILURE;
}
