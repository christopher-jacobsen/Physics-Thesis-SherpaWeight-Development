//
//  main.cpp
//  SherpaWeight

#include "common.h"
#include "SherpaWeightMain.h"

int main(int argc, const char * argv[])
{
    try
    {
        int  result = 0;
        SherpaWeightMain SherpaWeight;
        
        SherpaWeightMain::RunParameters param;
        result = SherpaWeight.ParseCommandLine( argc, argv, param );
        if (result != 0)
            return result;
        
        return SherpaWeight.Run( param );
    }
    catch (const std::exception & error)
    {
        LogMsgError( "std exception caught in main(): \"%hs\"", FMT_HS(error.what()) );
        return -3;
    }
    catch (...)
    {
        LogMsgError("Unknown exception caught in main()");
        return -3;
    }
}
