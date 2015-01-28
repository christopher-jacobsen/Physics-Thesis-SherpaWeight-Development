//
//  main.cpp
//  SMEX

#include "common.h"
#include "SMEX.h"

int main(int argc, const char * argv[])
{
    try
    {
        int  result = 0;
        SMEX smex;
        
        SMEX::RunParameters param;
        result = smex.ParseCommandLine( argc, argv, param );
        if (result != 0)
            return result;
        
        return smex.Run( param );
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
