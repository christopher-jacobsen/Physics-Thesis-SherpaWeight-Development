////////////////////////////////////////////////////////////////////////////////////////////////////
//  main.cpp
//  SherpaME
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "common.h"
#include "SherpaMEProgram.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    try
    {
        int  result = 0;
        SherpaMEProgram program;
        
        SherpaMEProgram::RunParameters param;
        result = program.ParseCommandLine( argc, argv, param );
        if (result != 0)
            return result;
        
        return program.Run( param );
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
