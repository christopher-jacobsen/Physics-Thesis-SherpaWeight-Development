////////////////////////////////////////////////////////////////////////////////////////////////////
//  main.cpp
//  SherpaME
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaMEProgram.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    try
    {
        int result = 0;

        SherpaMEProgram program;
        
        SherpaMEProgram::RunParameters param;
        result = program.ParseCommandLine( argc, argv, param );
        if (result != 0)
            return result;
        
        return program.Run( param );
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
