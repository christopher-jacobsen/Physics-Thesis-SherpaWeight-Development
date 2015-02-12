////////////////////////////////////////////////////////////////////////////////////////////////////
//  common.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef COMMON_H
#define COMMON_H

////////////////////////////////////////////////////////////////////////////////////////////////////
// common standard C includes

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cwchar>

////////////////////////////////////////////////////////////////////////////////////////////////////
// common standard C++ STL includes

#include <string>
#include <list>
#include <vector>
#include <array>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <stdexcept>
#include <system_error>

////////////////////////////////////////////////////////////////////////////////////////////////////
// message output methods

// TODO: add decl printf

inline void LogMsg( const char * prefix, const char * format, std::va_list args, std::FILE * fp ) throw()
{
    if (prefix) std::fputs( prefix, fp );
    std::vfprintf( fp, format, args );
    std::fputc( '\n', fp );
}

inline void LogMsgInfo( const char * format, ... ) throw()
{
    std::va_list args;
    va_start(args, format);
    LogMsg( nullptr, format, args, stdout );
    va_end(args);
}

inline void LogMsgWarning( const char * format, ... ) throw()
{
    std::va_list args;
    va_start(args, format);
    LogMsg( "Warning: ", format, args, stdout );
    va_end(args);
}

inline void LogMsgError( const char * format, ... ) throw()
{
    std::va_list args;
    va_start(args, format);
    LogMsg( "Error: ", format, args, stderr );
    va_end(args);
}

/////////////////////////////////////////////////////////////////////////////
// exception helpers

[[noreturn]] inline void ThrowError( const char * what )
{
    LogMsgError(what);
    throw std::runtime_error( what );
}

[[noreturn]] inline void ThrowError( std::errc code, const char * what )
{
    LogMsgError(what);
    throw std::system_error( make_error_code(code), what );
}

[[noreturn]] inline void ThrowError( const std::exception & error )
{
    LogMsgError(error.what());
    throw error;
}

/////////////////////////////////////////////////////////////////////////////
// string format helpers

#define FMT_I(x)    static_cast<int>(x)
#define FMT_HI(x)   static_cast<short>(x)
#define FMT_LI(x)   static_cast<long>(x)
#define FMT_LLI(x)  static_cast<long long>(x)

#define FMT_U(x)    static_cast<unsigned int>(x)
#define FMT_HU(x)   static_cast<unsigned short>(x)
#define FMT_LU(x)   static_cast<unsigned long>(x)
#define FMT_LLU(x)  static_cast<unsigned long long>(x)

#define FMT_F(x)    static_cast<double>(x)
#define FMT_LF(x)   static_cast<long double>(x)

#define FMT_HC(x)   static_cast<char>(x)
#define FMT_LC(x)   static_cast<wchar_t>(x)

#define FMT_HS(x)   static_cast<const char *>(x)
#define FMT_LS(x)   static_cast<const wchar_t *>(x)

#define FMT_P(x)    static_cast<const void *>(x)

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif  // COMMON_H
