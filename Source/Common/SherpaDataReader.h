////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaDataReader.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPA_DATA_READER_H
#define SHERPA_DATA_READER_H

#include "common.h"
#include <ATOOLS/Org/Data_Reader.H>

////////////////////////////////////////////////////////////////////////////////////////////////////

class DefaultDataReader : public ATOOLS::Data_Reader
{
public:
    DefaultDataReader()
        : Data_Reader(" ",";","!","=")
    {
        AddComment("#");
        AddWordSeparator("\t");
    }

    DefaultDataReader( const std::string & inputPath )
        : DefaultDataReader()
    {
        SetInputPath( inputPath );
    }

    DefaultDataReader( const std::string & inputPath, const std::string & inputFile )
        : DefaultDataReader( inputPath )
    {
        SetInputFile( inputFile );
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPA_DATA_READER_H
