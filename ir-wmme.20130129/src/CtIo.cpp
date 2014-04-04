/* Copyright (c) 2012  Gerald Knizia
 * 
 * This file is part of the IR/WMME program
 * (See http://www.theochem.uni-stuttgart.de/~knizia)
 * 
 * IR/WMME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * IR/WMME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
 */

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

// #include "CtCommon.h"
// #include "CtConstants.h"
// #define USE_INDENT_STREAM
#ifdef USE_INDENT_STREAM
   #define INDENTSTREAM_IMPL
   #include "CxIndentStream.h"
   namespace ct {
      static fmt::FIndentStream1
         s_OutputStreamAdp(std::cout,true,1);
      std::ostream
         &xout = s_OutputStreamAdp.stream,
         &xerr = std::cerr;
   } // namespace ct
#else
   namespace ct {
      std::ostream
         &xout = std::cout,
         &xerr = std::cerr;
   } // namespace ct
#endif // USE_INDENT_STREAM

#include "CtIo.h"


namespace ct {

char const
   *pResultFmt = " %-32s%18.12f\n",
   *pResultFmtAnnoted = " %-32s%18.12f  (%s)\n",
   *pResultFmtI = " %-32s%5i\n",
   *pResultFmtIAnnoted = " %-32s%5i  (%s)\n",
   *pTimingFmt = " Time for %s:%40t%10.2f sec\n";

int
   Verbosity = 0;

void FatalError( std::string const &Message,
    char const *pFromWhere, int nLine )
{
   std::stringstream str;
   str << Message;
   if ( pFromWhere ) {
      str << " at " << pFromWhere << ":" << nLine;
   };
   throw std::runtime_error(str.str());
};



// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version )
{
//     out << fmt::unind();
    out << "\n*** " << Name;
    if ( Version != "" )
        out << " [Ver. " << Version << "]";
    out << " ***\n" << std::endl;
//     out << fmt::eind();
};

// makes a string lower case. Will break with localized multibyte characters,
// but we don't such.
std::string tolower( std::string const &in )
{
    std::string
        Result(in);
    std::string::iterator
        it;
    for( it = Result.begin(); it != Result.end(); ++it )
        *it = std::tolower(*it);
    return Result;
};

// returns string 'in' with all spaces (but not tabs) removed.
std::string stripwhitespace( std::string const &in )
{
    std::string
        Result;
    Result.reserve( in.size() );
    for ( size_t i = 0; i < in.size(); ++i )
        if ( ' ' != in[i] )
            Result.push_back( in[i] );
    return Result;
};

// will load an entire file into the stringstream str. Returns false if failed.
bool LoadFileIntoMemory( TArray<char> &pFileContent,
        std::string const &FileName, uint *pFileLength )
{
    // read the entire file into an stringstream object.
    std::ifstream
        File( FileName.c_str() );
    std::size_t
        FileLength;
    if ( false == File.good() )
        return false;
    File.seekg( 0, std::ios::end );
    FileLength = File.tellg();
    if ( 0 != pFileLength )
        *pFileLength = FileLength;
    pFileContent.resize(2 + FileLength);
    memset( &pFileContent[0], 0, 2 + FileLength );
    File.seekg( 0, std::ios::beg );
    File.read( &pFileContent[0], FileLength );
    return true;
};


} // namespace ct
