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

#ifndef CT_IO_H
#define CT_IO_H

namespace ct {

extern char const
   *pResultFmt,
   *pResultFmtAnnoted,
   *pResultFmtI,
   *pResultFmtIAnnoted,
   *pTimingFmt;

} // namespace ct


#include <string>
#include "CxPodArray.h"

namespace ct {
   extern std::ostream
      &xerr, &xout;
   extern int
      Verbosity;
}

namespace ct {

// makes a string lower case. Will break with multibyte characters.
std::string tolower( std::string const &in );
// returns string 'in' with all spaces (but not tabs) removed.
std::string stripwhitespace( std::string const &in );

// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version="" );

// will load an entire file into the stringstream str. Returns false if failed.
// if 0 != FileLength, *FileLength will receive the number of loaded chars.
bool LoadFileIntoMemory( TArray<char> &pFileContent,
        std::string const &FileName, uint *pFileLength = 0 );
} // namespace ct.


#endif // CT_IO_H
