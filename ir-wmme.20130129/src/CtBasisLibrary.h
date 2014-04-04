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

#ifndef CT_BASIS_LIBRARY_H
#define CT_BASIS_LIBRARY_H

#include <vector>
#include <string>
#include <iosfwd>

#include "CtBasisShell.h"

namespace ct {

struct FBasisSet;

struct FBasisSetLibraryImpl;

struct FBasisSetLibrary
{
public:
   // imports a .libmol file into memory. If given, a ``xxx loaded'' message
   // is prited to pxout.
   void ImportMolproLib(std::string const &FileName, std::ostream *pxout=0);

   // return false if failed. otherwise: *ADD*s basis functions to the
   // array, does not replace them.
   //   - nAtomIdx: index of atom in atom set allowed to be broken when not used
   //     anyway.
   //   - nElement: nuclear charge of atom.
   void LoadBasisFunctions(std::vector<FBasisShell> &Shells,
            int iElement, std::string const &BasisDesc,
            FVector3 const &vAtomPos, int iAtomIdx) const;


   FBasisSetLibrary();
   ~FBasisSetLibrary();
protected:
   FBasisSetLibraryImpl
      *p;
private:
   FBasisSetLibrary(FBasisSetLibrary const&); // not implemented
   void operator = (FBasisSetLibrary const&); // not implemented
};

extern FBasisSetLibrary
   g_BasisSetLibrary; // <- no point in having more than one of these things around.

} // namespace ct

#endif // CT_BASIS_LIBRARY_H
