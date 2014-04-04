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

// This file is released under the GNU General Public License ("GPL", version 2)
// as part of the CT8K program. Copying, modification, creating derivative works
// and redistribution of the file is allowed, but _only_ subject to the terms
// of that GPL. You should have received a version of this license along
// with this source code. The program comes "as is", without any kind of
// warranty.
//
// Authors/Copyright holders:
//  - Gerald Knizia, 2006 (tag: cgk, contact: cgk.d@gmx.net)

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <utility> // for std::pair
#include <memory>
#include <set>
#include <map>
#include <list>
#include <set>
#include <sstream>
#include <iostream> // for cerr.
#include <boost/format.hpp>

using boost::str;
using boost::format;

// #include "CtCommon.h"
#include "CtBasisLibrary.h"
#include "CtIo.h"
#include "CtAtomSet.h" // for ElementNumberFromName

namespace ct {

typedef std::vector<double>
   FScalarArray;
typedef std::vector<int>
   FIntArray;
typedef std::list<std::string>
   FStringList;

typedef unsigned int
   uint;

struct FBasisSetLibraryImpl
{
protected:
   FBasisSetLibraryImpl
      *p;

   struct FBasisKey {
      std::string const
         *pName; // pointer into this->BasisNames.
      uint
         iElement; // element number
      bool operator < ( FBasisKey const &other ) const {
         if ( pName < other.pName ) return true;
         if ( pName > other.pName ) return false;
         return iElement < other.iElement;
      }
   };
   FBasisKey MakeKey(std::string const &Name, uint iElement, bool AssertExists = true) const;
   void TryMakeBasis(std::string const &BasisDesc, int iElement);

   typedef std::multimap<FBasisKey, FAtomShellPtr>
      FBasisFnMap;
   FBasisFnMap
      m_BasisFns;

   typedef std::set<std::string>
      FBasisNameSet;
   FBasisNameSet
      m_BasisNames; // list of all names of imported entries (AO sets, Aux sets, Guesses, ECPs...)
public:
   // imports a .libmol file into memory (does not yet write any files)
   void ImportMolproLib(std::string const &FileName, std::ostream *pxout);

   // return false if failed. otherwise: *ADD*s basis functions to the
   // array, does not replace them.
   //   - nAtomIdx: index of atom in atom set allowed to be broken when not used
   //     anyway.
   //   - nElement: nuclear charge of atom.
   void LoadBasisFunctions( std::vector<FBasisShell> &Shells,
            int iElement, std::string const &BasisDesc,
            FVector3 const &vAtomPos, int iAtomIdx ) const;
};


FBasisSetLibrary
   g_BasisSetLibrary;

FBasisSetLibraryImpl::FBasisKey FBasisSetLibraryImpl::MakeKey(std::string const &Name, uint iElement, bool AssertExists) const
{
   FBasisNameSet::const_iterator
      itName = m_BasisNames.find(Name);
   FBasisKey
      r;
   if ( itName == m_BasisNames.end() ) {
      if ( AssertExists ) {
         std::cerr << "Recognized names:\n";
         _for_each( itName, m_BasisNames )
            std::cerr << "   " << *itName << "\n";
         std::cerr << "\n";
         throw std::runtime_error("FBasisSetLibraryImpl: Basis/ECP entry '" + Name + "' not found.");
      }
      r.pName = 0;
      r.iElement = 0;
      return r;
   }
   r.pName = &*itName;
   r.iElement = iElement;
   return r;
};

void FBasisSetLibraryImpl::ImportMolproLib(std::string const &FileName, std::ostream *pxout)
{
   // how this files look:
   //   Note: Lines with * are comments.
      // He s STO-3G STO3G : 3 1 1.3
      // STO-3G
      // 6.3624214 1.158923 0.31364979 0.15432897 0.53532814 0.44463454
      // Li s STO-3G STO3G : 6 2 1.3 4.6
      // STO-3G
      // 16.119575 2.9362007 0.7946505 0.6362897 0.1478601 0.0480887 0.15432897
      // 0.53532814 0.44463454 -0.09996723 0.39951283 0.70011547
      // Li p STO-3G STO3G : 3 1 1.3
      // STO-3G
      // 0.6362897 0.1478601 0.0480887 0.15591627 0.60768372 0.39195739
      // B s cc-pVDZ VDZ : 9 3 1.9 1.9 9.9
      // cc-pVDZ
      // 4570 685.9 156.5 44.47 14.48 5.131 1.898 0.3329 0.1043 0.000696
      // 0.005353 0.027134 0.10138 0.272055 0.448403 0.290123 0.014322 -0.003486
      // -0.000139 -0.001097 -0.005444 -0.021916 -0.059751 -0.138732 -0.131482
      // 0.539526 0.580774 1
   // interpretation: First Line:
   //      ElementName OrbitalType N*[AlternativeNameI] :
   //            #PrimOrbitals #OrbitalsTheyAreContractedTo N*[ContractFrom.ContractTo]
   //      Any string, usually reference to be cited for the basis set
   //      M*[PrimitiveOrbitalExponents] L*[PrimitiveOrbitalCoefficients]
   //   for each contraction contraction length coefficients are supplied.
   // and data sets like this are also seen:
      // H  P 6-311G2P :   2  0
      // G92 6-311G polarization
      //    .15000000D+01   .37500000D+00
      // H  P 6-311G3P :   3  0
      // G92 6-311G polarization
      //    .30000000D+01   .75000000D+00   .18750000D+00
      // H  S 6-31G** 6-31G* 6-31G :   4  1      1.3
      // G92 6-31G**
       // .18731137D+02   .28253944D+01   .64012169D+00   .16127776D+00   .33494604D-01
      // .23472695D+00   .81375733D+00
   // all/some primitive orbitals are probably left uncontracted. Therefore no
   // contraction coefficients are supplied for these (only exponents).

   // Also I originally tought that the names of the basis set after the first
   // one were alternative Names of the set (like cc-pVDZ = VDZ). This is
   // however not how it works, at least not in all cases: It seems that names
   // which are supplied actually mean that the currently described elements
   // have to be inserted into the basis sets of ALL names provided, which
   // may or may not be identical basis sets (for example, for the s and p
   // orbitals of the 931G basis set, also the starred names are listed, which
   // means that these basis sets share these orbitals, altough in general they
   // are different).

   // so.. let's begin the mess.
   // read the entire file into a stringstream object (we need to modify it).
   TArray<char>
      pFileContent;
   if (!LoadFileIntoMemory(pFileContent, FileName))
      throw std::runtime_error( "FBasisSetLibraryImpl: Failed to open file '"+FileName+"' for basis library import." );

   // okay, now this is very bad. This FORTRAN stuff (sometimes) denotes
   // exponents in scientific notation not as "23423e-3" but as "23423D-03". We
   // do a lame attempt to convert that. This might break in some situations.
   // FIXME: correct this somehow.
   for (uint i = 0; i < pFileContent.size() - 2; ++i) {
      if (pFileContent[i]=='D' &&
          (pFileContent[i+1]=='+' || pFileContent[i+1]=='-') &&
          pFileContent[i+2]=='0')
         pFileContent[i] = 'E';
   }

   FBasisNameSet
      AllBasisNames;
   std::stringstream
      str(&pFileContent[0], std::stringstream::in);
   // ^- NOTE: "in" means from stream, into local data, not
   // into stream (file naming convention)

   try {
      while(str.good())
      {
         // clear exception mask, now stream will not throw() when something
         // unexpected happens.
         str.exceptions(std::ios::goodbit);
         std::string
            s;
         std::getline( str, s );
         if (s.size() == 0 || s[0] == '*' || s[0] == '!')
            // empty or comment line, throw it away and go on with the next
            continue;
         if (!str.good()) // eof, bad etc.
            break;
         str.exceptions(std::ios::failbit);
         // ^- when something fails, throw an exception. this will happen
         // if the actual file format does not match the one I had in mind
         // when coding this.

         // expected format: ElementName[w]OrbitalType[w]AlternativeNames[w] :
         // using namespace std;
         // cout << "** read line: '" << s << "'" << endl;
         std::stringstream
            line(s, std::stringstream::in);
         line.exceptions(std::ios::badbit | std::ios::failbit);
         std::string
            Element, Type;
         line >> Element >> Type;
         if ( Type.size() != 1 )
            throw std::runtime_error( "Parsing error, cannot interpret orbital type '" + Type + "'" );

         int
            AngMom;
         std::vector<std::pair<int,int> >
            Cos;
         std::vector<double>
            Exps;
         char
            cAngMom = ::tolower(Type[0]);
         for (AngMom = 0; AngMom < 9; ++ AngMom)
            if (cAngMom == "spdfghikl"[AngMom])
               break;
         if (AngMom == 9)
            throw std::runtime_error((format("Failed to understand angular momentum '%s'.") % cAngMom).str());
//          std::cout << "Element " << Element << " Type " << Type << std::endl;

         FStringList
            BasisNames; // all names of basis sets in which the
                     // current entry is to be inserted.
         for (line >> s; s != ":"; line >> s) {
//             std::cout << "Alternative Name:" << s << std::endl;
            BasisNames.push_back( tolower(stripwhitespace(s)) );
            AllBasisNames.insert(stripwhitespace(s));
         }

         // expected format: #prim orbitals #contractions (#contr.)*[a.b]
         // denoting indices of begin and end of a contraction with the
         // following exponents/contraction coefficients.
         int
            nExp,
            nCo,
            nCoeff(0), // total number of contraction coefficients to read (in all contractions).
            nHighestExpInCo(0); // 1-based index.
         line >> nExp >> nCo;
//          std::cout << "#Prim " << nExp << " #Co " << nCo << std::endl;
         Cos.reserve(nCo);
         for (int i = 0; i < nCo; ++ i){
            std::pair<int,int>
               iCo;
            char Dot;
            line >> iCo.first >> Dot >> iCo.second;
            iCo.first -= 1; // convert to 0-based [begin,end).
            if (Dot != '.')
               throw std::runtime_error("GTO-Contraction read format error.");
//             std::cout << "  Co: #" << iCo.first << "-#" << iCo.second << std::endl;
            if (iCo.second <= iCo.first || iCo.second > nExp)
               throw std::runtime_error("GTO-Contraction logical error.");
            nCoeff += iCo.second - iCo.first;
            nHighestExpInCo = std::max(nHighestExpInCo, iCo.second);
            Cos.push_back(iCo);
         }

         std::string
            EntryComment;
         do{ // read name, maybe skip comments.
            getline(str, EntryComment);
         } while (EntryComment.size() != 0 && EntryComment[0] == '*');
         // cout << "Entry Comment: " << EntryComment << endl;

         // now read exponents and contraction coefficients;
         // (this will break if comments are present in between)
         Exps.resize(nExp);
         std::vector<double>
            Coeffs(nCoeff, 0);
         for ( int i = 0; i < nExp; ++ i ){
            str >> Exps[i];
            // cout << "Exp: " << Exps.back() << endl;
         }

         // read in contraction coefficients.
         if ( nCo != 0 )
            for ( int i = 0; i < nCoeff; ++ i ){
               double Coeff;
               str >> Coeff;
               // cout << "Coeff: " << Coeff << endl;
               Coeffs[i] = Coeff;
            }
         // copy over the contraction coefficients to the contractions.
         std::vector<double>
            CoMatrix(Exps.size() * nCo, 0.);
         int
            iCoeff = 0;
         for ( int i = 0; i < nCo; ++ i ){
            for (int j = Cos[i].first; j != Cos[i].second; ++ j)
               CoMatrix[j + i*Exps.size()] = Coeffs[iCoeff + j - Cos[i].first];
            iCoeff += (Cos[i].second - Cos[i].first);
         }

         // in some files some primitive orbitals are left uncontracted.
         // but these are not stored as 1-GTO contractions but the
         // coefficients for these are just not present in the file.
         // Make 1-GTO contractions for them.
         int
            nAdditionalCo = nExp - nHighestExpInCo;
         if (0 != nAdditionalCo)
         {   // generate 1-GTO-each contractions manually.
            int nCoExplicit = nCo;
            nCo += nAdditionalCo;
            CoMatrix.resize(Exps.size() * nCo, 0.);
            int
               iStart = nHighestExpInCo;
               // ^- 0 based index, the rhs one is 1-based.
            for ( int i = 0; i < nAdditionalCo; ++ i ) {
               int iCo = i + nCoExplicit;
               CoMatrix[(i + iStart) + Exps.size() * iCo] = 1.;
            }
         }

         // import all names of the basis function
         FStringList::const_iterator
            itName;
         _for_each(itName, BasisNames)
            m_BasisNames.insert(*itName);

         // make the actual basis function and link it to all the names.
         FAtomShellPtr
            pBfn(new FAtomShell(AngMom, &Exps[0], Exps.size(), &CoMatrix[0],
                                CoMatrix.size()/Exps.size()));
         int
            iElement = ElementNumberFromName(Element);
         _for_each(itName, BasisNames)
            m_BasisFns.insert( FBasisFnMap::value_type(MakeKey(*itName, iElement), pBfn) );

         // chew the EOL marker if present, leave loop otherwise.
         str.exceptions(std::ios::goodbit);
         str.ignore(0xbad, '\n');
      };
   } catch (std::ios_base::failure &e){
      // this is not exactly something i would usually
      // call "error handling" but i really hate this string-
      // fiddling stuff and we can't really do anything better
      // about it anyway.
      std::cerr << "PARSER EXCEPTION:" << e.what() << std::endl;
      throw std::runtime_error( "Parsing of LIBMOL file FAILED because the actual syntax did not match the expected one. Last entry successfully processed: ");
   } catch (std::exception &e){
      std::cerr << "Exception during LibmolFile parsing: " << e.what() << std::endl;
      throw;
   };


   // if provided, write some imporant looking comments about what
   // we loaded to the standard output.
   if (pxout) {
      std::ostream
         &xout = *pxout;
      std::size_t
         iDirSep = FileName.rfind('/');
      if ( iDirSep == std::string::npos )
         iDirSep = 0;
      else
         iDirSep += 1;
      xout << format(" LOADED %-25s") % FileName.substr(iDirSep);
      if ( 1 ) {
         xout << "[";
         FBasisNameSet::const_iterator
            itSet;
         uint
            nLen = 0;
         _for_each(itSet, AllBasisNames) {
            if ( nLen >= 40 ) {
               xout << ",...";
               break;
            }

            if (itSet != AllBasisNames.begin())
               xout << ", ";
            xout << *itSet;
            nLen += itSet->size();
         }
         xout << "]";
      }
      xout << std::endl;
   }
};


// converts something like "sto3g>>spd" into "sto3g" and ['s','p','d'].
// First two parameters are Out parameters.
void ParseBasisDesc( std::string &BasisName, std::set<char> &OrbitalTypes,
   FBasisDesc const &BasisDesc )
{
   OrbitalTypes.clear();
   std::size_t
      nSeparatorIdx = BasisDesc.find(">>");
   // separator present in string?
   if ( nSeparatorIdx == std::string::npos ){
      // no->primitive declaration, "import all orbitals from [...]".
      BasisName = tolower(stripwhitespace(BasisDesc));
   } else { // yes. find out which orbitals to import.
      BasisName = tolower(stripwhitespace(BasisDesc.substr(0,nSeparatorIdx)));
      for ( std::size_t i = nSeparatorIdx+2; i < BasisDesc.size(); ++i ){
         char c = ::tolower( BasisDesc[i] );
         if ( ( c >= 'a' ) && ( c <= 'z' ) )
            // ^- we should consider using more b- and q-orbitals.
            OrbitalTypes.insert(c);
      }
   }
};

bool starts_with(std::string const &s, char const *p)
{
   std::size_t
      N = s.size();
   for ( std::size_t i = 0; i < N && *p && s[i] == *p; ++ i, ++p ){
   }
   return *p == 0;
};

int AngMomFromChar(char c) {
   switch (c) {
      case 's': return 0;
      case 'p': return 1;
      case 'd': return 2;
      case 'f': return 3;
      case 'g': return 4;
      case 'h': return 5;
      case 'i': return 6;
      case 'k': return 7;
      default: {
         std::stringstream str;
         str << "Angular momentum not recognized: '" << c << "'. Should be one of spdfghik.";
         throw std::runtime_error(str.str());
      }
   };
}

void FBasisSetLibraryImpl::TryMakeBasis(std::string const &BasisDesc, int iElement)
{
   if (starts_with(BasisDesc, ".et[")) {
      // syntax: .ET[spdf,<center>, <ratio>, <powmin>, <powmax>]
      char
         AngMoms[20];
      double
         fCenter, fRatio;
      int
         iPowMin, iPowMax, iScan;
      iScan = std::sscanf(BasisDesc.c_str(), ".et[%10[spdfghik],%lf,%lf,%i,%i]",
         &AngMoms[0], &fCenter, &fRatio, &iPowMax, &iPowMin);
      if (!((iScan == 4 && iPowMin > 0) || iScan == 5)) {
         std::stringstream str;
         str << "Failed to understand even tempered basis declaration '" << BasisDesc << "'."
            << " Expected: .et[<spd...>, <center>, <ratio>, <powmax>(, <powmin>)].";
         throw std::runtime_error(str.str());
      }
      if (iScan == 4)
         iPowMin = iPowMax;
      iPowMin *= -1;

      // make the basis function and store it.
      std::vector<double>
         Exps;
      Exps.reserve(iPowMax - iPowMin + 1);
      for (int i = iPowMax; i >= iPowMin; -- i)
         Exps.push_back(fCenter * std::pow(fRatio, (double)i));

      m_BasisNames.insert(BasisDesc);

      for (char *pAm = &AngMoms[0]; *pAm != 0; ++ pAm) {
         FAtomShellPtr
            pBfn(new FAtomShell(AngMomFromChar(*pAm), &Exps[0], Exps.size()));
         m_BasisFns.insert(FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn));
      }
   } else if ( starts_with(BasisDesc, "aug-" ) ||
            starts_with(BasisDesc, "daug-" ) ||
            starts_with(BasisDesc, "taug-" ) ) {
//       xout << "  in .aug generation code." << std::endl;
      uint
         iAugLevel = 1;
      if ( starts_with(BasisDesc, "daug-" ) )
         iAugLevel = 2;
      if ( starts_with(BasisDesc, "taug-" ) )
         iAugLevel = 3;
      // look up the parent basis.
      typedef FBasisFnMap::const_iterator
         FBfnIt;
      std::string
         ParentDesc(BasisDesc.substr((iAugLevel==1)? 4 : 5)); // includes the '-'.
      std::pair<FBfnIt,FBfnIt>
         itBfs = m_BasisFns.equal_range(MakeKey(ParentDesc, iElement));
      if ( itBfs.first == m_BasisFns.end() )
         // not there.
         return;

      m_BasisNames.insert(BasisDesc);

      uint
         AmMask = 0; // bitmask of angular momenta already handled.
      for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
         // note: this code assumes that all functions within a angular
         // momentum are stored inside one single contracted shell!
         FAtomShellCptr
            pFn = it->second;
         if ( (AmMask & (1 << pFn->AngMom)) != 0 ) {
            std::stringstream str;
            str << "sorry, diffuse augmentation code got confused. Encountered multiple shells with angular momentum l=" << pFn->AngMom << ".";
            throw std::runtime_error(str.str());
         }
         AmMask |= 1 << pFn->AngMom;

         // make exponents.
         TArray<double>
            Exps = pFn->Exponents;
         if ( Exps.empty() )
            throw std::runtime_error("Encountered a basis function without exponents during augmentation.");
         double
            fCen = Exps.back(),
            fRatio = 1./2.5;
         if ( Exps.size() > 1 ) {
            fRatio = Exps[Exps.size()-1]/Exps[Exps.size()-2];
         }
         for ( uint iAug = 0; iAug < iAugLevel; ++ iAug )
            Exps.push_back(fCen * std::pow(fRatio, (double)(1+iAug)));
         // make new contraction matrix.
         TArray<double>
            CoMatrix(Exps.size() * (pFn->nCo() + iAugLevel), 0.);
         // copy original data. Note that the number of exponents has changed,
         // so we need a strided copy.
         for (uint iCo = 0; iCo < pFn->nCo(); ++ iCo)
            for (uint iExp = 0; iExp < pFn->nExp(); ++ iExp)
               CoMatrix[iExp + iCo * Exps.size()] = pFn->CoMatrix[iExp + pFn->nExp() * iCo];
         // make coefficients for new primitives.
         for (uint iAug = 0; iAug < iAugLevel; ++ iAug)
            CoMatrix[Exps.size() * (pFn->nCo() + iAug) + iAug] =
               1./RawGaussNorm(Exps[pFn->nExp() + iAug + 1], pFn->AngMom);
         FAtomShellPtr
            pBfn(new FAtomShell(pFn->AngMom, &Exps[0], Exps.size(),
               &CoMatrix[0], CoMatrix.size()/Exps.size(), FAtomShell::TYPE_Unnormalized));
         m_BasisFns.insert(FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn));
      }
   }
};


void FBasisSetLibraryImpl::LoadBasisFunctions( std::vector<FBasisShell> &Shells,
   int iElement, std::string const &BasisDesc_,
   ct::FVector3 const &vAtomPos, int iAtomIdx ) const
{
   // todo: tokenize BasisDesc at ';'s here (and strip whitespace).
   // do the same as we he have below for the individual parts.

   std::string BasisDesc = BasisDesc_;
   if ( BasisDesc == "cc-pvtz-ig" ) BasisDesc = "cc-pvtz";
   if ( BasisDesc == "def2-tzvp-ig" ) BasisDesc = "def2-tzvp";
   if ( BasisDesc == "def2-sv(p)-jfit-ig" ) BasisDesc = "def2-sv(p)-jfit";

   typedef FBasisFnMap::const_iterator
      FBfnIt;
   FBasisKey const
      &Key = MakeKey(BasisDesc, iElement, false);
   std::pair<FBfnIt,FBfnIt>
      itBfs = m_BasisFns.equal_range(Key);
   if ( Key.pName == 0 || itBfs.first == itBfs.second ) {
      // basis not yet officially registered, but maybe it is one we
      // can make by modifying another one (e.g., aug-cc-pVTZ from cc-pVTZ)
      const_cast<FBasisSetLibraryImpl*>(this)->TryMakeBasis(BasisDesc, iElement);
      itBfs = m_BasisFns.equal_range(MakeKey(BasisDesc, iElement));
   }

   if ( itBfs.first == m_BasisFns.end() ) {
      std::stringstream str;
      str << "Requested basis '" << BasisDesc << "' for element '"
          << ElementNameFromNumber(iElement) << "' was not found.";
      throw std::runtime_error(str.str());
   }

   // note that we do NOT clear Shells but rather append only.
   for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
      FAtomShellCptr
         pFn = it->second;
      Shells.push_back(FBasisShell(vAtomPos, iAtomIdx, pFn));
   }
};


FBasisSetLibrary::FBasisSetLibrary()
   : p(new FBasisSetLibraryImpl())
{}

FBasisSetLibrary::~FBasisSetLibrary() {
   delete p;
}

void FBasisSetLibrary::ImportMolproLib(std::string const &FileName, std::ostream *pxout) {
   return p->ImportMolproLib(FileName, pxout);
}

void FBasisSetLibrary::LoadBasisFunctions(std::vector<FBasisShell> &Shells,
   int iElement, std::string const &BasisDesc,
   FVector3 const &vAtomPos, int iAtomIdx) const
{
   return p->LoadBasisFunctions(Shells, iElement, BasisDesc, vAtomPos, iAtomIdx);
}

} // namespace ct







// kate: indent-width 3; indent-mode normal;
