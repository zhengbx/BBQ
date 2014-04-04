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

#include <sstream> // for PrintAligned
#include <ostream> // for PrintAligned
#include <boost/format.hpp> // for PrintAligned

#include <math.h>
#include "Ir.h"
#include "CtBasisShell.h"

namespace ct {

static unsigned DoubleFactR(int l) {
   unsigned r = 1;
   while (l > 1) {
      r *= l;
      l -= 2;
   }
   return r;
}

double InvRawRawGaussNorm(double fExp, unsigned l)
{
   // purple book (6.6.14). (factor 4*pi/(2*l+1) converting between Ylm and Slm omitted)
   return pow(2*fExp/M_PI,.75) * sqrt(pow(4.*fExp,l) / DoubleFactR(2*l-1));
}

double RawGaussNorm(double fExp, unsigned l)
{
//    return 1./InvRawGaussNorm(fExp, l);
   return pow(M_PI/(2*fExp),.75) * sqrt(DoubleFactR(2*l-1)/pow(4.*fExp,l));
}


FAtomShell::FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_)
{
   Init(AngMom_, pExp_, nExp_, pCo_, nCo_, InitFlags_);
}

FAtomShell::FAtomShell(unsigned AngMom_, double const fExp_,  unsigned InitFlags_)
{
   double fCo = 1.;
   Init(AngMom_, &fExp_, 1, &fCo, 1, InitFlags_);
}

FAtomShell::FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, unsigned InitFlags_)
{
   TArray<double>
      CoMatrix_(nExp_ * nExp_, 0.);
   for (unsigned i = 0; i != nExp_; ++ i)
      CoMatrix_[(nExp_+1)*i] = 1.;
   Init(AngMom_, pExp_, nExp_, &CoMatrix[0], nExp_, InitFlags_);
}

void FAtomShell::Init(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_=0)
{
   AngMom = AngMom_;
   Exponents.resize(nExp_);
   CoMatrix.resize(nExp_ * nCo_);
   for (unsigned iExp = 0; iExp < nExp_; ++ iExp) {
      double
         fNorm = 1.;
      if (0 == (InitFlags_ & TYPE_Unnormalized))
         fNorm = 1./RawGaussNorm(pExp_[iExp], AngMom);
      Exponents[iExp] = pExp_[iExp];
      for (unsigned iCo = 0; iCo < nCo_; ++ iCo)
         CoMatrix[iExp + iCo*nExp_] = fNorm * pCo_[iExp + iCo*nExp_];
   }
}

ir::FRawShell FBasisShell::MakeIrShell()
{
   return ir::FRawShell(pAs->AngMom, &pAs->Exponents[0], pAs->nExp(), &pAs->CoMatrix[0], pAs->nCo(), &vCenter[0]);
};


void FBasisShell::PrintAligned(std::ostream &xout, uint Indent) const
{
   using boost::format;
   std::streampos
      p0 = xout.tellp(),
      p1;
//    xout << format("%3i " % iCenter;
   xout << format("%3i: %c   %8.4f %8.4f %8.4f    ")
          % nCo() % "spdfghiklm"[l()] % vCenter[0] % vCenter[1] % vCenter[2];
   p1 = xout.tellp();

   for (uint iExp = 0; iExp < pAs->Exponents.size(); ++iExp){
      if ( iExp != 0 ){
         xout << "\n";
         for ( uint i = 0; i < Indent + p1 - p0; ++ i )
            xout << " ";
      }
      xout << format("%16.7f  ") % pAs->Exponents[iExp];

      double
         fRenorm = RawGaussNorm(pAs->Exponents[iExp], l());
      std::stringstream
         str;
      for ( uint iCo = 0; iCo < nCo(); ++ iCo ){
         double
            fCo = pAs->CoMatrix[nExp() * iCo + iExp];
         if ( fCo != 0. )
            str << format(" %9.5f") % (fCo*fRenorm);
         else
            str << format(" %9s") % "  - - - -";
      }
      std::string
         s = str.str();
      while( !s.empty() && (s[s.size()-1] == ' ' || s[s.size()-1] == '-' ) )
         s.resize(s.size() - 1);
      xout << s;
   }
}



} // namespace ct
