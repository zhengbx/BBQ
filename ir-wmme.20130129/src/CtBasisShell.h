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

#ifndef CT_BASIS_SHELL_H
#define CT_BASIS_SHELL_H

#include "CxTypes.h"
#include "CxPodArray.h"
#include "CxVec3.h"
#include "Ir.h" // only for making ir:FRawShell instances.

#include <iosfwd>

namespace ct {

typedef TVector3<double>
   FVector3;

// A generally contracted shell of basis functions of an atom:
//
// - A shell contains all contracted radial functions which are defined by the
//   fixed linear combination in CoMatrix
//
// - A shell contains all spherical components of the solid harmonic S(l,m).
//   E.g., for a d-shell, all five components d2- d1- d0 d1+ and d2+ are
//   included, for each radial function. (e.g., a d shell with two contractions
//   has ten functions).
//
// - This object DOES NOT contain a basis function center. These objects are
//   intended to be managed by a basis set library and shared across FBasisShell
//   objects (see below) of the same atom type. E.g., two C atoms sitting on
//   different centers can share the same FAtomShell object.
struct FAtomShell : public FIntrusivePtrDest {
   TArray<double>
      // one exponent for each primitive.
      Exponents,
      // nExp x nCo contraction matrix. Contractions are stored
      // with respect to unnormalized Gaussians.
      CoMatrix;
   unsigned
      AngMom;

   unsigned l() const { return AngMom; }
   unsigned nExp() const { return Exponents.size(); }
   unsigned nCo() const { return CoMatrix.size() / Exponents.size(); }
   unsigned nSh() const { return 2*l()+1; }
   unsigned nFn() const { return nCo() * nSh(); }
   double fCo(unsigned iExp, unsigned iCo) const { assert(iExp < nExp() && iCo < nCo()); return CoMatrix[iExp + nExp() * iCo]; }
   double fExp(unsigned iExp) const { assert(iExp < nExp()); return Exponents[iExp]; }

   enum FInitType {
      // if set, input co matrix is considered as referring to raw Gaussians
      // instead of normalized Gaussians.
      TYPE_Unnormalized = 0x01
   };

   // create an empty shell. Fill components yourself.
   FAtomShell() {};
   // create a generally contracted shell. By defauly pCo is a nExp x nCo matrix
   // given in library format (i.e., referring to normalized primitive Gaussians).
   FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_=0);
   // create a single normalized primitive.
   FAtomShell(unsigned AngMom_, double const fExp_, unsigned InitFlags_=0);
   // create a number of primitives.
   FAtomShell(unsigned AngMom_, double const *pExp_, unsigned nExp_, unsigned InitFlags_=0);
private:
   void Init(unsigned AngMom_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_, unsigned InitFlags_);
};

typedef boost::intrusive_ptr<FAtomShell>
   FAtomShellPtr;
typedef boost::intrusive_ptr<FAtomShell const>
   FAtomShellCptr;

// A shell of generally contracted molecular basis functions. Contrary to
// FAtomShell, this object is centered in space and represents an actual basis
// function used in the current calculation.
struct FBasisShell
{
   FAtomShellCptr
      pAs;
   FVector3
      vCenter;
   int
      // center index; not used by the low-level integral driver, but may be
      // used by the program or by mid-level functions to identify atoms.
      iCenter;

   FBasisShell() {};
   FBasisShell(FVector3 vCenter_, int iCenter_, FAtomShellCptr pAs_)
      : pAs(pAs_), vCenter(vCenter_), iCenter(iCenter_)
   {}

   ir::FRawShell MakeIrShell();

   unsigned l() const { return pAs->l(); }
   unsigned nExp() const { return pAs->nExp(); }
   unsigned nCo() const { return pAs->nCo(); }
   unsigned nSh() const { return pAs->nSh(); }
   unsigned nFn() const { return pAs->nFn(); }
   double fCo(unsigned iExp, unsigned iCo) const { return pAs->fCo(iExp, iCo); }
   double fExp(unsigned iExp) const { return pAs->fExp(iExp); }

   void PrintAligned(std::ostream &xout, uint Indent) const;
};

// calculate and return the integral Sqrt[<mu|mu>] for mu being a raw primitive Gaussian.
// Used for conversion of contraction coefficients between raw Gaussians (integral
// driver format) and normalized Gaussians (library format).
//
// This function accounts for the radial part. The angular normalization is
// included in the coefficients of Slm(r).
double RawGaussNorm(double fExp, unsigned l);


} // namespace ct

#endif // CT_BASIS_SHELL_H
