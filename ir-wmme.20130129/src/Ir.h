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

#ifndef IR_H
#define IR_H

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

#ifdef _DEBUG
   #include <ostream> // for debug printing of shell information.
   #include <stdio.h> // for printf
#endif

#include "CxDefs.h"
#include "CxMemoryStack.h"
#include "IrBoysFn.h"

namespace ir {
   // low level shell data structure the integral drivers work with
   // (for simplifying talking to existing fortran and C++ programs!)
   struct FRawShell
   {
      unsigned
         // angular momentum, number of primitive exponents, number of contractions
         l, nExp, nCo;
      double const
         // pExp[iExp], iExp = 0..nExp-1: Exponent of primitive #i.
         *pExp,
         // pCo[iExp + nExp*iCo]: coefficient of primitive #iExp in contraction #iCo
         //
         // Note: Contraction coefficients are stored with respect to raw,
         // unnormalized Gaussians. Use the RawGaussNorm() function
         // to convert between contraction coefficients in terms of normalized
         // Gaussians (as used in libraries, for example), and raw coefficients
         // (as used in most programs)
         *pCo,
         // pCenter[i]: components 0..2: x,y,z position in space.
         *vCen;
      double const
         // Can be 0. If provided: Screening information.
         // [0]: largest range of any contracted function in the shell,
         // [1+iExp]: range of primitive #iExp
         // [1+nExp+iCo]: range of contraction #iCo
         *pRange;

      FRawShell() {}
      FRawShell(unsigned l_, double const *pExp_, unsigned nExp_, double const *pCo_, unsigned nCo_,
                double const *pvCen_, double const *pRange_=0)
         : l(l_), nExp(nExp_), nCo(nCo_), pExp(pExp_), pCo(pCo_),
           vCen(pvCen_), pRange(pRange_)
      {}

      // return number of spherical components per contraction
      inline unsigned nSh() const { return (2*l+1); }
      // return number of functions represented by the shell
      inline unsigned nFn() const { return (2*l+1) * nCo; }
      // return contraction coefficient of primitive iExp in contraction iCo.
      inline double fCo(unsigned iExp, unsigned iCo) const { assert(iCo < nCo && iExp < nExp); return pCo[iExp + nExp*iCo]; }

      inline double MaxCoRange() const { assert(pRange); return pRange[0]; }
      inline double CoRange(unsigned iCo) const { assert(pRange && iCo < nCo); return pRange[1+nExp+iCo]; }
      inline double ExpRange(unsigned iExp) const { assert(pRange && iExp < nExp); return pRange[1+nExp]; }
   };


   // integral kernels represent scalar functions G(r_12), to be used in matrix
   // elements like \int A(r_1) G(r_1 - r_2) C(r_2) d^3r_1 d^3r_2.
   struct FIntegralKernel
   {
      virtual double MaxRange() = 0;
      // evaluate Gm(rho,T) for m = 0...MaxM  (MaxM included), where T = rho |P-Q|^2.
      // All output data is multiplied by Factor.
      virtual void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const = 0;
      virtual ~FIntegralKernel() = 0;
   };

   // Kernel for G(r_12) = delta(r_12)
   struct FOverlapKernel : public FIntegralKernel
   {
      double MaxRange(); // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FOverlapKernel();
   };

   // Kernel for G(r_12) = 1/r_12
   struct FCoulombKernel : public FIntegralKernel
   {
      double MaxRange(); // override
      void EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const; // override
      ~FCoulombKernel();
   };

   extern FOverlapKernel g_IrOverlapKernel;
   extern FCoulombKernel g_IrCoulombKernel;


   void EvalInt2e3c(double *pOut, uint *Strides,
      FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, uint nC,
      double Prefactor, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   void EvalInt2e2c(double *pOut, uint StrideA, uint StrideC, FRawShell const *pA, FRawShell const *pC,
      double Prefactor, bool Add, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);
   void EvalInt2e2c_LaplaceC(double *pOut, uint StrideA, uint StrideC, FRawShell const *pA, FRawShell const *pC,
      double Prefactor, bool Add, unsigned LaplaceOrder, FIntegralKernel const *pKernel, ct::FMemoryStack &Mem);




#ifdef _DEBUG
void IrPrintMatrixGen(std::ostream &xout, double *pData, unsigned nRows, unsigned iRowSt, unsigned nCols, unsigned iColSt, std::string const &Caption);
void operator << (std::ostream &xout, FRawShell const &rs);
#endif // _DEBUG

} // namespace ir

#endif // IR_H
