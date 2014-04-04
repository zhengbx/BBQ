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

#include <algorithm> // for std::swap
#include <cmath>
#include "CxVec3.h"
#include "IrAmrr.h"
#include "Ir.h"
// using ct::FMemoryStack;
using namespace ct;
using std::size_t;
typedef ct::TVector3<double>
   FVec3;

// #include <stdio.h>    // FIXME: remove this
// #include "CtMatrix.h" // FIXME: remove this

// references:
//   [1]: PCCP 8 3072 (2006), doi: 10.1039/b605188j
//   [2]: PCCP 6 5119 (2004), doi: 10.1039/b413539c

namespace ir {

template<class T>
inline T sqr(T x) {
   return x*x;
}

// return x^(3/2).
inline double pow15(double x) {
   return x * std::sqrt(x);
}

inline double DistSq3(double const *pA, double const *pB) {
   return sqr(pA[0] - pB[0]) + sqr(pA[1] - pB[1]) + sqr(pA[2] - pB[2]);
}

inline void SubVec3(double *pOut, double const *pA, double const *pB) {
   pOut[0] = pA[0] - pB[0];
   pOut[1] = pA[1] - pB[1];
   pOut[2] = pA[2] - pB[2];
}

// pOut += f * pIn
static void Add2(double *IR_RP pOut, double const *IR_RP pIn, double f, size_t n)
{
   size_t i = 0;
   for ( ; i < (n & ~3); i += 4 ) {
      pOut[i]   += f * pIn[i];
      pOut[i+1] += f * pIn[i+1];
      pOut[i+2] += f * pIn[i+2];
      pOut[i+3] += f * pIn[i+3];
   }
   for ( ; i != n; ++ i ) {
      pOut[i] += f * pIn[i];
   }
}



static bool IsWithinRange(FRawShell const *pA, FRawShell const *pB) {
   if (!pA->pRange || !pB->pRange)
      return true; // shells have no screening data.
   return ir::sqr(pA->MaxCoRange() + pB->MaxCoRange()) <= DistSq3(pA->vCen, pB->vCen);
}

static bool IsPrimitiveWithinRange(FRawShell const *pA, uint iExpA, FRawShell const *pB, uint iExpB, double fDistSqAB)
{
   if (!pA->pRange || !pB->pRange)
      return true; // shells have no screening data.
   return ir::sqr(pA->ExpRange(iExpA) + pB->ExpRange(iExpB)) <= fDistSqAB;
}

static bool IsPrimitiveWithinRange(FRawShell const *pA, uint iExpA, FRawShell const *pB, uint iExpB)
{
   return IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, DistSq3(pA->vCen, pB->vCen));
}


struct FGaussProduct
{
   double
      ExpA, ExpB, // exponents of the primitives
      Eta, // ExpA + ExpB
      InvEta, // 1/(ExpA+ExpB)
      Exp, // exponent of the primitive product.
      DistSq; // squared distance between A and B
   double
      vCen[3],
      vAmB[3];
   FGaussProduct(double const *vA, double ExpA_, double const *vB, double ExpB_) {
      ExpA = ExpA_; ExpB = ExpB_;
      Eta = ExpA + ExpB;
      InvEta = 1./Eta;
      Exp = ExpA * ExpB * InvEta;
      DistSq = 0;
      for ( uint i = 0; i < 3; ++ i ) {
         vCen[i] = InvEta * (ExpA * vA[i] + ExpB * vB[i]);
         vAmB[i] = vA[i] - vB[i];
         DistSq += sqr(vAmB[i]);
      }
   };

   double Sab() {
      return std::exp(-Exp * DistSq);
   }
};


static uint *MakeFnOffsets(FRawShell const *pCs, uint nC, FMemoryStack &Mem)
{
   uint *piFnC;
   Mem.Alloc(piFnC, nC+1);
   piFnC[0] = 0;
   for ( uint iC = 0; iC < nC; ++ iC )
      piFnC[iC+1] = piFnC[iC] + pCs[iC].nFn();
   return piFnC;
};

// accumulate matrix nSize  at pIn, corresponding to primitive iExpC,
// to nSize x nCo at pOut.
static void Contract1(double *pOut, double *pIn, uint nSize, FRawShell const *pC, uint iExpC)
{
   for ( uint iCoC = 0; iCoC < pC->nCo; ++ iCoC ) {
      double fCoC = pC->fCo(iExpC,iCoC);
      if (fCoC != 0)
         Add2(&pOut[nSize*iCoC], pIn, fCoC, nSize);
   }
}


void EvalInt2e3c(double *pOut, uint *Strides,
    FRawShell const *pA, FRawShell const *pB, FRawShell const *pCs, uint nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   uint
      StrideA = Strides[0], StrideB = Strides[1], StrideC = Strides[2];
   if (pA->l < pB->l) { // <- OsrrC only implemented for la >= lb.
      std::swap(pA, pB);
      std::swap(StrideA, StrideB);
   }

   // count number of C functions and find largest lc for memory purposes.
   uint
      *piFnC = MakeFnOffsets(pCs, nC, Mem),
      nFnC_Total = piFnC[nC],
      lc_Max = 0;
   for (uint iC = 0; iC < nC; ++ iC)
      lc_Max = std::max(lc_Max, (uint)pCs[iC].l);
   // allocate intermediates
   uint
      nCartX_AB = nCartX(pA->l + pB->l),
      nCartX_Am1 = nCartX(pA->l-1),
      nCartX_B = nCartX(pB->l),
      nCartX_ABmA = nCartX_AB - nCartX_Am1,
      nShA_CartXB = pA->nSh() * nCartX_B;
   // intermediates for primitive integrals
   double
      dbl = 0., // dummy.
      // Kernel derivatives (-d/dT)^m of (00|00) integral
      *pGm = Mem.AllocN(pA->l + pB->l + lc_Max + 1, dbl),
      // (a0|0) intermediate
      *p_A00 = Mem.AllocN(nCartX_AB, dbl),
      *p_A0C_sh_mem = Mem.AllocN(nCartX_ABmA * (2*lc_Max+1), dbl),
      *pMemOsrrB = Mem.AllocN(nCartX_AB * nCartX(lc_Max), dbl);
   // intermediates for contractions
   double
      // intermediates (a0|c) with AB primitives and C contracted, a = la..lab
      *p_A0C_ppc = Mem.AllocN(nCartX_ABmA * nFnC_Total, dbl),
      // intermediates (a0|c) with A,C contracted, a = la..lab.
      *p_A0C_cpc = Mem.AllocN(nCartX_ABmA * nFnC_Total * pA->nCo, dbl),
      // intermediates (a0|c) with A,B,C all contracted, a = la..lab.
      *p_A0C_ccc = Mem.ClearAllocN(nCartX_ABmA * nFnC_Total * pA->nCo * pB->nCo, dbl),
      // intermediates (xa|c) with A,B,C contracted and (xa| = nCartX(lb) x (2*la+1)
      *p_xAC_ccc = Mem.AllocN(nShA_CartXB * nFnC_Total * pA->nCo * pB->nCo, dbl);

//    printf("size on p_A0C_ccc = %i  expected: %i\n", p_xAC_ccc-p_A0C_ccc, nCartX_ABmA * nFnC_Total * pA->nCo * pB->nCo);
   FVec3
      vAmB = FVec3(pA->vCen) - FVec3(pB->vCen);
   double
//       fRangeKernel = sqr(pKernel->MaxRange()),
      fDistSqAB = LengthSq(vAmB);

   for (uint iExpB = 0; iExpB < pB->nExp; ++ iExpB)
   {
      memset(p_A0C_cpc, 0, nCartX_ABmA * nFnC_Total * pA->nCo * sizeof(*p_A0C_cpc));
      for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
      {
         // skip if Dist(A,B) < Range(A) + Range(B)
         if (!IsPrimitiveWithinRange(pA, iExpA, pB, iExpB, fDistSqAB))
            continue;

         FGaussProduct
            OvAB(pA->vCen, pA->pExp[iExpA], pB->vCen, pB->pExp[iExpB]);
            // ^- P == OvAB.vCen
         double
            Sab = std::exp(-OvAB.Exp * fDistSqAB), // [1] (6)
            PmA[3];
         SubVec3(PmA, OvAB.vCen, pA->vCen);

         memset(p_A0C_ppc, 0, nCartX_ABmA * nFnC_Total * sizeof(*p_A0C_ppc));
         for (uint iC = 0; iC < nC; ++ iC) {
            FRawShell const
               *pC = &pCs[iC];
            uint
               TotalL = pA->l + pB->l + pC->l;
            for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
            {
               FGaussProduct
                  OvPC(OvAB.vCen, OvAB.Eta, pC->vCen, pC->pExp[iExpC]);
               double
                  *PmC = OvPC.vAmB,
                  rho = OvPC.Exp, // [1] (3)
                  T = rho * OvPC.DistSq,
                  Factor = pow15(M_PI * OvPC.InvEta) * Sab * Prefactor; // [1] (7)

               // make I[m] = (00|0)^m, m = 0..TotalL (inclusive)
               pKernel->EvalGm(pGm, rho, T, TotalL, Factor);

               // make (a0|0)^m for a = 0..lab with lab + la+lb.
               OsrrA(p_A00, pGm + pC->l, (pA->l + pB->l), PmA[0], PmA[1], PmA[2],
                  PmC[0], PmC[1], PmC[2], rho, OvAB.InvEta);

               // make (a0|c) for a = la..lab, c = 0..lc.
               double
                  *p_A0C_sh;
               if (pC->l == 0) {
                  p_A0C_sh = p_A00 + nCartX_Am1;
               } else {
                  p_A0C_sh = p_A0C_sh_mem;
                  OsrrB_3c_shc(p_A0C_sh, p_A00, pMemOsrrB, pA->l, (pA->l + pB->l), pC->l,
                     PmC[0], PmC[1], PmC[2], OvPC.InvEta, rho/pC->pExp[iExpC]);
               }
               // (a0|c) with solid harmonic c is ready now. Just need to add it to
               // its contractions.
               Contract1(&p_A0C_ppc[nCartX_ABmA * piFnC[iC]], p_A0C_sh,
                  nCartX_ABmA*pC->nSh(), pC, iExpC);
            } // c exponents
         } // c shells
         // p_A0C_ppc should be done now. Contract A and B.
         Contract1(p_A0C_cpc, p_A0C_ppc, nCartX_ABmA * nFnC_Total, pA, iExpA);
      } // a exponents
      Contract1(p_A0C_ccc, p_A0C_cpc, (nCartX_ABmA * nFnC_Total * pA->nCo), pB, iExpB);
   } // b exponents

   // transform A to solid harmonics by factoring nCartX(lab) into nCartX(lb) x Slm(A).
   ShTrA_XY(p_xAC_ccc, p_A0C_ccc, pA->l, (pA->l + pB->l), nFnC_Total * pA->nCo * pB->nCo);

   // we now have nCartX(lb) x nShA x nFnC_Total x nCoA x nCoB at p_xAC_ccc.
   // we still need to move the angular momentum from a to b and to write the
   // output integrals to their final destination.
   for (uint iCoB = 0; iCoB < pB->nCo; ++ iCoB)
      for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA)
         for (uint iFnC = 0; iFnC < nFnC_Total; ++ iFnC) {
            uint
               iFnA = iCoA * pA->nSh(), iFnB = iCoB * pB->nSh();
            OsrrC(
               &pOut[iFnA*StrideA + iFnB*StrideB + iFnC*StrideC], StrideA, StrideB,
               &p_xAC_ccc[nShA_CartXB * (iFnC + (nFnC_Total * (iCoA + pA->nCo * iCoB)))],
               vAmB[0], vAmB[1], vAmB[2], pB->l, pA->nSh() );
         };

   Mem.Free(pBaseOfMemory);
}

/*
// evaluate (ab|c) for given range of shells [pA[0]..pA[nA]) x [pB[0]..pB[nB]) x [pC[0]..pC[nC])
// Output integral [ia,ib,ic] is stored at pOut[Strides[0]*ia + Strides[1]*ib + Strides[2]*ic],
// and integrals are multiplied by Prefactor.
void EvalInt2e3c(double *pOut, uint *Strides,
    FRawShell const *pAs, uint nA, FRawShell const *pBs, uint nB, FRawShell const *pCs, uint nC,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   for ( uint iA = 0; iA < nA; ++ iA ) {
      for ( uint iB = 0; iB < nB; ++ iB ) {
         // ...
      }
   }
};






void EvalInt2e3c(double *pOut, uint *Strides,
    FRawShell const &A, FRawShell const &B, FRawShell const &C,
    double Prefactor, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   return EvalInt2e3c(pOut, Strides, &A,1, &B,1, &C,1, pKernel, Mem);
};*/




// output: contracted kernels Fm(rho,T), format: (TotalL+1) x nCoA x nCoC
void Int2e2c_EvalCoKernels(double *pCoFmT, uint TotalL,
    FRawShell const *pA, FRawShell const *pC,
    double PrefactorExt, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   double
      t = DistSq3(pA->vCen, pC->vCen),
      *pFmT;
   Mem.Alloc(pFmT, TotalL + 1); // FmT for current primitive.

   // loop over primitives (that's all the per primitive stuff there is)
   for (uint iExpC = 0; iExpC < pC->nExp; ++ iExpC)
   for (uint iExpA = 0; iExpA < pA->nExp; ++ iExpA)
   {
      double
         Alpha = pA->pExp[iExpA],
         Gamma = pC->pExp[iExpC],
         InvEta = 1./(Alpha + Gamma),
         Rho = (Alpha * Gamma)*InvEta, // = (Alpha * Gamma)*/(Alpha + Gamma)
         Prefactor = (M_PI*InvEta)*std::sqrt(M_PI*InvEta); // = (M_PI/(Alpha+Gamma))^{3/2}

      Prefactor *= PrefactorExt;
      if(pC->l) Prefactor *= std::pow( 1.0/(2*Gamma), (int)pC->l); // <- use Hermites with D Ax := [1/(2 alpha)] \partial/\partial A_i.
      if(pA->l) Prefactor *= std::pow(-1.0/(2*Alpha), (int)pA->l); // <- -1 because \partial_A R \propto -\partial_B R!

      // calculate derivatives (D/Dt)^m exp(-rho t) with t = (A-C)^2.
      pKernel->EvalGm(pFmT, Rho, Rho*t, TotalL, Prefactor);

      // convert from Gm(rho,T) to Fm(rho,T) by absorbing powers of rho
      // (those would normally be present in the R of the MDRR)
      double
         RhoPow = 1.;
      for ( uint i = 0; i < TotalL + 1; ++ i ){
         pFmT[i] *= RhoPow;
         RhoPow *= 2*Rho;
      }

      // contract (lamely). However, normally either nCo
      // or nExp, or TotalL (or even all of them at the same time)
      // will be small, so I guess it's okay.
      for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC)
      for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
         double CoAC = pC->pCo[iExpC + pC->nExp*iCoC] *
                       pA->pCo[iExpA + pA->nExp*iCoA];
         Add2(&pCoFmT[(TotalL+1)*(iCoA + pA->nCo*iCoC)],
               pFmT, CoAC, (TotalL+1));
      }
   }

   Mem.Free(pFmT);
}

// write (2*la+1) x (2*lc+1) x nCoA x nCoC matrix to final destination.
static void Scatter2e2c(double *IR_RP pOut, unsigned StrideA, unsigned StrideC,
   double const *IR_RP pIn, unsigned la, unsigned lc, unsigned nCoA, unsigned nCoC, bool Add)
{
   unsigned nShA = 2*la+1, nShC = 2*lc+1;
   if ( Add ) {
      for (unsigned iCoC = 0; iCoC < nCoC; ++ iCoC)
         for (unsigned iCoA = 0; iCoA < nCoA; ++ iCoA)
            for (unsigned iShC = 0; iShC < nShC; ++ iShC)
               for (unsigned iShA = 0; iShA < nShA; ++ iShA)
                  pOut[(iShA + nShA*iCoA)*StrideA + (iShC + nShC*iCoC)*StrideC]
                     += pIn[iShA + nShA * (iShC + nShC * (iCoA + nCoA * iCoC))];
   } else {
      for (unsigned iCoC = 0; iCoC < nCoC; ++ iCoC)
         for (unsigned iCoA = 0; iCoA < nCoA; ++ iCoA)
            for (unsigned iShC = 0; iShC < nShC; ++ iShC)
               for (unsigned iShA = 0; iShA < nShA; ++ iShA)
                  pOut[(iShA + nShA*iCoA)*StrideA + (iShC + nShC*iCoC)*StrideC]
                      = pIn[iShA + nShA * (iShC + nShC * (iCoA + nCoA * iCoC))];
   }
}


// evaluate 2-electron 2-center integrals <a|krn|c>.
// if add is given: increment the output instead of overwriting it.
//
// evaluate 2-electron 2-center integrals <a|krn * laplace|c>
// note: to obtain the kinetic energy operator, pass an overlap kernel
//       and supply -.5 as Prefactor (ekin = -.5 laplace).
// if add is given: increment the output instead of overwriting it.
void EvalInt2e2c_LaplaceC( double *pOut, uint StrideA, uint StrideC,
    FRawShell const *pA, FRawShell const *pC, double Prefactor, bool Add,
    unsigned LaplaceOrder, FIntegralKernel const *pKernel, FMemoryStack &Mem )
{
   if (pA->l > pC->l) { // <- isn't this the wrong way around?
      std::swap(pA, pC);
      std::swap(StrideA, StrideC);
   }
   FVec3
      R;
   SubVec3(R, pA->vCen, pC->vCen);
   uint
      lc = pC->l, la = pA->l,
      TotalL = la + lc + 2*LaplaceOrder,
      TotalCo = pA->nCo * pC->nCo;
   double
      *pCoFmT, *pDataR_LapC, *pDataR, *pR1, *pFinal;
   Mem.ClearAlloc(pCoFmT, (TotalL+1) * TotalCo);
   Int2e2c_EvalCoKernels(pCoFmT, TotalL, pA, pC, Prefactor, pKernel, Mem);

   Mem.Alloc(pDataR_LapC, nCartY(TotalL));
   Mem.Alloc(pDataR, nCartY(la+lc) * TotalCo);
   Mem.Alloc(pR1, nCartY(la)*(2*lc+1) * TotalCo);
   Mem.Alloc(pFinal, (2*la+1)*(2*lc+1) * TotalCo);

   for (uint iCoC = 0; iCoC < pC->nCo; ++ iCoC)
   for (uint iCoA = 0; iCoA < pA->nCo; ++ iCoA) {
      // note: if skipping stuff here due to screening, the output must
      // be wiped unless Add == true!
      double
         *pFmT = &pCoFmT[(TotalL+1)*(iCoA + pA->nCo*iCoC)],
         *pDataR_ = &pDataR[nCartY(la+lc) * (iCoA + pA->nCo*iCoC)];
      if (LaplaceOrder == 0)
         ShellMdrr(pDataR_, pFmT, R[0], R[1], R[2], TotalL);
      else {
         ShellMdrr(pDataR_LapC, pFmT, R[0], R[1], R[2], TotalL);
         // note: a simple way of getting higher derivatives is to
         // just apply this function multiple times.
         ShellLaplace(pDataR_, pDataR_LapC, LaplaceOrder, la+lc);
      }
   }
   ShTrA_YY(pR1, pDataR, lc, (la + lc), TotalCo);
   ShTrA_YY(pFinal, pR1, la, la, (2*lc + 1)*TotalCo);
   // now: (2*la+1) x (2*lc+1) x nCoA x nCoC
   Scatter2e2c(pOut, StrideA, StrideC, pFinal, la, lc, pA->nCo, pC->nCo, Add);

   Mem.Free(pCoFmT);
};


void EvalInt2e2c( double *pOut, uint StrideA, uint StrideC,
    FRawShell const *pA, FRawShell const *pC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem )
{
   if (1) {
      return EvalInt2e2c_LaplaceC(pOut, StrideA, StrideC, pA, pC,
         Prefactor, Add, 0, pKernel, Mem);
   } else {
      assert(!Add);
      double ExpB = 0., CoB = 1., CenB[3] = {0., 0., 0.};
      FRawShell
         ShB(0, &ExpB, 1, &CoB, 1, &CenB[0], 0);
      uint Strides3[3] = {StrideA, 1, StrideC};
      return EvalInt2e3c(pOut, &Strides3[0], pA, &ShB, pC, 1, Prefactor, pKernel, Mem);
   }
}




} // namespace ir
