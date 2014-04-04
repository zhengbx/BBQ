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

#include <cmath>
#include <stdexcept>
#include <boost/format.hpp>
using boost::format;
#include "Ir.h"
#include "CtAtomSet.h"
#include "CtBasisSet.h"
#include "CtBasisLibrary.h"
#include "CtMatrix.h"

#ifdef _DEBUG
   #include "CtIo.h"
#endif // _DEBUG

namespace ct {


FBasisSet::FBasisSet( FAtomSet const &AtomSet, FBasisContext Context_ )
   : Context(Context_)
{
   LoadFromAtomSet(AtomSet, Context);
}


std::string BasisContextName( FBasisContext Context )
{
   switch( Context ) {
      case BASIS_Orbital: return "ORBITAL";
      case BASIS_JFit: return "JFIT";
      case BASIS_JkFit: return "JKFIT";
      case BASIS_Mp2Fit: return "MP2FIT";
      case BASIS_CcsdFit: return "EXTFIT";
      case BASIS_F12RI: return "OPTRI";

      default:
         assert(0);
         return "[unknown basis context]";
   };
}


void FindDefaultBasis( std::string &Out, FBasisContext Context, std::string const &In )
{
   assert( Context != BASIS_Orbital );
   // TODO: do special element/context/basis-type specific stuff here.

   Out = In + "-" + BasisContextName(Context);
}


void FBasisSet::LoadFromAtomSet( FAtomSet const &AtomSet, FBasisContext Context )
{
   std::string
      BasisName;
   bool
      AllEqual = true;

   BasisName.reserve(64);
   for ( uint iAt = 0; iAt < AtomSet.Atoms.size(); ++ iAt ){
      FAtom const
         &Atom = AtomSet.Atoms[iAt];

      // get name of basis we are supposed to be loading.
      FBasisDescs::const_iterator
         itDesc = Atom.BasisDesc.find(Context);
      if ( itDesc != Atom.BasisDesc.end() )
         // basis for context explicitly set
         BasisName = itDesc->second;
      else {
         // basis for current context not set---go find one.
         if ( Context == BASIS_Guess ) {
            // should actually check for ECP-aware guess basis. But I guess
            // this will be okay for H-Kr.
            BasisName = "cc-pVTZ";
         } else {
            itDesc = Atom.BasisDesc.find(BASIS_Orbital);
            if ( itDesc == Atom.BasisDesc.end() ) {
               std::stringstream str;
               str << format("No basis set assigned to atom %i at (%.4f, %.4f, %.4f).") % iAt % Atom.vPos[0] % Atom.vPos[1] % Atom.vPos[2];
               throw std::runtime_error(str.str());
            }
            FindDefaultBasis( BasisName, Context, itDesc->second );
         }
      }

      if ( this->Name.empty() ) this->Name = BasisName;
      if ( this->Name != BasisName ) AllEqual = false;

      // we're importing molpro data... Molpro truncates basis names
      // at 32 characters. Do so here, too.
      BasisName.resize( std::min(BasisName.size(), 32ul) );

      // also, convert everything to lowercase. We do that when importing sets.
      for ( uint i = 0; i < BasisName.size(); ++ i )
         BasisName[i] = ::tolower(BasisName[i]);

      // ask basis library to add the functions.
      g_BasisSetLibrary.LoadBasisFunctions( this->Shells,
         Atom.AtomicNumber, BasisName,
         Atom.vPos, iAt );
   };

   if ( !AllEqual || this->Name.empty() )
      this->Name = BasisContextName(Context);

   // TODO:
   //   - prevent re-allocation of this->Shells during filling
   //   - make local copy of FGaussBfn objects and re-link this->Shells
   //     to these objects. This would reduce the memory spread of the data
   //     and prevent TLB misses during integration.

   Finalize();
}


void FBasisSet::Print( std::ostream &xout ) const
{
   xout << "Basis set '" << BasisContextName(Context) << "'\n" << std::endl;
   xout << "  Offs   NC/AM        Center/Position        Exponents           Contractions\n";
   xout << " -----------------------------------------   -----------------   ----------------------" << std::endl;
   uint
      nOff = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh ) {
      std::streampos p0 = xout.tellp(), p1;
      xout << format("%5i   ") % nOff;
      p1 = xout.tellp();
      Shells[iSh].PrintAligned(xout, p1-p0);
      xout << std::endl;
      nOff += Shells[iSh].nFn();
   }

   xout << std::endl;
};

uint FBasisSet::nPrimitiveGtos() const
{
   uint r = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      r += Shells[iSh].nExp() * Shells[iSh].nSh();
   return r;
};

uint FBasisSet::nFn() const
{
   uint r = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      r += Shells[iSh].nFn();
   return r;
};


void MakeIntMatrix( FMatrixView &Dest, FRawBasis const &BasisRow,
   FRawBasis const &BasisCol, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor, bool Add )
{
   bool
      MatrixSymmetric = (&BasisRow == &BasisCol);

   assert(Dest.nRows == BasisRow.nFn());
   assert(Dest.nCols == BasisCol.nFn());

   double
      *pIntResult;
   Mem.Alloc(pIntResult, BasisRow.nMaxFnPerShell * BasisCol.nMaxFnPerShell);
   for (uint iShellA = 0; iShellA < BasisRow.Shells.size(); ++ iShellA){
      ir::FRawShell const
         *pShellA = &BasisRow.Shells[iShellA];
      for (uint iShellB = (!MatrixSymmetric)? 0 : iShellA;
           iShellB < BasisCol.Shells.size(); ++ iShellB)
      {
         ir::FRawShell const
            *pShellB = &BasisCol.Shells[iShellB];
         uint
            nSizeA = pShellA->nFn(),
            nSizeB = pShellB->nFn();
         assert(nSizeA == BasisRow.nFn(iShellA) &&
                nSizeB == BasisCol.nFn(iShellB));
//          for (uint i = 0; i < BasisRow.nMaxFnPerShell * BasisCol.nMaxFnPerShell; ++i)
//             pIntResult[i] = 0.00777777; // remove this.
         Krn2i.EvalInt2e2c(pIntResult, 1, nSizeA, pShellA, pShellB, Prefactor, false, Mem);

         // fill data we gathered into the matrices
         for ( uint j_ = 0; j_ < nSizeB; ++ j_ )
            for ( uint i_ = 0; i_ < nSizeA; ++ i_ ) {
               uint
                  i = BasisRow.iFn(iShellA) + i_,
                  j = BasisRow.iFn(iShellB) + j_;
               FScalar const
                  &r = pIntResult[i_ + nSizeA * j_];
               if (!Add) {
                  Dest(i,j) = r;
                  if ( MatrixSymmetric )
                     Dest(j,i) = r;
               } else {
                  Dest(i,j) += r;
                  if ( MatrixSymmetric && iShellA != iShellB )
                     Dest(j,i) += r;
               }
            }
      }
   }
   Mem.Free(pIntResult);
}

void MakeIntMatrix( FMatrixView &Dest, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor, bool Add )
{
   return MakeIntMatrix(Dest, *RowBasis.pRawBasis, *ColBasis.pRawBasis, Krn2i, Mem, Prefactor, Add);
}

FKrn2i_Direct::FKrn2i_Direct(ir::FIntegralKernel const *pIrKernel_, uint LaplaceOrder)
   : m_pIrKernel(pIrKernel_), m_LaplaceOrder(LaplaceOrder)
{}

FKrn2i_PointMultipoles::FKrn2i_PointMultipoles(ir::FIntegralKernel const *pIrKernel_, FPoint const *pPoints, double const *pCoeffs, uint nPoints)
   : m_pIrKernel(pIrKernel_), m_Exp(1e20)
{
   m_PointShells.reserve(nPoints);
   // count number of coefficients and find max angular momentum.
   uint
      nCoeff = 0;
   m_MaxL = 0;
   for (uint iPt = 0; iPt < nPoints; ++ iPt) {
      nCoeff += 2*pPoints[iPt].l + 1;
      m_MaxL = std::max(m_MaxL, pPoints[iPt].l);
   }
   // store the point centers and their strenghts here.
   m_Data.resize(3*nPoints + nCoeff);
   // copy centers
   for (uint iPt = 0; iPt < nPoints; ++ iPt)
      for (uint i = 0; i < 3; ++ i)
         m_Data[3*iPt + i] = pPoints[iPt].vPos[i];
   // copy coefficients and make RawShell objects.
   std::size_t
      iCoeff = 0;
   for (uint iPt = 0; iPt < nPoints; ++ iPt) {
      FPoint const &Pt = pPoints[iPt];
      m_PointShells.push_back(ir::FRawShell(Pt.l, &m_Exp, 1,
         &m_Data[3*nPoints+iCoeff], 1, &m_Data[3*iPt]));

      double
         fNorm = RawGaussNorm(m_Exp, Pt.l);
//       xout << format("iPt = %i  m_Exp = %10.4e  fNrm = %10.4e  fChg = %.4f") %
//                iPt % m_Exp % fNorm % pCoeffs[iCoeff] << std::endl;
      // convert from norm to self-overlap. We want the factor which
      // would make <c|c> = 1, and it occurs on both sides.
      // What the sqrt(8) is for.. I am very confused about it.
      // It works for point charges, but might be wrong for higher multipoles.
      fNorm = fNorm * fNorm * sqrt(8.);
      if (Pt.l != 0)
         throw std::runtime_error("case Pt.l != 0 not tested. Norm factor likely wrong!.");
      for (uint iSh = 0; iSh < 2*Pt.l + 1; ++ iSh)
         m_Data[3*nPoints + iCoeff + iSh] = pCoeffs[iCoeff + iSh] / fNorm;
      iCoeff += 2*Pt.l+1;
   };
   assert(iCoeff == nCoeff);
}

void FKrn2i_PointMultipoles::EvalInt2e2c(double *pOut, uint StrideA, uint StrideB, ir::FRawShell const *pA, ir::FRawShell const *pB, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
   if (!Add)
      for (uint iB = 0; iB != pB->nFn(); ++ iB)
         for (uint iA = 0; iA != pA->nFn(); ++ iA)
            pOut[StrideA*iA + StrideB*iB] = 0.;
   double
      *pIntData = Mem.AllocN(pA->nFn() * pB->nFn() * (2*m_MaxL+1), (double)0.);
   for (uint iPt = 0; iPt < m_PointShells.size(); ++ iPt) {
      uint
         Strides3[3] = {1, pA->nFn(), pA->nFn()*pB->nFn()};
      // todo: do this with a specialized routine which takes ZetaC == +inf
      // into account. Can be done much cheaper then.
      ir::FRawShell const
         *pShPt = &m_PointShells[iPt];
      EvalInt2e3c(pIntData, Strides3,
         pA, pB, pShPt, 1, Prefactor, m_pIrKernel, Mem);
      for (uint iB = 0; iB != pB->nFn(); ++ iB)
         for (uint iA = 0; iA != pA->nFn(); ++ iA)
            for (uint iCo = 0; iCo < pShPt->nSh(); ++ iCo)
               pOut[iA*StrideA + iB*StrideB] += pIntData[iA*Strides3[0] + iB*Strides3[1] + iCo*Strides3[2]];
   }
   Mem.Free(pIntData);
};


void FKrn2i_Direct::EvalInt2e2c(double *pOut, uint StrideA, uint StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const
{
   ir::EvalInt2e2c_LaplaceC(pOut, StrideA, StrideC, pA, pC, Prefactor, Add,
      m_LaplaceOrder, m_pIrKernel, Mem);
}


void FBasisSet::MakeAtomOffsets( uint *&pAtomShellOffsets, uint *&pAtomBfnOffsets, uint nAtoms_, FMemoryStack &Mem ) const
{
   // pAtomShellOffsets: maps atom id to first basis function shell
   //    (Shells[i]) of the atom
   // pAtomBfnOffsets: maps atom id to first basis function (AO index)
   //    of the atom
   uint
      nAtoms = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      if ( Shells[iSh].iCenter > 0 )
         nAtoms = std::max(nAtoms, 1 + static_cast<unsigned>(Shells[iSh].iCenter));
   assert( nAtoms <= nAtoms_ );
   // ^- might differ if there are some atoms without basis functions
   nAtoms = nAtoms_;

   Mem.Alloc(pAtomShellOffsets, nAtoms+1);
   Mem.Alloc(pAtomBfnOffsets, nAtoms+1);
   *pAtomShellOffsets = 0;
   *pAtomBfnOffsets = 0;
   // note: this code assumes that this->shells is ordered according
   // to the atom set!
   for ( uint iAt = 0; iAt < nAtoms; ++ iAt ){
      uint
         iSh = pAtomShellOffsets[iAt],
         iBf = pAtomBfnOffsets[iAt];
      while ( iSh < Shells.size() && Shells[iSh].iCenter == (signed)iAt ) {
         iBf += Shells[iSh].nFn();
         ++ iSh;
      }
//       _xout0("iAt: " << iAt << "   iSh: "<< iSh << "   iBf: " << iBf);
      pAtomShellOffsets[iAt+1] = iSh;
      pAtomBfnOffsets[iAt+1] = iBf;
      assert(iBf <= nFn());
   };
}

void FBasisSet::MakeShellOffsets( uint *&pShellOffsets, FMemoryStack &Mem ) const
{
   Mem.Alloc(pShellOffsets, Shells.size() + 1);
   *pShellOffsets = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      pShellOffsets[iSh+1] = pShellOffsets[iSh] + Shells[iSh].nFn();
};


void FBasisSet::RebuildInterfacingData()
{
   pRawBasis = MakeRawBasis();
};

FRawBasisPtr FBasisSet::MakeRawBasis()
{
   FRawBasisPtr
      r(new FRawBasis());
   r->Shells.reserve(Shells.size());
   r->ShellOffsets.reserve(Shells.size() + 1);
   r->ShellOffsets.push_back(0);
   r->nMaxFnPerShell = 0;
   for (unsigned iSh = 0; iSh < Shells.size(); ++ iSh){
      r->Shells.push_back(Shells[iSh].MakeIrShell()); // three copies...
      unsigned
         nShFn = Shells[iSh].nFn();
      r->ShellOffsets.push_back(r->ShellOffsets.back() + nShFn);
      r->nMaxFnPerShell = std::max((unsigned)r->nMaxFnPerShell, nShFn);
   }
   assert(r->ShellOffsets.back() == this->nFn());
   return r;
};


void FBasisSet::Finalize()
{
   RebuildInterfacingData();
};




} // namespace ct
