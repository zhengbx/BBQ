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

#ifndef CT8K_BASISSET_H
#define CT8K_BASISSET_H

#include "Ir.h"
#include "CxPodArray.h"
#include "CtAtomSet.h"
#include "CtBasisShell.h"

namespace ct {

// IR driver-level interface to a basis set. Semi-lightweight.
// supposed to be generated from the program's native basis set
// format (FBasisSet here) on demand for interfacing between IR
// and the program.
struct FRawBasis : public FIntrusivePtrDest
{
   typedef TArray<ir::FRawShell>
      FShellArray;
   FShellArray
      Shells;
   TArray<unsigned>
      // offset of Shell[i] with respect to underlying basis in terms of
      // individual functions. Contains one additional element giving the
      // total number of functions.
      ShellOffsets;
   unsigned
      // largest number of functions in a shell occuring in the basis.
      // May be useful for allocating intermediates.
      nMaxFnPerShell;
   unsigned nFn() const { return ShellOffsets.back(); }
   unsigned nSh() const { return Shells.size(); }
   unsigned nFn(unsigned iSh) const { return ShellOffsets[iSh+1]-ShellOffsets[iSh]; }
   unsigned iFn(unsigned iSh) const { return ShellOffsets[iSh]; }
};
typedef boost::intrusive_ptr<FRawBasis>
   FRawBasisPtr;

// program native format of a basis set.
struct FBasisSet : FIntrusivePtrDest
{
   typedef std::vector<FBasisShell>
      FBasisShellArray;
   FBasisShellArray
      Shells;  // Gaussian shells centered on different atoms.
   FBasisContext
      Context; // for informative purposes only.
   std::string
      Name; // for informative purposes only.

   void Print(std::ostream &out) const;
   uint nFn() const;
   uint nPrimitiveGtos() const; // for decoration purposes only.

   FBasisSet(FAtomSet const &AtomSet, FBasisContext Context_);

public:
   FRawBasisPtr
      // Only valid if RebuildInterfacingData has been called after the basis
      // was last modified. If constructed regularly, this is done
      // automatically.
      pRawBasis;

   // rebuilds the data for low-level interfaces (i.e., the data in this section)
   // from this->Shells
   void RebuildInterfacingData();

   FRawBasisPtr MakeRawBasis();
private:
   // makes *this a basis set corresponding to the basis descriptions in
   // AtomSet. Attempts to load the neccessary external data. Throws
   // std::runtime_error if failed.
   void LoadFromAtomSet( FAtomSet const &AtomSet, FBasisContext Context );

   void MakeAtomOffsets( uint *&pAtomShellOffsets, uint *&pAtomBfnOffsets, uint nAtoms_, FMemoryStack &Mem ) const;
   void MakeShellOffsets( uint *&pShellOffsets, FMemoryStack &Mem ) const;
   void Finalize();
};

inline std::ostream &operator << ( std::ostream &out, FBasisSet const &BasisSet ) {
   BasisSet.Print(out);
   return out;
}

typedef boost::intrusive_ptr<FBasisSet>
   FBasisSetPtr;

struct FMatrixView;

// generalized 2-index integral kernel: effectively a functor which
// can calculate shell doublets for two shell objects. Used for making
// one-electron matrices.
struct FKrn2i
{
   virtual void EvalInt2e2c(double *pOut, uint StrideA, uint StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const = 0;
};

void MakeIntMatrix( FMatrixView &Out, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1., bool Add=false);
void MakeIntMatrix( FMatrixView &Out, FRawBasis const &RowBasis,
   FRawBasis const &ColBasis, FKrn2i const &Krn2i, FMemoryStack &Mem, double Prefactor=1., bool Add=false);

// evaluate <a|krn Laplace^LaplaceOrder|b> by directly calling
// a IR integral kernel.
struct FKrn2i_Direct : public FKrn2i
{
   FKrn2i_Direct(ir::FIntegralKernel const *pIrKernel_, uint LaplaceOrder=0);
   void EvalInt2e2c(double *pOut, uint StrideA, uint StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
protected:
   ir::FIntegralKernel const
      *m_pIrKernel;
   uint
      m_LaplaceOrder;
};

// evaluate <a| v(r) |b> where v(r) is sum of the fields emitted via pIrKernel
// for the sum of the supplied spherical point multipoles.
// Can be used, for example, to make nuclear attraction integrals.
struct FKrn2i_PointMultipoles : public FKrn2i
{
   // each point should provide 2*l+1 coefficients; one for each of its spherical
   // components. Slm order is equal to the order used for basis functions.
   // (in particular, dipole components are x,y,z, not x,z,y.).
   struct FPoint {
      FVector3 vPos;
      uint l;
   };
   FKrn2i_PointMultipoles(ir::FIntegralKernel const *pIrKernel_, FPoint const *pPoints, double const *pCoeffs, uint nPoints);
   void EvalInt2e2c(double *pOut, uint StrideA, uint StrideC, ir::FRawShell const *pA, ir::FRawShell const *pC, double Prefactor, bool Add, ct::FMemoryStack &Mem) const;
protected:
   ir::FIntegralKernel const
      *m_pIrKernel;
   TArray<ir::FRawShell>
      m_PointShells;
   TArray<double>
      m_Data;
   double
      m_Exp;
   uint
      m_MaxL;
};



} // namespace ct

#endif // CT8K_BASISSET_H
