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

#include <stdexcept>
#include <iostream>
#include <string.h> // for memcpy

#include <fstream> // for importing stuff (cp2k data)
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// #include "CtCommon.h"
#include "CtIo.h"
#include "CtAtomSet.h"
#include "CtBasisLibrary.h"
#include "CtBasisSet.h"
// #include "CtInt1e.h"
#include "CtMatrix.h"
// #include "AicDrv.h" // for making integrals

#include "CtTiming.h" // for embdedded HF timer.
#include "CtConstants.h"

#include "CxNumpyArray.h"
// #include "AicGridOps.h"

using boost::format;
using boost::str;

namespace ct {

void FormIntIJA(double *pIntFai, FMatrixView OrbA, FMatrixView OrbI, FMatrixView Jcd, FRawBasis const *pOrbBasis, FRawBasis const *pFitBasis, ir::FIntegralKernel &IntKernel, FMemoryStack &Mem)
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   uint
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn(),
      nOrbA = OrbA.nCols,
      nOcc = OrbI.nCols;
   if ( nAo != OrbI.nRows || nAo != OrbA.nRows )
      throw std::runtime_error("MakeFittingIntegrals: Input orbitals not consistent with orbital basis set.");
   double
      *pDF_NFi;
   Mem.Alloc(pDF_NFi, nAo * nFit * nOcc);

   for ( uint iShF = 0; iShF != pFitBasis->Shells.size(); ++ iShF ){
      ir::FRawShell const &ShF = pFitBasis->Shells[iShF];
      uint nFnF = ShF.nFn();

      double
         *pMNF;
      Mem.Alloc(pMNF, nAo * nAo * nFnF );

      for ( uint iShB = 0; iShB < pOrbBasis->Shells.size(); ++ iShB ){
         ir::FRawShell const &ShB = pOrbBasis->Shells[iShB];
         uint nFnB = ShB.nFn();
         for ( uint iShA = iShB; iShA < pOrbBasis->Shells.size(); ++ iShA ) {
            ir::FRawShell const &ShA = pOrbBasis->Shells[iShA];
            uint
               nFnA = ShA.nFn(),
               Strides[3] = {1, nFnA, nFnA * nFnB};
            double
               *pIntData;
            Mem.Alloc(pIntData, nFnA * nFnB * nFnF );

            ir::EvalInt2e3c(pIntData, Strides, &ShA, &ShB, &ShF,1, 1.0, &IntKernel, Mem);

            for ( uint iF = 0; iF < nFnF; ++ iF )
               for ( uint iB = 0; iB < nFnB; ++ iB )
                  for ( uint iA = 0; iA < nFnA; ++ iA ) {
                     double
                        f = pIntData[iA + nFnA * (iB + nFnB * iF)];
                     // assign to (\mu\nu| and (\nu\mu|. (int has perm symmetry).
                     pMNF[(pOrbBasis->iFn(iShA) + iA) + nAo * (pOrbBasis->iFn(iShB) + iB) + nAo * nAo * iF] = f;
                     pMNF[(pOrbBasis->iFn(iShB) + iB) + nAo * (pOrbBasis->iFn(iShA) + iA) + nAo * nAo * iF] = f;
                  }

            Mem.Free(pIntData);
         }
      }

      // D[\nu A i] = (i\nu|A) = C(\mu,i) (\mu [\nu|A)]
      FMatrixView
         DF_NFI(pDF_NFi + nAo * pFitBasis->iFn(iShF), nAo * nFnF, nOcc, 1, nFit * nAo);
      Mxm( DF_NFI,
           FMatrixView(pMNF, nAo * nFnF, nAo, nAo, 1),
           OrbI );

      Mem.Free(pMNF);
   }

   // Jcd version:
   //   Solve[Jcd[AB] D[\nu B i]] -> D[\nu A I]
   for ( uint i = 0; i < nOcc; ++ i ){
      FMatrixView
         D_NF(pDF_NFi + nAo * nFit * i, nAo, nFit),
         K_Fa(pIntFai + nFit * nOrbA * i, nFit, nOrbA);
//       TriangularSolve(Transpose(D_NF), Jcd);
//       Mxm(K_Fa, Transpose(D_NF), OrbA);
// ^- hm... I think I thought this was a good idea for some reason,
//          but I cannot remember why.
      Mxm(K_Fa, Transpose(D_NF), OrbA);
      TriangularSolve(K_Fa, Jcd);
   }

   Mem.Free(pDF_NFi);
   Mem.Free(pBaseOfMemory);
}




void MakeFittingIntegrals(double *pIntFai, FMatrixView OrbA, FMatrixView OrbI, FBasisSet const *pOrbBasis, FBasisSet const *pFitBasis, FMemoryStack &Mem)
{
   uint
      nAo = pOrbBasis->nFn(),
      nFit = pFitBasis->nFn();

   xout << boost::format(" ORBITAL BASIS %s:%35t%5i FUNCTIONS") % pOrbBasis->Name % nAo << "\n"
        << boost::format(" FITTING BASIS %s:%35t%5i FUNCTIONS") % pFitBasis->Name % nFit << "\n"
        << std::endl;

   // Make fitting coefficients J^{-1/2}.
   FTimer TimerJmh;
   FStackMatrix
      Jmh(nFit, nFit, &Mem);
   MakeIntMatrix(Jmh, *pFitBasis, *pFitBasis, FKrn2i_Direct(&ir::g_IrCoulombKernel), Mem);
   CalcCholeskyFactors(Jmh);
   xout << format(pTimingFmt) % "J^{-1/2} of fitting basis" % (double)TimerJmh; xout.flush();

   // make (F|ai) the lame way: assume we can keep (ai|A) in memory.
   FTimer TimerIJA;
   FormIntIJA(pIntFai, OrbA, OrbI, Jmh,
      &*pOrbBasis->pRawBasis, &*pFitBasis->pRawBasis, ir::g_IrCoulombKernel, Mem);
   xout << format(pTimingFmt) % "3-index integrals" % (double)TimerIJA; xout.flush();
};

#if 0
void MakeGridValues(double *pOut, FMatrixView Grid, FMatrixView Orb, uint GridDxOrder, FBasisSet const *pBasis, FMemoryStack &Mem)
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   uint
      nGridPt_Total = Grid.nCols,
      nComp = 1,
      nGridPt = 128,
      nAo = pBasis->nFn(),
      nOrb = Orb.nCols;
   if ( GridDxOrder == 1 ) nComp = 4;
   if ( GridDxOrder == 2 ) nComp = 10;
   std::cout << "nAo = " << nAo << "  Orb.nRows =" << Orb.nRows << std::endl;
   if ( nAo != Orb.nRows )
      throw std::runtime_error("MakeGridValues: Input orbitals not consistent with orbital basis set.");

   // convert basis format for low level driver input...
   FBasisSet::FGaussShellArray const
      &ShellsA = pBasis->Shells;
   FShellData
      *pShData;
   Mem.Alloc(pShData, ShellsA.size());
   for ( uint iSh = 0; iSh != ShellsA.size(); ++ iSh ){
      FGaussBfn const &BfA = *ShellsA[iSh].pFn;
      pShData[iSh] = FShellData(BfA.l, BfA.Contractions.size(), BfA.Exponents.size(), 0,
         &BfA.Exponents[0], &BfA.CoMatrix[0], ShellsA[iSh].vCenter);
   }

   double
      *pOrbValAo, *pOrbsCompressed;
   Mem.Alloc(pOrbValAo, nGridPt * nComp * nAo);
   Mem.Alloc(pOrbsCompressed, nAo * nOrb);
   assert(Grid.nRowSt == 1 && Grid.nRows == 3);

   for ( uint iGridPt = 0; iGridPt < nGridPt_Total; iGridPt += nGridPt ) {
      if ( iGridPt + nGridPt > nGridPt_Total )
         nGridPt = nGridPt_Total - iGridPt;
      void
         *pMem1 = Mem.Alloc(0);
      FORTINT
         *pMap, nMap = 0, iFnBase = 0;
      FShellData
         *pLastBf = pShData;
      Mem.Alloc(pMap, nAo);
      for ( FShellData *pFirstBf = pShData; pFirstBf != pShData + ShellsA.size(); pFirstBf = pLastBf ) {
         pLastBf = pFirstBf;
         uint nFn = 0;
         while ( pLastBf != pShData + ShellsA.size() && pLastBf->vCenter == pFirstBf->vCenter ) {
            nFn += pLastBf->nFn();
            ++pLastBf;
         }
         uint
            nBfStride = nGridPt*nComp;
         EvalShellGroupOnGrid(&pOrbValAo[nBfStride*nMap], nGridPt, nBfStride, // comp stride, bf stride
            pMap, nMap, iFnBase, &Grid(0,iGridPt), Grid.nColSt, nGridPt,
            pFirstBf, pLastBf, pFirstBf->vCenter,
            GridDxOrder, 1e-10, 40, Mem);
         iFnBase += nFn;
      }
      assert(iFnBase == nAo && nMap <= nAo);
      assert(pLastBf == pShData + ShellsA.size());

      FMatrixView
         OrbCmpr(pOrbsCompressed, nMap, nOrb),
         OrbValAo(pOrbValAo, nGridPt * nComp, nMap),
         Out(&pOut[iGridPt], nGridPt*nComp, nOrb, 1, nGridPt_Total*nComp);
      // compress input orbital matrix: delete basis functions not occurring in nMap.
      for ( uint iOrb = 0; iOrb != nOrb; ++ iOrb )
         for ( FORTINT iMap_ = 0; iMap_ != nMap; ++ iMap_ )
            OrbCmpr(iMap_, iOrb) = Orb(pMap[iMap_], iOrb);
      Mxm(Out, OrbValAo, OrbCmpr);
      Mem.Free(pMem1);
   };

   Mem.Free(pBaseOfMemory);
}
#endif



void WriteMatrixToFile_Rect(std::string const &FileName, std::string const Desc, double *pData, uint nRows, uint nCols, char const *pNumFmt)
{
    if ( pNumFmt == 0 )
        pNumFmt = "%21.14f";
    std::ofstream
        File(FileName.c_str());
    File << format("! %s, %i x %i\n") % Desc % nRows % nCols;
    if ( !File.good() )
        throw std::runtime_error("failed to open file '" + FileName + "' for writing.");
    for ( uint iRow = 0; iRow < nRows; ++ iRow ) {
        for ( uint iCol = 0; iCol < nCols; ++ iCol ) {
            File << " " << format(pNumFmt) % pData[iRow + nRows*iCol];
        };
        File << "\n";
    }
}


void WriteMatrixtoFile2(std::string const &FileName, ct::FMatrixView const &M, std::string const &MatrixName, std::string const &MatrixFormat, int Verbosity = 0)
{
   if ( MatrixFormat == "npy" ) {
      if ( M.nRowSt != 1 || M.nColSt != M.nRows )
         throw std::runtime_error("Sorry, WriteNpy can currently only export continuous matrices in default layout.");
      ct::WriteNpy(FileName, M.pData, ct::MakeShape(M.nRows, M.nCols));
   }
   else if ( MatrixFormat == "list" )
      WriteMatrixToFile(FileName, M, MatrixName);
   else if ( MatrixFormat == "rect" )
      WriteMatrixToFile_Rect(FileName, MatrixName, M.pData, M.nRows, M.nCols, " %22.15e");
   else
      throw std::runtime_error("output matrix format not recognized: '" + MatrixFormat + "'");
   if ( Verbosity >= 1 )
      ct::xout << format(" * wrote %s matrix to '%s'.") % MatrixName % FileName << std::endl;
};

// Set Out := L^T In R.
void BasisChange2( FMatrixView &Out, FMatrixView const &L,
    FMatrixView const &In, FMatrixView const &R, FMemoryStack &Mem)
{
    assert( L.nRows == In.nRows && R.nRows == In.nCols );
//     Out = FMatrixView( 0, L.nCols, R.nCols );
//     Mem.Alloc(Out.pData, Out.GetStridedSize());
    assert( Out.nRows == L.nCols && Out.nCols == R.nCols );

    if ( L.nRows * L.nCols <=  R.nRows * R.nCols ) {
        // T1 := L^T In
        FStackMatrix
            T1(L.nCols, In.nCols, &Mem);
        Mxm(T1, Transpose(L), In);
        Mxm(Out, T1, R);
    } else {
        // T1 := In * R
        FStackMatrix
            T1(In.nRows, R.nCols, &Mem);
        Mxm(T1, In, R);
        Mxm(Out, Transpose(L), T1);
    }

    if ( L.pData == R.pData && In.IsSymmetric(1e-10) )
       Symmetrize(Out);
}

} // namespace ct


int main_integral_export(int argc, char *argv[])
{
   using namespace ct;
   std::cout << " :::: IR/WMME [v20130129] -- Write Molecular Matrix Elements (Gerald Knizia, 2013) ::::"
             << std::endl;


   po::options_description options_desc("Options");
   std::vector< std::string >
      FileNames_BasisLibs;
   std::string
      BasisName_Orbital,
      BasisName_Orbital2,
      BasisName_Fit,
      FileName_Overlap,
      FileName_Overlap12,
      FileName_Overlap2,
      FileName_CoreH,
      FileName_2eFit,
      FileName_Atoms,
      FileName_AtomsAu,
      FileName_AtomsAng,
      FileName_InputOrbitals,
      FileName_GridCoords,
      FileName_GridValues,
      Name_OrbitalTrafo,
      MatrixFormat;
   int
      GridDxOrder = 0;
   options_desc.add_options()
      ("help",
          "print this help message")
      ("basis-lib", po::value<std::vector<std::string> >(&FileNames_BasisLibs)->composing(),
          "file names of .libmol basis set libraries to load (Molpro basis lib format). Can occur multiple times. "
          "all basis sets used in the input must be defined in one of those files.")
      ("save-overlap", po::value<std::string>(&FileName_Overlap)->default_value(""),
          "if given, save the AO overlap matrix S <mu|nu> to this file.")
      ("save-coreh", po::value<std::string>(&FileName_CoreH)->default_value(""),
          "if given, save the core Hamiltonian (1e terms) to this file.")
      ("save-fint2e", po::value<std::string>(&FileName_2eFit)->default_value(""),
          "if given, save the full (F|mu nu) 2e fitting integral matrix to this file. "
          "Fitting metric is already absorbed in F; i.e., the 2e integral "
          "(mu nu|eta xi) \\approx \\sum_F (mu nu|F)(F|eta xi).")
      ("basis-orb", po::value<std::string>(&BasisName_Orbital)->default_value("MINAO"),
          "name of the orbital basis to use for r,s. (e.g., cc-pVTZ, def2-SVP)")
      ("basis-fit", po::value<std::string>(&BasisName_Fit)->default_value("univ-JKFIT"),
          "name of the fitting basis to use for F (e.g., def2-TZVPP-JKFIT, def2-TZVPP-MP2FI)")
      ("orb-trafo", po::value<std::string>(&Name_OrbitalTrafo)->default_value("None"),
          "can be 'None' or 'Smh'. In the latter case, all output integrals (except "
          "for overlap matrices) are stored in terms of symmetrically orthogonalized "
          "basis functions for mu,nu instead of raw basis functions.")
      ("matrix-format", po::value<std::string>(&MatrixFormat)->default_value("list"),
          "Can be 'npy', 'list' or 'rect'. Defines file format of output matrices. "
          "npy is a binary format that can be loaded by numpy.load(..).")
      ("atoms-ang", po::value<std::string>(&FileName_AtomsAng)->default_value(""),
          "File name of input atoms (xyz format). Input coordinates are given in "
          "angstroms (1.0 == 10^-{10} meter). Normally .xyz files use these units.")
      ("atoms-au", po::value<std::string>(&FileName_AtomsAu)->default_value(""),
          "File name of input atoms (xyz format). Input coordinates are given in "
          "atomic units (1.0 == 1.0 bohr radius)")
      ("basis-orb-2", po::value<std::string>(&BasisName_Orbital2)->default_value(""),
          "name of a second orbital basis; if given, --save-overlap-2 and --save-overlap-12 "
          "can be used to calculate the overlap between two orbital basis sets.")
      ("save-overlap-2", po::value<std::string>(&FileName_Overlap2)->default_value(""),
          "if given, save the overlap matrix in the second orbital basis to this file")
      ("save-overlap-12", po::value<std::string>(&FileName_Overlap12)->default_value(""),
          "if given, save the overlap matrix between the first and the second basis to this file")
#if 0
      ("eval-orbitals", po::value<std::string>(&FileName_InputOrbitals)->default_value(""),
          "if given, evaluate given orbitals (nAo x nOrb matrix) on grid. nAo must be comaptible with basis-orb. If given, --grid-coords must also be given.")
      ("grid-coords", po::value<std::string>(&FileName_GridCoords)->default_value(""),
          "read grid coordinates(3 x nGrid matrix) from this file.")
      ("save-grid-values", po::value<std::string>(&FileName_GridValues)->default_value(""),
          "save values of orbitals on grid to this file.")
      ("eval-orbitals-dx", po::value<int>(&GridDxOrder)->default_value(0),
          "order of derivatives to evaluate for orbitals on grid: 0: values only ... 2: up to 2nd derivatives.")
#endif
   ;
//       Args =    [("--eval-orbitals-dx=%s" % DerivativeOrder)]
//       Inputs =  [("--grid-coords", "GRID", Grid)]
//       Outputs = [("--grid-values", "ORBS_ON_GRID")]

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, options_desc), vm);
   po::notify(vm);

   if (vm.count("help") || argc == 1) {
      xout << "\nThis program calculates molecular matrix elements (``integrals'') over Gaussian"
              "\nbasis functions and writes them to disk. See README for an explanation of the"
              "\nterms, the intended usage, and more information on input and output formats."
              "\n"
           << std::endl;
      xout << options_desc << std::endl;
      return 1;
   }

   // find file name of input atomic coordinates.
   double
      fInputAtomScale = 1.0;
   if ( FileName_AtomsAu != "" ) {
      FileName_Atoms = FileName_AtomsAu;
      fInputAtomScale = ToAng;
   } else if ( FileName_AtomsAng != "" ) {
      FileName_Atoms = FileName_AtomsAng;
      fInputAtomScale = 1.0;
   } else {
      throw std::runtime_error("either --atoms-au or --atoms-ang arguments must be given.");
   }

   // check combinations of other input options.
   if ( FileNames_BasisLibs.empty() )
      throw std::runtime_error("you must specify at least one basis set library file to load (--basis-lib missing)");
   if ( BasisName_Fit.empty() && !FileName_2eFit.empty() )
      throw std::runtime_error("if calculating 2e integrals, you must specify a fitting basis (--basis-fit)");
   bool
      TransformToSmh = false;
   if ( Name_OrbitalTrafo != "None" and Name_OrbitalTrafo != "Smh" )
      throw std::runtime_error("orbital transformation not recognized (--orb-trafo)");
   if ( Name_OrbitalTrafo == "Smh" )
      TransformToSmh = true;
   if ( MatrixFormat != "npy" && MatrixFormat != "list" && MatrixFormat != "rect" )
      throw std::runtime_error("Output matrix format must be 'npy', 'list' or 'rect' (--matrix-format)");
   if ( FileName_InputOrbitals == "" ) {
      if ( FileName_GridCoords != "" || FileName_GridValues != "" || GridDxOrder != 0 )
         throw std::runtime_error("no input orbitals provided. Grid options are to be used in conjunction with --eval-orbitals.");
   } else {
      if ( FileName_GridCoords == "" || FileName_GridValues == "" )
         throw std::runtime_error("when --eval-orbitals is given, --grid-coords and --save-grid-values must also be given.");
      if ( GridDxOrder != 0 && GridDxOrder != 1 && GridDxOrder != 2 )
         throw std::runtime_error("grid derivatives can only be evaluated up to 2nd order.");
   }
   if ( BasisName_Orbital2 == "" && (FileName_Overlap2 != "" || FileName_Overlap12 != "") )
      throw std::runtime_error("overlaps <2|2> or <1|2> between two orbital basis sets were requested, but no second orbital basis was specified. Use --basis-orb-2 to set one.");

//    xout << fmt::ind();

   MajorProgramIntro(xout, "BASIS LIBRARY");
   for ( uint iBasisLib = 0; iBasisLib < FileNames_BasisLibs.size(); ++ iBasisLib )
      g_BasisSetLibrary.ImportMolproLib(FileNames_BasisLibs[iBasisLib], &xout);

   boost::intrusive_ptr<FAtomSet>
      pAtoms = new FAtomSet;
   FBasisDescs
      DefaultBases;
   DefaultBases[BASIS_Orbital] = BasisName_Orbital;
   DefaultBases[BASIS_JkFit] = BasisName_Fit;
   if ( BasisName_Orbital2 != "" )
      DefaultBases[BASIS_Guess] = BasisName_Orbital2;
      // ^- slight abuse of the orbital context

   pAtoms->AddAtomsFromXyzFile(FileName_Atoms, DefaultBases);
   if ( fInputAtomScale != 1.0 ) {
      for ( uint iAtom = 0; iAtom < pAtoms->size(); ++ iAtom )
         (*pAtoms)[iAtom].vPos *= fInputAtomScale;
   }

   MajorProgramIntro(xout, "GEOMETRY & ENVIRONMENT");
   xout << " ATOMIC CONFIGURATION\n\n";
   xout << *pAtoms;

   MajorProgramIntro(xout, "EXPORT OF INTEGRALS");

   // instanciate the basis sets.
   FBasisSetPtr
      pOrbBasis = new FBasisSet(*pAtoms, BASIS_Orbital),
      pFitBasis;
   if ( FileName_2eFit != "" )
      pFitBasis = new FBasisSet(*pAtoms, BASIS_JkFit);
   uint
      nAo = pOrbBasis->nFn(),
      nFit = 0;
   if ( pFitBasis )
      nFit = pFitBasis->nFn();

   std::size_t
      RequiredMem = 4*nAo*nAo + (1500<<20);
   if ( FileName_2eFit != "" )
      RequiredMem += (2*nFit*nAo*nAo) + nFit*nFit;
   FMemoryStack2
      Mem(RequiredMem);

   // we'll just do the transformation always, and push a indentity matrix
   // if raw output integrals are requested. Can simply use MakeFittingIntegrals
   // then.
   FStackMatrix
      S(nAo, nAo, &Mem),
      Orbs(nAo, nAo, &Mem);
   pAtoms->MakeOverlapMatrix(S, *pOrbBasis, *pOrbBasis, Mem);
//    MakeIntMatrix(S, *pOrbBasis, *pOrbBasis, ir::g_IrOverlapKernel, Mem);
//    S.Print(xout, "OVERLAP");
   if ( !TransformToSmh ) {
      Orbs.SetIdentity();
   } else {
      Move(Orbs, S);
      CalcSmhMatrix(Orbs, Mem, FSmhOptions(1e-15,1e-15,(" S^{-1/2} for orbital basis ("+BasisName_Orbital+")").c_str(), &ct::xout));
//       S.SetIdentity();
   }
   if ( FileName_Overlap != "" )
      WriteMatrixtoFile2(FileName_Overlap, S, "Overlap", MatrixFormat, 1);

   if ( FileName_CoreH != "" ) {
      FStackMatrix
         CoreH_Ao(nAo, nAo, &Mem),
         CoreH_Mo(nAo, nAo, &Mem);
      pAtoms->MakeCoreHamiltonMatrix(CoreH_Ao, *pOrbBasis, *pOrbBasis, Mem);
      BasisChange2(CoreH_Mo, Orbs, CoreH_Ao, Orbs, Mem);
      WriteMatrixtoFile2(FileName_CoreH, CoreH_Mo, "CoreH", MatrixFormat, 1);
   }

   if ( FileName_2eFit != "" ) {
      FStackMatrix
         IntFai(nFit, nAo*nAo, &Mem);
      MakeFittingIntegrals(IntFai.pData, Orbs, Orbs, &*pOrbBasis, &*pFitBasis, Mem);
      WriteMatrixtoFile2(FileName_2eFit, IntFai, "2e-(F|rs)", MatrixFormat, 1);
   }

   if ( BasisName_Orbital2 != "" && (FileName_Overlap2 != "" || FileName_Overlap12 != "") ) {
      // make overlap/stuff matrices for overlap with second basis.
      FBasisSetPtr
         pOrbBasis2 = new FBasisSet(*pAtoms, BASIS_Guess);
      uint
         nAo2 = pOrbBasis2->nFn();
      if ( FileName_Overlap2 != "" ) {
         FStackMatrix
            S2(nAo2, nAo2, &Mem);
         pAtoms->MakeOverlapMatrix(S2, *pOrbBasis2, *pOrbBasis2, Mem);
         WriteMatrixtoFile2(FileName_Overlap2, S2, "Overlap/<basis2|basis2>", MatrixFormat, 1);
      }
      if ( FileName_Overlap12 != "" ) {
         FStackMatrix
            S12(nAo, nAo2, &Mem);
         pAtoms->MakeOverlapMatrix(S12, *pOrbBasis, *pOrbBasis2, Mem);
         WriteMatrixtoFile2(FileName_Overlap12, S12, "Overlap/<basis1|basis2>", MatrixFormat, 1);
      }
   }

#if 0
   if ( FileName_InputOrbitals != "" ) {
      FArrayNpy
         Grid, Orbs;
      ReadNpy(Grid, FileName_GridCoords);
      ReadNpy(Orbs, FileName_InputOrbitals);
      if ( Grid.Rank() != 2 || Grid.Shape[0] != 3 )
         throw std::runtime_error("input grid array must have rank 2 and shape (3,*).");
      if ( Orbs.Rank() != 2 || Orbs.Shape[0] != nAo )
         throw std::runtime_error("input orbital value array must have rank 2 and shape (nAo,*), where nAo is the number of functions in --basis-orb.");

      uint const
         nDerivativeComps_[] = {1,4,10},
         nGrid = Grid.Shape[1],
         nComps = nDerivativeComps_[GridDxOrder];
      TArray<double>
         GridValues(nAo * nComps * nGrid);
      FMatrixView
         Orb_(&Orbs.Data[0], Orbs.Shape[0], Orbs.Shape[1], Orbs.Strides[0], Orbs.Strides[1]),
         Grid1_(&Grid.Data[0], Grid.Shape[0], Grid.Shape[1], Grid.Strides[0], Grid.Strides[1]);
      FStackMatrix
         Grid_(Grid.Shape[0], Grid.Shape[1], &Mem);
      Move(Grid_, Grid1_); // <- copy over the array to get rid of funky strides.
      GridValues.clear_data();
      MakeGridValues(&GridValues[0], Grid_, Orb_, GridDxOrder, pOrbBasis.get(), Mem);
      FMatrixView
         GridValues_(&GridValues[0], nComps*Grid.Shape[1], Orbs.Shape[1]);
      WriteMatrixtoFile2(FileName_GridValues, GridValues_, "Orbitals on grid: (nGridPt x nDerivComp) x nOrb", MatrixFormat, 1);
   }
#endif
   xout << "\n";
   xout << format(pResultFmt) % "Nuclear repulsion energy" % pAtoms->NuclearRepulsionEnergy();

   xout.flush();
   return 0;
};


int main(int argc, char *argv[])
{
//    return mainWheHeehee();
   return main_integral_export(argc, argv);
};
