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
#include <cmath>
#include "Ir.h"

#ifdef _DEBUG
   #include <boost/format.hpp>
   #include <iostream>
#endif


namespace ir {

// those things do not actually carry data, and they are very likely to be needed.
FOverlapKernel
   g_IrOverlapKernel;
FCoulombKernel
   g_IrCoulombKernel;


void FCoulombKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   IrBoysFn(pOut, T, MaxM, Factor*(2*M_PI)/rho);
}


double FCoulombKernel::MaxRange()
{ return 1e30; }


void FOverlapKernel::EvalGm(double *pOut, double rho, double T, unsigned MaxM, double Factor) const
{
   pOut[0] = Factor * std::exp(-T);
   for (unsigned i = 1; i <= MaxM; ++i)
      pOut[i] = pOut[i-1];
}

double FOverlapKernel::MaxRange()
{ return 0.; } // it's a delta function.


FIntegralKernel::~FIntegralKernel() {}
FCoulombKernel::~FCoulombKernel() {}
FOverlapKernel::~FOverlapKernel() {}


#ifdef _DEBUG
// some printing stuff which may be helpful in debug mode.
// Not included in default mode because it induces unnecessary dependencies, and
// the host programs probably have own means of printing this kind of stuff anyway.

template<class T>
void IrPrintRow(std::ostream &xout, T const *pData, int iColSt, unsigned nCols, std::string const &FmtF, char const *pEndLn="\n")
{
   for (unsigned iCol = 0; iCol < nCols; ++ iCol)
      if ( FmtF != "" )
         xout << boost::format("%14.6f") % pData[iColSt * iCol];
      else
         xout << boost::format(FmtF) % pData[iColSt * iCol];
   xout << pEndLn;
};

void IrPrintMatrixGen(std::ostream &xout, double *pData, unsigned nRows, unsigned iRowSt, unsigned nCols, unsigned iColSt, std::string const &Caption)
{
   using boost::format;
   xout << format("  Matrix %s, %i x %i.\n") % Caption % nRows % nCols;
   std::string
      FmtS = "%14s",
      FmtI = "%11i   ",
      FmtF = "%14.6f";
   xout << format(FmtS) % "";
   for (unsigned i = 0; i < nCols; ++ i)
      xout << " " << format(FmtI) % i;
   xout << std::endl;
   for (unsigned iRow = 0; iRow < nRows; ++ iRow) {
      xout << format(FmtI) % iRow;
      IrPrintRow(xout, &pData[iRow*iRowSt], iColSt, nCols, FmtF);
   }
   xout << std::endl;
};

void operator << ( std::ostream &xout, FRawShell const &rs )
{
   using boost::format;
   xout << format("IrRawShell: l = %i  cen = [%+8.3f,%+8.3f,%+8.3f]  nCo = %i  nExp = %i\n")
      % rs.l % rs.vCen[0] % rs.vCen[1] % rs.vCen[2] % rs.nCo % rs.nExp;
   xout << format(   "      exps:"); IrPrintRow(xout, rs.pExp, 1, rs.nExp, "%12.4e");
   if (rs.pRange) {
      xout << format("    rg/exp:"); IrPrintRow(xout, rs.pRange+1, 1, rs.nExp, "%12.4f");
   }
   for (unsigned iCo = 0; iCo < rs.nCo; ++ iCo) {
      xout << format("     %2i-co:") % iCo;
      IrPrintRow(xout, &rs.pCo[rs.nExp*iCo], rs.nExp, rs.nExp, "%12.6f", "");
      if (rs.pRange)
         xout << format(" | rg/co: %12.4f") % rs.CoRange(iCo);
      else
         xout << "\n";
   }
   xout.flush();
};

#endif // _DEBUG


} // namespace ir


#include <sstream>
void AicAssertFail( char const *pExpr, char const *pFile, int iLine )
{
   std::stringstream
      str;
   str << "assertion failed: '" << pExpr << "'";
   throw std::runtime_error(str.str());
};
