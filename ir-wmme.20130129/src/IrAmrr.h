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

/* IrAmrr.h v20130123 CET [storm, Gerald Knizia] */
#ifndef IR_RR_H
#define IR_RR_H

// IrAmrr -- Angular Momentum Recurrence Relations.

#include "CxDefs.h" // for assert and AIC_RP
#define IR_RP AIC_RP // restricted pointer

namespace ir {
   unsigned const
      MaxLa = 4,
      MaxLc = 5;

   // number of cartesian components with angular momentum <= l
   inline unsigned nCartX(int l) { return static_cast<unsigned>((l+1)*(l+2)*(l+3)/6); }
   // number of cartesians components with angular momentum == l
   inline unsigned nCartY(int l) { return static_cast<unsigned>((l+1)*(l+2)/2); }
   // number of solid harmonic components with angular momentum <= l
   inline unsigned nSlmX(int l) { return l*l; }
   // number of solid harmonic components with angular momentum == l
   inline unsigned nSlmY(int l) { return 2*l*1; }

   typedef unsigned short
      cart_index_t;

   void OsrrA(double *IR_RP pOut, double *IR_RP pGm, unsigned lab, double PmAx, double PmAy, double PmAz, double PmQx, double PmQy, double PmQz, double rho, double InvEta);
   void ShTrN(double *IR_RP pOut, double const *IR_RP pIn, unsigned N, unsigned l);
   void OsrrB_3c_shc(double *IR_RP pOut, double const *IR_RP pIn, double *IR_RP pMem, int la, unsigned lab, unsigned lc, double fPmQx, double fPmQy, double fPmQz, double InvEtaABC, double riz);
   void ShTrN_Indirect(double *IR_RP pOut, unsigned so, double const *IR_RP pIn, unsigned si, unsigned la, cart_index_t const *ii, unsigned N, unsigned M);
   void ShTrA_XY(double *IR_RP pOut, double const *IR_RP pIn, unsigned la, unsigned lab, unsigned M);
   void ShTrA_YY(double *IR_RP pOut, double const *IR_RP pIn, unsigned la, unsigned lab, unsigned M);
   void OsrrC(double *IR_RP pOut, unsigned sa, unsigned sb, double const *IR_RP pIn, double AmBx, double AmBy, double AmBz, unsigned lb, unsigned nCount);
   void ShellMdrr(double *IR_RP pOut, double const *IR_RP pIn, double Rx, double Ry, double Rz, unsigned lab);
   void ShellLaplace(double *IR_RP pOut, double const *IR_RP pIn, unsigned LaplaceOrder, unsigned lab);

   extern unsigned char iCartPow[56][3];
}


#endif // IR_RR_H
