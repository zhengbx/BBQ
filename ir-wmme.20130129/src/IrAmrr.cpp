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

/* IrAmrr.cpp v20130123 CET [storm, Gerald Knizia] */
#include <stddef.h> // for size_t
#include "IrAmrr.h" // for cart_index_t

namespace ir {

// consntants for direct expansion of solid harmonics in terms of cartesians
// these ones are good for L=0..5
static const double sd0 = 5.e-01;
static const double sd1 = 1.7320508075688772;
static const double sd10 = 7.3950997288745202e-01;
static const double sd11 = 4.4370598373247123;
static const double sd12 = 5.5901699437494734e-01;
static const double sd13 = 3.3541019662496838;
static const double sd14 = 2.9580398915498081;
static const double sd15 = 2.0916500663351889;
static const double sd16 = 6.2749501990055672;
static const double sd17 = 4.8412291827592718e-01;
static const double sd18 = 9.6824583655185437e-01;
static const double sd19 = 5.809475019311126;
static const double sd1a = 2.5617376914898995;
static const double sd1b = 5.1234753829797981;
static const double sd1c = 5.2291251658379723e-01;
static const double sd1d = 1.0458250331675947;
static const double sd1e = 4.1833001326703778;
static const double sd1f = 1.5687375497513918;
static const double sd2 = 8.660254037844386e-01;
static const double sd20 = 1.2549900398011134e+01;
static const double sd21 = 8.8741196746494246;
static const double sd22 = 2.2185299186623562;
static const double sd23 = 1.3311179511974137e+01;
static const double sd24 = 3.5078038001005702;
static const double sd25 = 7.0156076002011396;
static const double sd26 = 7.0156076002011403e-01;
static const double sd27 = 1.8750000000000002;
static const double sd28 = 3.7500000000000004;
static const double sd29 = 5.;
static const double sd2a = 1.0246950765959596e+01;
static const double sd3 = 6.1237243569579458e-01;
static const double sd4 = 2.4494897427831779;
static const double sd5 = 1.5;
static const double sd6 = 7.9056941504209488e-01;
static const double sd7 = 2.3717082451262845;
static const double sd8 = 3.872983346207417;
static const double sd9 = 1.9364916731037085;
static const double sda = 3.75e-01;
static const double sdb = 7.5e-01;
static const double sdc = 3.;
static const double sdd = 1.1180339887498947;
static const double sde = 6.7082039324993676;
static const double sdf = 3.1622776601683791;

// Calculate (a0|0)^0 from [0]^m, for a = 0..lab, m = 0..lab (both inclusive).
// PmA is P - A, PmQ is P - Q, InvEta is 1/Eta. See PCCP 6 5119 (2004) eq. 11.
// Output is ordered as AngCompsFull(lab).
void OsrrA(double *IR_RP pOut, double *IR_RP pGm, unsigned lab, double PmAx, double PmAy, double PmAz, double PmQx, double PmQy, double PmQz, double rho, double InvEta)
{
   // m + lab == 0
   pOut[0] = pGm[0];
   if (lab == 0) return;

   // m + lab == 1
   double rie = -rho*InvEta;
   double rPmQx = rie*PmQx, rPmQy = rie*PmQy, rPmQz = rie*PmQz;
   pOut[1] = PmAx*pGm[0] + pGm[1]*rPmQx;
   pOut[2] = PmAy*pGm[0] + pGm[1]*rPmQy;
   pOut[3] = PmAz*pGm[0] + pGm[1]*rPmQz;
   if (lab == 1) return;

   // m + lab == 2
   double iz2 = .5*InvEta;
   double r_100_1 = PmAx*pGm[1] + pGm[2]*rPmQx;
   double r_010_1 = PmAy*pGm[1] + pGm[2]*rPmQy;
   double r_001_1 = PmAz*pGm[1] + pGm[2]*rPmQz;
   pOut[4] = PmAx*pOut[1] + iz2*(pGm[0] + pGm[1]*rie) + rPmQx*r_100_1;
   pOut[5] = PmAy*pOut[2] + iz2*(pGm[0] + pGm[1]*rie) + rPmQy*r_010_1;
   pOut[6] = PmAz*pOut[3] + iz2*(pGm[0] + pGm[1]*rie) + rPmQz*r_001_1;
   pOut[7] = PmAy*pOut[1] + rPmQy*r_100_1;
   pOut[8] = PmAz*pOut[1] + rPmQz*r_100_1;
   pOut[9] = PmAz*pOut[2] + rPmQz*r_010_1;
   if (lab == 2) return;

   // m + lab == 3
   double r_100_2 = PmAx*pGm[2] + pGm[3]*rPmQx;
   double r_010_2 = PmAy*pGm[2] + pGm[3]*rPmQy;
   double r_001_2 = PmAz*pGm[2] + pGm[3]*rPmQz;
   double r_200_1 = PmAx*r_100_1 + iz2*(pGm[1] + pGm[2]*rie) + rPmQx*r_100_2;
   double r_020_1 = PmAy*r_010_1 + iz2*(pGm[1] + pGm[2]*rie) + rPmQy*r_010_2;
   double r_002_1 = PmAz*r_001_1 + iz2*(pGm[1] + pGm[2]*rie) + rPmQz*r_001_2;
   double r_110_1 = PmAy*r_100_1 + rPmQy*r_100_2;
   pOut[10] = PmAx*pOut[4] + 2*iz2*(pOut[1] + r_100_1*rie) + rPmQx*r_200_1;
   pOut[11] = PmAy*pOut[5] + 2*iz2*(pOut[2] + r_010_1*rie) + rPmQy*r_020_1;
   pOut[12] = PmAz*pOut[6] + 2*iz2*(pOut[3] + r_001_1*rie) + rPmQz*r_002_1;
   pOut[13] = PmAx*pOut[5] + rPmQx*r_020_1;
   pOut[14] = PmAx*pOut[6] + rPmQx*r_002_1;
   pOut[15] = PmAy*pOut[4] + rPmQy*r_200_1;
   pOut[16] = PmAy*pOut[6] + rPmQy*r_002_1;
   pOut[17] = PmAz*pOut[4] + rPmQz*r_200_1;
   pOut[18] = PmAz*pOut[5] + rPmQz*r_020_1;
   pOut[19] = PmAz*pOut[7] + rPmQz*r_110_1;
   if (lab == 3) return;

   // m + lab == 4
   double r_100_3 = PmAx*pGm[3] + pGm[4]*rPmQx;
   double r_010_3 = PmAy*pGm[3] + pGm[4]*rPmQy;
   double r_001_3 = PmAz*pGm[3] + pGm[4]*rPmQz;
   double r_200_2 = PmAx*r_100_2 + iz2*(pGm[2] + pGm[3]*rie) + rPmQx*r_100_3;
   double r_020_2 = PmAy*r_010_2 + iz2*(pGm[2] + pGm[3]*rie) + rPmQy*r_010_3;
   double r_002_2 = PmAz*r_001_2 + iz2*(pGm[2] + pGm[3]*rie) + rPmQz*r_001_3;
   double r_300_1 = PmAx*r_200_1 + 2*iz2*(r_100_1 + r_100_2*rie) + rPmQx*r_200_2;
   double r_030_1 = PmAy*r_020_1 + 2*iz2*(r_010_1 + r_010_2*rie) + rPmQy*r_020_2;
   double r_003_1 = PmAz*r_002_1 + 2*iz2*(r_001_1 + r_001_2*rie) + rPmQz*r_002_2;
   double r_120_1 = PmAx*r_020_1 + rPmQx*r_020_2;
   double r_102_1 = PmAx*r_002_1 + rPmQx*r_002_2;
   double r_210_1 = PmAy*r_200_1 + rPmQy*r_200_2;
   double r_012_1 = PmAy*r_002_1 + rPmQy*r_002_2;
   pOut[20] = PmAx*pOut[10] + 3*iz2*(pOut[4] + r_200_1*rie) + rPmQx*r_300_1;
   pOut[21] = PmAy*pOut[11] + 3*iz2*(pOut[5] + r_020_1*rie) + rPmQy*r_030_1;
   pOut[22] = PmAz*pOut[12] + 3*iz2*(pOut[6] + r_002_1*rie) + rPmQz*r_003_1;
   pOut[23] = PmAy*pOut[10] + rPmQy*r_300_1;
   pOut[24] = PmAx*pOut[11] + rPmQx*r_030_1;
   pOut[25] = PmAz*pOut[10] + rPmQz*r_300_1;
   pOut[26] = PmAx*pOut[12] + rPmQx*r_003_1;
   pOut[27] = PmAz*pOut[11] + rPmQz*r_030_1;
   pOut[28] = PmAy*pOut[12] + rPmQy*r_003_1;
   pOut[29] = PmAx*pOut[13] + iz2*(pOut[5] + r_020_1*rie) + rPmQx*r_120_1;
   pOut[30] = PmAx*pOut[14] + iz2*(pOut[6] + r_002_1*rie) + rPmQx*r_102_1;
   pOut[31] = PmAy*pOut[16] + iz2*(pOut[6] + r_002_1*rie) + rPmQy*r_012_1;
   pOut[32] = PmAy*pOut[14] + rPmQy*r_102_1;
   pOut[33] = PmAz*pOut[13] + rPmQz*r_120_1;
   pOut[34] = PmAz*pOut[15] + rPmQz*r_210_1;
   if (lab == 4) return;

   // m + lab == 5
   double r_100_4 = PmAx*pGm[4] + pGm[5]*rPmQx;
   double r_010_4 = PmAy*pGm[4] + pGm[5]*rPmQy;
   double r_001_4 = PmAz*pGm[4] + pGm[5]*rPmQz;
   double r_200_3 = PmAx*r_100_3 + iz2*(pGm[3] + pGm[4]*rie) + rPmQx*r_100_4;
   double r_020_3 = PmAy*r_010_3 + iz2*(pGm[3] + pGm[4]*rie) + rPmQy*r_010_4;
   double r_002_3 = PmAz*r_001_3 + iz2*(pGm[3] + pGm[4]*rie) + rPmQz*r_001_4;
   double r_300_2 = PmAx*r_200_2 + 2*iz2*(r_100_2 + r_100_3*rie) + rPmQx*r_200_3;
   double r_030_2 = PmAy*r_020_2 + 2*iz2*(r_010_2 + r_010_3*rie) + rPmQy*r_020_3;
   double r_003_2 = PmAz*r_002_2 + 2*iz2*(r_001_2 + r_001_3*rie) + rPmQz*r_002_3;
   double r_120_2 = PmAx*r_020_2 + rPmQx*r_020_3;
   double r_102_2 = PmAx*r_002_2 + rPmQx*r_002_3;
   double r_012_2 = PmAy*r_002_2 + rPmQy*r_002_3;
   double r_400_1 = PmAx*r_300_1 + 3*iz2*(r_200_1 + r_200_2*rie) + rPmQx*r_300_2;
   double r_040_1 = PmAy*r_030_1 + 3*iz2*(r_020_1 + r_020_2*rie) + rPmQy*r_030_2;
   double r_004_1 = PmAz*r_003_1 + 3*iz2*(r_002_1 + r_002_2*rie) + rPmQz*r_003_2;
   double r_310_1 = PmAy*r_300_1 + rPmQy*r_300_2;
   double r_130_1 = PmAx*r_030_1 + rPmQx*r_030_2;
   double r_301_1 = PmAz*r_300_1 + rPmQz*r_300_2;
   double r_103_1 = PmAx*r_003_1 + rPmQx*r_003_2;
   double r_031_1 = PmAz*r_030_1 + rPmQz*r_030_2;
   double r_013_1 = PmAy*r_003_1 + rPmQy*r_003_2;
   double r_220_1 = PmAx*r_120_1 + iz2*(r_020_1 + r_020_2*rie) + rPmQx*r_120_2;
   double r_202_1 = PmAx*r_102_1 + iz2*(r_002_1 + r_002_2*rie) + rPmQx*r_102_2;
   double r_022_1 = PmAy*r_012_1 + iz2*(r_002_1 + r_002_2*rie) + rPmQy*r_012_2;
   pOut[35] = PmAx*pOut[20] + 4*iz2*(pOut[10] + r_300_1*rie) + rPmQx*r_400_1;
   pOut[36] = PmAy*pOut[21] + 4*iz2*(pOut[11] + r_030_1*rie) + rPmQy*r_040_1;
   pOut[37] = PmAz*pOut[22] + 4*iz2*(pOut[12] + r_003_1*rie) + rPmQz*r_004_1;
   pOut[38] = PmAx*pOut[21] + rPmQx*r_040_1;
   pOut[39] = PmAx*pOut[22] + rPmQx*r_004_1;
   pOut[40] = PmAy*pOut[20] + rPmQy*r_400_1;
   pOut[41] = PmAy*pOut[22] + rPmQy*r_004_1;
   pOut[42] = PmAz*pOut[20] + rPmQz*r_400_1;
   pOut[43] = PmAz*pOut[21] + rPmQz*r_040_1;
   pOut[44] = PmAy*pOut[23] + iz2*(pOut[10] + r_300_1*rie) + rPmQy*r_310_1;
   pOut[45] = PmAz*pOut[25] + iz2*(pOut[10] + r_300_1*rie) + rPmQz*r_301_1;
   pOut[46] = PmAx*pOut[24] + iz2*(pOut[11] + r_030_1*rie) + rPmQx*r_130_1;
   pOut[47] = PmAz*pOut[27] + iz2*(pOut[11] + r_030_1*rie) + rPmQz*r_031_1;
   pOut[48] = PmAx*pOut[26] + iz2*(pOut[12] + r_003_1*rie) + rPmQx*r_103_1;
   pOut[49] = PmAy*pOut[28] + iz2*(pOut[12] + r_003_1*rie) + rPmQy*r_013_1;
   pOut[50] = PmAz*pOut[23] + rPmQz*r_310_1;
   pOut[51] = PmAz*pOut[24] + rPmQz*r_130_1;
   pOut[52] = PmAy*pOut[26] + rPmQy*r_103_1;
   pOut[53] = PmAx*pOut[31] + rPmQx*r_022_1;
   pOut[54] = PmAy*pOut[30] + rPmQy*r_202_1;
   pOut[55] = PmAz*pOut[29] + rPmQz*r_220_1;
   if (lab == 5) return;

   // m + lab == 6
   double r_100_5 = PmAx*pGm[5] + pGm[6]*rPmQx;
   double r_010_5 = PmAy*pGm[5] + pGm[6]*rPmQy;
   double r_001_5 = PmAz*pGm[5] + pGm[6]*rPmQz;
   double r_200_4 = PmAx*r_100_4 + iz2*(pGm[4] + pGm[5]*rie) + rPmQx*r_100_5;
   double r_020_4 = PmAy*r_010_4 + iz2*(pGm[4] + pGm[5]*rie) + rPmQy*r_010_5;
   double r_002_4 = PmAz*r_001_4 + iz2*(pGm[4] + pGm[5]*rie) + rPmQz*r_001_5;
   double r_300_3 = PmAx*r_200_3 + 2*iz2*(r_100_3 + r_100_4*rie) + rPmQx*r_200_4;
   double r_030_3 = PmAy*r_020_3 + 2*iz2*(r_010_3 + r_010_4*rie) + rPmQy*r_020_4;
   double r_003_3 = PmAz*r_002_3 + 2*iz2*(r_001_3 + r_001_4*rie) + rPmQz*r_002_4;
   double r_012_3 = PmAy*r_002_3 + rPmQy*r_002_4;
   double r_400_2 = PmAx*r_300_2 + 3*iz2*(r_200_2 + r_200_3*rie) + rPmQx*r_300_3;
   double r_040_2 = PmAy*r_030_2 + 3*iz2*(r_020_2 + r_020_3*rie) + rPmQy*r_030_3;
   double r_004_2 = PmAz*r_003_2 + 3*iz2*(r_002_2 + r_002_3*rie) + rPmQz*r_003_3;
   double r_310_2 = PmAy*r_300_2 + rPmQy*r_300_3;
   double r_130_2 = PmAx*r_030_2 + rPmQx*r_030_3;
   double r_301_2 = PmAz*r_300_2 + rPmQz*r_300_3;
   double r_103_2 = PmAx*r_003_2 + rPmQx*r_003_3;
   double r_031_2 = PmAz*r_030_2 + rPmQz*r_030_3;
   double r_013_2 = PmAy*r_003_2 + rPmQy*r_003_3;
   double r_022_2 = PmAy*r_012_2 + iz2*(r_002_2 + r_002_3*rie) + rPmQy*r_012_3;
   double r_500_1 = PmAx*r_400_1 + 4*iz2*(r_300_1 + r_300_2*rie) + rPmQx*r_400_2;
   double r_050_1 = PmAy*r_040_1 + 4*iz2*(r_030_1 + r_030_2*rie) + rPmQy*r_040_2;
   double r_005_1 = PmAz*r_004_1 + 4*iz2*(r_003_1 + r_003_2*rie) + rPmQz*r_004_2;
   double r_140_1 = PmAx*r_040_1 + rPmQx*r_040_2;
   double r_104_1 = PmAx*r_004_1 + rPmQx*r_004_2;
   double r_410_1 = PmAy*r_400_1 + rPmQy*r_400_2;
   double r_014_1 = PmAy*r_004_1 + rPmQy*r_004_2;
   double r_401_1 = PmAz*r_400_1 + rPmQz*r_400_2;
   double r_041_1 = PmAz*r_040_1 + rPmQz*r_040_2;
   double r_320_1 = PmAy*r_310_1 + iz2*(r_300_1 + r_300_2*rie) + rPmQy*r_310_2;
   double r_302_1 = PmAz*r_301_1 + iz2*(r_300_1 + r_300_2*rie) + rPmQz*r_301_2;
   double r_230_1 = PmAx*r_130_1 + iz2*(r_030_1 + r_030_2*rie) + rPmQx*r_130_2;
   double r_032_1 = PmAz*r_031_1 + iz2*(r_030_1 + r_030_2*rie) + rPmQz*r_031_2;
   double r_203_1 = PmAx*r_103_1 + iz2*(r_003_1 + r_003_2*rie) + rPmQx*r_103_2;
   double r_023_1 = PmAy*r_013_1 + iz2*(r_003_1 + r_003_2*rie) + rPmQy*r_013_2;
   double r_122_1 = PmAx*r_022_1 + rPmQx*r_022_2;
   pOut[56] = PmAx*pOut[35] + 5*iz2*(pOut[20] + r_400_1*rie) + rPmQx*r_500_1;
   pOut[57] = PmAy*pOut[36] + 5*iz2*(pOut[21] + r_040_1*rie) + rPmQy*r_050_1;
   pOut[58] = PmAz*pOut[37] + 5*iz2*(pOut[22] + r_004_1*rie) + rPmQz*r_005_1;
   pOut[59] = PmAy*pOut[35] + rPmQy*r_500_1;
   pOut[60] = PmAx*pOut[36] + rPmQx*r_050_1;
   pOut[61] = PmAz*pOut[35] + rPmQz*r_500_1;
   pOut[62] = PmAx*pOut[37] + rPmQx*r_005_1;
   pOut[63] = PmAz*pOut[36] + rPmQz*r_050_1;
   pOut[64] = PmAy*pOut[37] + rPmQy*r_005_1;
   pOut[65] = PmAy*pOut[40] + iz2*(pOut[20] + r_400_1*rie) + rPmQy*r_410_1;
   pOut[66] = PmAz*pOut[42] + iz2*(pOut[20] + r_400_1*rie) + rPmQz*r_401_1;
   pOut[67] = PmAx*pOut[38] + iz2*(pOut[21] + r_040_1*rie) + rPmQx*r_140_1;
   pOut[68] = PmAx*pOut[39] + iz2*(pOut[22] + r_004_1*rie) + rPmQx*r_104_1;
   pOut[69] = PmAz*pOut[43] + iz2*(pOut[21] + r_040_1*rie) + rPmQz*r_041_1;
   pOut[70] = PmAy*pOut[41] + iz2*(pOut[22] + r_004_1*rie) + rPmQy*r_014_1;
   pOut[71] = PmAy*pOut[44] + 2*iz2*(pOut[23] + r_310_1*rie) + rPmQy*r_320_1;
   pOut[72] = PmAz*pOut[45] + 2*iz2*(pOut[25] + r_301_1*rie) + rPmQz*r_302_1;
   pOut[73] = PmAz*pOut[47] + 2*iz2*(pOut[27] + r_031_1*rie) + rPmQz*r_032_1;
   pOut[74] = PmAy*pOut[39] + rPmQy*r_104_1;
   pOut[75] = PmAz*pOut[38] + rPmQz*r_140_1;
   pOut[76] = PmAz*pOut[40] + rPmQz*r_410_1;
   pOut[77] = PmAy*pOut[45] + rPmQy*r_302_1;
   pOut[78] = PmAx*pOut[47] + rPmQx*r_032_1;
   pOut[79] = PmAz*pOut[44] + rPmQz*r_320_1;
   pOut[80] = PmAx*pOut[49] + rPmQx*r_023_1;
   pOut[81] = PmAz*pOut[46] + rPmQz*r_230_1;
   pOut[82] = PmAy*pOut[48] + rPmQy*r_203_1;
   pOut[83] = PmAx*pOut[53] + iz2*(pOut[31] + r_022_1*rie) + rPmQx*r_122_1;
   if (lab == 6) return;

   // m + lab == 7
   double r_100_6 = PmAx*pGm[6] + pGm[7]*rPmQx;
   double r_010_6 = PmAy*pGm[6] + pGm[7]*rPmQy;
   double r_001_6 = PmAz*pGm[6] + pGm[7]*rPmQz;
   double r_200_5 = PmAx*r_100_5 + iz2*(pGm[5] + pGm[6]*rie) + rPmQx*r_100_6;
   double r_020_5 = PmAy*r_010_5 + iz2*(pGm[5] + pGm[6]*rie) + rPmQy*r_010_6;
   double r_002_5 = PmAz*r_001_5 + iz2*(pGm[5] + pGm[6]*rie) + rPmQz*r_001_6;
   double r_300_4 = PmAx*r_200_4 + 2*iz2*(r_100_4 + r_100_5*rie) + rPmQx*r_200_5;
   double r_030_4 = PmAy*r_020_4 + 2*iz2*(r_010_4 + r_010_5*rie) + rPmQy*r_020_5;
   double r_003_4 = PmAz*r_002_4 + 2*iz2*(r_001_4 + r_001_5*rie) + rPmQz*r_002_5;
   double r_400_3 = PmAx*r_300_3 + 3*iz2*(r_200_3 + r_200_4*rie) + rPmQx*r_300_4;
   double r_040_3 = PmAy*r_030_3 + 3*iz2*(r_020_3 + r_020_4*rie) + rPmQy*r_030_4;
   double r_004_3 = PmAz*r_003_3 + 3*iz2*(r_002_3 + r_002_4*rie) + rPmQz*r_003_4;
   double r_310_3 = PmAy*r_300_3 + rPmQy*r_300_4;
   double r_301_3 = PmAz*r_300_3 + rPmQz*r_300_4;
   double r_031_3 = PmAz*r_030_3 + rPmQz*r_030_4;
   double r_013_3 = PmAy*r_003_3 + rPmQy*r_003_4;
   double r_500_2 = PmAx*r_400_2 + 4*iz2*(r_300_2 + r_300_3*rie) + rPmQx*r_400_3;
   double r_050_2 = PmAy*r_040_2 + 4*iz2*(r_030_2 + r_030_3*rie) + rPmQy*r_040_3;
   double r_005_2 = PmAz*r_004_2 + 4*iz2*(r_003_2 + r_003_3*rie) + rPmQz*r_004_3;
   double r_140_2 = PmAx*r_040_2 + rPmQx*r_040_3;
   double r_104_2 = PmAx*r_004_2 + rPmQx*r_004_3;
   double r_410_2 = PmAy*r_400_2 + rPmQy*r_400_3;
   double r_014_2 = PmAy*r_004_2 + rPmQy*r_004_3;
   double r_401_2 = PmAz*r_400_2 + rPmQz*r_400_3;
   double r_041_2 = PmAz*r_040_2 + rPmQz*r_040_3;
   double r_320_2 = PmAy*r_310_2 + iz2*(r_300_2 + r_300_3*rie) + rPmQy*r_310_3;
   double r_302_2 = PmAz*r_301_2 + iz2*(r_300_2 + r_300_3*rie) + rPmQz*r_301_3;
   double r_032_2 = PmAz*r_031_2 + iz2*(r_030_2 + r_030_3*rie) + rPmQz*r_031_3;
   double r_023_2 = PmAy*r_013_2 + iz2*(r_003_2 + r_003_3*rie) + rPmQy*r_013_3;
   double r_600_1 = PmAx*r_500_1 + 5*iz2*(r_400_1 + r_400_2*rie) + rPmQx*r_500_2;
   double r_060_1 = PmAy*r_050_1 + 5*iz2*(r_040_1 + r_040_2*rie) + rPmQy*r_050_2;
   double r_006_1 = PmAz*r_005_1 + 5*iz2*(r_004_1 + r_004_2*rie) + rPmQz*r_005_2;
   double r_510_1 = PmAy*r_500_1 + rPmQy*r_500_2;
   double r_150_1 = PmAx*r_050_1 + rPmQx*r_050_2;
   double r_501_1 = PmAz*r_500_1 + rPmQz*r_500_2;
   double r_105_1 = PmAx*r_005_1 + rPmQx*r_005_2;
   double r_051_1 = PmAz*r_050_1 + rPmQz*r_050_2;
   double r_015_1 = PmAy*r_005_1 + rPmQy*r_005_2;
   double r_420_1 = PmAy*r_410_1 + iz2*(r_400_1 + r_400_2*rie) + rPmQy*r_410_2;
   double r_402_1 = PmAz*r_401_1 + iz2*(r_400_1 + r_400_2*rie) + rPmQz*r_401_2;
   double r_240_1 = PmAx*r_140_1 + iz2*(r_040_1 + r_040_2*rie) + rPmQx*r_140_2;
   double r_204_1 = PmAx*r_104_1 + iz2*(r_004_1 + r_004_2*rie) + rPmQx*r_104_2;
   double r_042_1 = PmAz*r_041_1 + iz2*(r_040_1 + r_040_2*rie) + rPmQz*r_041_2;
   double r_024_1 = PmAy*r_014_1 + iz2*(r_004_1 + r_004_2*rie) + rPmQy*r_014_2;
   double r_330_1 = PmAy*r_320_1 + 2*iz2*(r_310_1 + r_310_2*rie) + rPmQy*r_320_2;
   double r_303_1 = PmAz*r_302_1 + 2*iz2*(r_301_1 + r_301_2*rie) + rPmQz*r_302_2;
   double r_033_1 = PmAz*r_032_1 + 2*iz2*(r_031_1 + r_031_2*rie) + rPmQz*r_032_2;
   double r_312_1 = PmAy*r_302_1 + rPmQy*r_302_2;
   double r_132_1 = PmAx*r_032_1 + rPmQx*r_032_2;
   double r_123_1 = PmAx*r_023_1 + rPmQx*r_023_2;
   pOut[84] = PmAx*pOut[56] + 6*iz2*(pOut[35] + r_500_1*rie) + rPmQx*r_600_1;
   pOut[85] = PmAy*pOut[57] + 6*iz2*(pOut[36] + r_050_1*rie) + rPmQy*r_060_1;
   pOut[86] = PmAz*pOut[58] + 6*iz2*(pOut[37] + r_005_1*rie) + rPmQz*r_006_1;
   pOut[87] = PmAx*pOut[57] + rPmQx*r_060_1;
   pOut[88] = PmAx*pOut[58] + rPmQx*r_006_1;
   pOut[89] = PmAy*pOut[56] + rPmQy*r_600_1;
   pOut[90] = PmAy*pOut[58] + rPmQy*r_006_1;
   pOut[91] = PmAz*pOut[56] + rPmQz*r_600_1;
   pOut[92] = PmAz*pOut[57] + rPmQz*r_060_1;
   pOut[93] = PmAy*pOut[59] + iz2*(pOut[35] + r_500_1*rie) + rPmQy*r_510_1;
   pOut[94] = PmAz*pOut[61] + iz2*(pOut[35] + r_500_1*rie) + rPmQz*r_501_1;
   pOut[95] = PmAx*pOut[60] + iz2*(pOut[36] + r_050_1*rie) + rPmQx*r_150_1;
   pOut[96] = PmAz*pOut[63] + iz2*(pOut[36] + r_050_1*rie) + rPmQz*r_051_1;
   pOut[97] = PmAx*pOut[62] + iz2*(pOut[37] + r_005_1*rie) + rPmQx*r_105_1;
   pOut[98] = PmAy*pOut[64] + iz2*(pOut[37] + r_005_1*rie) + rPmQy*r_015_1;
   pOut[99] = PmAx*pOut[67] + 2*iz2*(pOut[38] + r_140_1*rie) + rPmQx*r_240_1;
   pOut[100] = PmAx*pOut[68] + 2*iz2*(pOut[39] + r_104_1*rie) + rPmQx*r_204_1;
   pOut[101] = PmAy*pOut[65] + 2*iz2*(pOut[40] + r_410_1*rie) + rPmQy*r_420_1;
   pOut[102] = PmAy*pOut[70] + 2*iz2*(pOut[41] + r_014_1*rie) + rPmQy*r_024_1;
   pOut[103] = PmAz*pOut[66] + 2*iz2*(pOut[42] + r_401_1*rie) + rPmQz*r_402_1;
   pOut[104] = PmAz*pOut[69] + 2*iz2*(pOut[43] + r_041_1*rie) + rPmQz*r_042_1;
   pOut[105] = PmAz*pOut[59] + rPmQz*r_510_1;
   pOut[106] = PmAz*pOut[60] + rPmQz*r_150_1;
   pOut[107] = PmAy*pOut[62] + rPmQy*r_105_1;
   pOut[108] = PmAx*pOut[69] + rPmQx*r_042_1;
   pOut[109] = PmAx*pOut[70] + rPmQx*r_024_1;
   pOut[110] = PmAy*pOut[66] + rPmQy*r_402_1;
   pOut[111] = PmAy*pOut[68] + rPmQy*r_204_1;
   pOut[112] = PmAz*pOut[65] + rPmQz*r_420_1;
   pOut[113] = PmAz*pOut[67] + rPmQz*r_240_1;
   pOut[114] = PmAz*pOut[71] + rPmQz*r_330_1;
   pOut[115] = PmAy*pOut[72] + rPmQy*r_303_1;
   pOut[116] = PmAx*pOut[73] + rPmQx*r_033_1;
   pOut[117] = PmAy*pOut[77] + iz2*(pOut[45] + r_302_1*rie) + rPmQy*r_312_1;
   pOut[118] = PmAx*pOut[78] + iz2*(pOut[47] + r_032_1*rie) + rPmQx*r_132_1;
   pOut[119] = PmAx*pOut[80] + iz2*(pOut[49] + r_023_1*rie) + rPmQx*r_123_1;
   if (lab == 7) return;

   // m + lab == 8
   double r_100_7 = PmAx*pGm[7] + pGm[8]*rPmQx;
   double r_010_7 = PmAy*pGm[7] + pGm[8]*rPmQy;
   double r_001_7 = PmAz*pGm[7] + pGm[8]*rPmQz;
   double r_200_6 = PmAx*r_100_6 + iz2*(pGm[6] + pGm[7]*rie) + rPmQx*r_100_7;
   double r_020_6 = PmAy*r_010_6 + iz2*(pGm[6] + pGm[7]*rie) + rPmQy*r_010_7;
   double r_002_6 = PmAz*r_001_6 + iz2*(pGm[6] + pGm[7]*rie) + rPmQz*r_001_7;
   double r_300_5 = PmAx*r_200_5 + 2*iz2*(r_100_5 + r_100_6*rie) + rPmQx*r_200_6;
   double r_030_5 = PmAy*r_020_5 + 2*iz2*(r_010_5 + r_010_6*rie) + rPmQy*r_020_6;
   double r_003_5 = PmAz*r_002_5 + 2*iz2*(r_001_5 + r_001_6*rie) + rPmQz*r_002_6;
   double r_400_4 = PmAx*r_300_4 + 3*iz2*(r_200_4 + r_200_5*rie) + rPmQx*r_300_5;
   double r_040_4 = PmAy*r_030_4 + 3*iz2*(r_020_4 + r_020_5*rie) + rPmQy*r_030_5;
   double r_004_4 = PmAz*r_003_4 + 3*iz2*(r_002_4 + r_002_5*rie) + rPmQz*r_003_5;
   double r_310_4 = PmAy*r_300_4 + rPmQy*r_300_5;
   double r_301_4 = PmAz*r_300_4 + rPmQz*r_300_5;
   double r_031_4 = PmAz*r_030_4 + rPmQz*r_030_5;
   double r_500_3 = PmAx*r_400_3 + 4*iz2*(r_300_3 + r_300_4*rie) + rPmQx*r_400_4;
   double r_050_3 = PmAy*r_040_3 + 4*iz2*(r_030_3 + r_030_4*rie) + rPmQy*r_040_4;
   double r_005_3 = PmAz*r_004_3 + 4*iz2*(r_003_3 + r_003_4*rie) + rPmQz*r_004_4;
   double r_140_3 = PmAx*r_040_3 + rPmQx*r_040_4;
   double r_104_3 = PmAx*r_004_3 + rPmQx*r_004_4;
   double r_410_3 = PmAy*r_400_3 + rPmQy*r_400_4;
   double r_014_3 = PmAy*r_004_3 + rPmQy*r_004_4;
   double r_401_3 = PmAz*r_400_3 + rPmQz*r_400_4;
   double r_041_3 = PmAz*r_040_3 + rPmQz*r_040_4;
   double r_320_3 = PmAy*r_310_3 + iz2*(r_300_3 + r_300_4*rie) + rPmQy*r_310_4;
   double r_302_3 = PmAz*r_301_3 + iz2*(r_300_3 + r_300_4*rie) + rPmQz*r_301_4;
   double r_032_3 = PmAz*r_031_3 + iz2*(r_030_3 + r_030_4*rie) + rPmQz*r_031_4;
   double r_600_2 = PmAx*r_500_2 + 5*iz2*(r_400_2 + r_400_3*rie) + rPmQx*r_500_3;
   double r_060_2 = PmAy*r_050_2 + 5*iz2*(r_040_2 + r_040_3*rie) + rPmQy*r_050_3;
   double r_006_2 = PmAz*r_005_2 + 5*iz2*(r_004_2 + r_004_3*rie) + rPmQz*r_005_3;
   double r_510_2 = PmAy*r_500_2 + rPmQy*r_500_3;
   double r_150_2 = PmAx*r_050_2 + rPmQx*r_050_3;
   double r_501_2 = PmAz*r_500_2 + rPmQz*r_500_3;
   double r_105_2 = PmAx*r_005_2 + rPmQx*r_005_3;
   double r_051_2 = PmAz*r_050_2 + rPmQz*r_050_3;
   double r_015_2 = PmAy*r_005_2 + rPmQy*r_005_3;
   double r_420_2 = PmAy*r_410_2 + iz2*(r_400_2 + r_400_3*rie) + rPmQy*r_410_3;
   double r_402_2 = PmAz*r_401_2 + iz2*(r_400_2 + r_400_3*rie) + rPmQz*r_401_3;
   double r_240_2 = PmAx*r_140_2 + iz2*(r_040_2 + r_040_3*rie) + rPmQx*r_140_3;
   double r_204_2 = PmAx*r_104_2 + iz2*(r_004_2 + r_004_3*rie) + rPmQx*r_104_3;
   double r_042_2 = PmAz*r_041_2 + iz2*(r_040_2 + r_040_3*rie) + rPmQz*r_041_3;
   double r_024_2 = PmAy*r_014_2 + iz2*(r_004_2 + r_004_3*rie) + rPmQy*r_014_3;
   double r_330_2 = PmAy*r_320_2 + 2*iz2*(r_310_2 + r_310_3*rie) + rPmQy*r_320_3;
   double r_303_2 = PmAz*r_302_2 + 2*iz2*(r_301_2 + r_301_3*rie) + rPmQz*r_302_3;
   double r_033_2 = PmAz*r_032_2 + 2*iz2*(r_031_2 + r_031_3*rie) + rPmQz*r_032_3;
   double r_700_1 = PmAx*r_600_1 + 6*iz2*(r_500_1 + r_500_2*rie) + rPmQx*r_600_2;
   double r_070_1 = PmAy*r_060_1 + 6*iz2*(r_050_1 + r_050_2*rie) + rPmQy*r_060_2;
   double r_007_1 = PmAz*r_006_1 + 6*iz2*(r_005_1 + r_005_2*rie) + rPmQz*r_006_2;
   double r_160_1 = PmAx*r_060_1 + rPmQx*r_060_2;
   double r_106_1 = PmAx*r_006_1 + rPmQx*r_006_2;
   double r_610_1 = PmAy*r_600_1 + rPmQy*r_600_2;
   double r_016_1 = PmAy*r_006_1 + rPmQy*r_006_2;
   double r_601_1 = PmAz*r_600_1 + rPmQz*r_600_2;
   double r_061_1 = PmAz*r_060_1 + rPmQz*r_060_2;
   double r_520_1 = PmAy*r_510_1 + iz2*(r_500_1 + r_500_2*rie) + rPmQy*r_510_2;
   double r_502_1 = PmAz*r_501_1 + iz2*(r_500_1 + r_500_2*rie) + rPmQz*r_501_2;
   double r_250_1 = PmAx*r_150_1 + iz2*(r_050_1 + r_050_2*rie) + rPmQx*r_150_2;
   double r_052_1 = PmAz*r_051_1 + iz2*(r_050_1 + r_050_2*rie) + rPmQz*r_051_2;
   double r_205_1 = PmAx*r_105_1 + iz2*(r_005_1 + r_005_2*rie) + rPmQx*r_105_2;
   double r_025_1 = PmAy*r_015_1 + iz2*(r_005_1 + r_005_2*rie) + rPmQy*r_015_2;
   double r_340_1 = PmAx*r_240_1 + 2*iz2*(r_140_1 + r_140_2*rie) + rPmQx*r_240_2;
   double r_304_1 = PmAx*r_204_1 + 2*iz2*(r_104_1 + r_104_2*rie) + rPmQx*r_204_2;
   double r_430_1 = PmAy*r_420_1 + 2*iz2*(r_410_1 + r_410_2*rie) + rPmQy*r_420_2;
   double r_034_1 = PmAy*r_024_1 + 2*iz2*(r_014_1 + r_014_2*rie) + rPmQy*r_024_2;
   double r_403_1 = PmAz*r_402_1 + 2*iz2*(r_401_1 + r_401_2*rie) + rPmQz*r_402_2;
   double r_043_1 = PmAz*r_042_1 + 2*iz2*(r_041_1 + r_041_2*rie) + rPmQz*r_042_2;
   double r_142_1 = PmAx*r_042_1 + rPmQx*r_042_2;
   double r_124_1 = PmAx*r_024_1 + rPmQx*r_024_2;
   double r_412_1 = PmAy*r_402_1 + rPmQy*r_402_2;
   double r_331_1 = PmAz*r_330_1 + rPmQz*r_330_2;
   double r_313_1 = PmAy*r_303_1 + rPmQy*r_303_2;
   double r_133_1 = PmAx*r_033_1 + rPmQx*r_033_2;
   pOut[120] = PmAx*pOut[84] + 7*iz2*(pOut[56] + r_600_1*rie) + rPmQx*r_700_1;
   pOut[121] = PmAy*pOut[85] + 7*iz2*(pOut[57] + r_060_1*rie) + rPmQy*r_070_1;
   pOut[122] = PmAz*pOut[86] + 7*iz2*(pOut[58] + r_006_1*rie) + rPmQz*r_007_1;
   pOut[123] = PmAy*pOut[84] + rPmQy*r_700_1;
   pOut[124] = PmAx*pOut[85] + rPmQx*r_070_1;
   pOut[125] = PmAz*pOut[84] + rPmQz*r_700_1;
   pOut[126] = PmAx*pOut[86] + rPmQx*r_007_1;
   pOut[127] = PmAz*pOut[85] + rPmQz*r_070_1;
   pOut[128] = PmAy*pOut[86] + rPmQy*r_007_1;
   pOut[129] = PmAy*pOut[89] + iz2*(pOut[56] + r_600_1*rie) + rPmQy*r_610_1;
   pOut[130] = PmAz*pOut[91] + iz2*(pOut[56] + r_600_1*rie) + rPmQz*r_601_1;
   pOut[131] = PmAx*pOut[87] + iz2*(pOut[57] + r_060_1*rie) + rPmQx*r_160_1;
   pOut[132] = PmAx*pOut[88] + iz2*(pOut[58] + r_006_1*rie) + rPmQx*r_106_1;
   pOut[133] = PmAz*pOut[92] + iz2*(pOut[57] + r_060_1*rie) + rPmQz*r_061_1;
   pOut[134] = PmAy*pOut[90] + iz2*(pOut[58] + r_006_1*rie) + rPmQy*r_016_1;
   pOut[135] = PmAy*pOut[93] + 2*iz2*(pOut[59] + r_510_1*rie) + rPmQy*r_520_1;
   pOut[136] = PmAx*pOut[95] + 2*iz2*(pOut[60] + r_150_1*rie) + rPmQx*r_250_1;
   pOut[137] = PmAz*pOut[94] + 2*iz2*(pOut[61] + r_501_1*rie) + rPmQz*r_502_1;
   pOut[138] = PmAx*pOut[97] + 2*iz2*(pOut[62] + r_105_1*rie) + rPmQx*r_205_1;
   pOut[139] = PmAz*pOut[96] + 2*iz2*(pOut[63] + r_051_1*rie) + rPmQz*r_052_1;
   pOut[140] = PmAy*pOut[98] + 2*iz2*(pOut[64] + r_015_1*rie) + rPmQy*r_025_1;
   pOut[141] = PmAx*pOut[99] + 3*iz2*(pOut[67] + r_240_1*rie) + rPmQx*r_340_1;
   pOut[142] = PmAx*pOut[100] + 3*iz2*(pOut[68] + r_204_1*rie) + rPmQx*r_304_1;
   pOut[143] = PmAy*pOut[102] + 3*iz2*(pOut[70] + r_024_1*rie) + rPmQy*r_034_1;
   pOut[144] = PmAy*pOut[88] + rPmQy*r_106_1;
   pOut[145] = PmAz*pOut[87] + rPmQz*r_160_1;
   pOut[146] = PmAz*pOut[89] + rPmQz*r_610_1;
   pOut[147] = PmAy*pOut[94] + rPmQy*r_502_1;
   pOut[148] = PmAx*pOut[96] + rPmQx*r_052_1;
   pOut[149] = PmAz*pOut[93] + rPmQz*r_520_1;
   pOut[150] = PmAx*pOut[98] + rPmQx*r_025_1;
   pOut[151] = PmAz*pOut[95] + rPmQz*r_250_1;
   pOut[152] = PmAy*pOut[97] + rPmQy*r_205_1;
   pOut[153] = PmAy*pOut[100] + rPmQy*r_304_1;
   pOut[154] = PmAx*pOut[102] + rPmQx*r_034_1;
   pOut[155] = PmAz*pOut[99] + rPmQz*r_340_1;
   pOut[156] = PmAx*pOut[104] + rPmQx*r_043_1;
   pOut[157] = PmAz*pOut[101] + rPmQz*r_430_1;
   pOut[158] = PmAy*pOut[103] + rPmQy*r_403_1;
   pOut[159] = PmAy*pOut[110] + iz2*(pOut[66] + r_402_1*rie) + rPmQy*r_412_1;
   pOut[160] = PmAx*pOut[108] + iz2*(pOut[69] + r_042_1*rie) + rPmQx*r_142_1;
   pOut[161] = PmAx*pOut[109] + iz2*(pOut[70] + r_024_1*rie) + rPmQx*r_124_1;
   pOut[162] = PmAz*pOut[114] + iz2*(pOut[71] + r_330_1*rie) + rPmQz*r_331_1;
   pOut[163] = PmAy*pOut[115] + iz2*(pOut[72] + r_303_1*rie) + rPmQy*r_313_1;
   pOut[164] = PmAx*pOut[116] + iz2*(pOut[73] + r_033_1*rie) + rPmQx*r_133_1;
   if (lab == 8) return;
   // If you get here we ran out of angular momenta. Regenerate with larger Lab.
   assert(0);
}

// Cartesian -> Solid harmonic transforms: Transform a matrix
// N x nCartY(l) to N x (2*l+1).
static void ShTrC5(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      double z0 = pIn[N*0], z1 = pIn[N*1], z3 = pIn[N*3], z5 = pIn[N*5], z7 = pIn[N*7], z8 = pIn[N*8], z9 = pIn[N*9], za = pIn[N*10], zb = pIn[N*11], zc = pIn[N*12], zd = pIn[N*13], ze = pIn[N*14], zf = pIn[N*15], z10 = pIn[N*16], z12 = pIn[N*18], z13 = pIn[N*19], z14 = pIn[N*20];
      pOut[N*0] = pIn[N*4]*sd8 + sd17*z0 + sd17*z3 + sd18*z9 - sd19*z12 - sd19*za;
      pOut[N*1] = pIn[N*6]*sd8 + sd17*z1 + sd17*z5 + sd18*zb - sd19*z13 - sd19*zc;
      pOut[N*2] = -sd1a*z7 + sd1a*z8 + sd1b*zd - sd1b*ze;
      pOut[N*3] = -sd1c*z0 + sd1d*z9 + sd1e*za + sd1f*z3 - sd20*z12;
      pOut[N*4] = -sd21*z10 + sd21*zf;
      pOut[N*5] = sd1c*z1 - sd1d*zb - sd1e*zc - sd1f*z5 + sd20*z13;
      pOut[N*6] = sd22*z7 + sd22*z8 - sd23*z14;
      pOut[N*7] = sd24*z5 - sd25*zb + sd26*z1;
      pOut[N*8] = pIn[N*2] + sd27*z7 + sd27*z8 + sd28*z14 - sd29*zd - sd29*ze;
      pOut[N*9] = sd24*z3 - sd25*z9 + sd26*z0;
      pOut[N*10] = pIn[N*17]*sd2a - sd1b*z10 - sd1b*zf;
      pOut += 1;
      pIn += 1;
   }
   return;
}

static void ShTrC4(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      double z0 = pIn[N*0], z1 = pIn[N*1], z3 = pIn[N*3], z4 = pIn[N*4], z5 = pIn[N*5], z7 = pIn[N*7], z9 = pIn[N*9], za = pIn[N*10], zb = pIn[N*11], zd = pIn[N*13], ze = pIn[N*14];
      pOut[N*0] = pIn[N*2] + sda*z0 + sda*z1 + sdb*z9 - sdc*za - sdc*zb;
      pOut[N*1] = pIn[N*12]*sde - sdd*z3 - sdd*z4;
      pOut[N*2] = pIn[N*6]*sdf - sd7*z5 - sd7*zd;
      pOut[N*3] = sd10*z0 + sd10*z1 - sd11*z9;
      pOut[N*4] = pIn[N*8]*sdf - sd7*z7 - sd7*ze;
      pOut[N*5] = -sd12*z0 + sd12*z1 + sd13*za - sd13*zb;
      pOut[N*6] = sd14*z3 - sd14*z4;
      pOut[N*7] = sd15*z5 - sd16*zd;
      pOut[N*8] = -sd15*z7 + sd16*ze;
      pOut += 1;
      pIn += 1;
   }
   return;
}

static void ShTrC3(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      double z0 = pIn[N*0], z1 = pIn[N*1], z3 = pIn[N*3], z5 = pIn[N*5], z7 = pIn[N*7], z8 = pIn[N*8];
      pOut[N*0] = pIn[N*4]*sd4 - sd3*z0 - sd3*z3;
      pOut[N*1] = pIn[N*6]*sd4 - sd3*z1 - sd3*z5;
      pOut[N*2] = pIn[N*2] - sd5*z7 - sd5*z8;
      pOut[N*3] = sd6*z0 - sd7*z3;
      pOut[N*4] = pIn[N*9]*sd8;
      pOut[N*5] = -sd6*z1 + sd7*z5;
      pOut[N*6] = sd9*z7 - sd9*z8;
      pOut += 1;
      pIn += 1;
   }
   return;
}

static void ShTrC2(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      double z0 = pIn[N*0], z1 = pIn[N*1];
      pOut[N*0] = pIn[N*2] - sd0*z0 - sd0*z1;
      pOut[N*1] = pIn[N*3]*sd1;
      pOut[N*2] = pIn[N*4]*sd1;
      pOut[N*3] = sd2*z0 - sd2*z1;
      pOut[N*4] = pIn[N*5]*sd1;
      pOut += 1;
      pIn += 1;
   }
   return;
}

static void ShTrC1(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      pOut[N*0] = pIn[N*0];
      pOut[N*1] = pIn[N*1];
      pOut[N*2] = pIn[N*2];
      pOut += 1;
      pIn += 1;
   }
   return;
}

static void ShTrC0(double *IR_RP pOut, double const *IR_RP pIn, unsigned N)
{
   for (unsigned i = 0; i < N; ++ i) {
      pOut[N*0] = pIn[N*0];
      pOut += 1;
      pIn += 1;
   }
   return;
}

// Cartesian -> Solid harmonic transform: Transform a matrix
// N x nCartY(l) to N x (2*l+1).
void ShTrN(double *IR_RP pOut, double const *IR_RP pIn, unsigned N, unsigned l)
{
   switch(l) {
      case 0: return ShTrC0(pOut, pIn, N);
      case 1: return ShTrC1(pOut, pIn, N);
      case 2: return ShTrC2(pOut, pIn, N);
      case 3: return ShTrC3(pOut, pIn, N);
      case 4: return ShTrC4(pOut, pIn, N);
      case 5: return ShTrC5(pOut, pIn, N);
   }
   assert(0);
}

// Calculate (a0|c)^m from (a0|c-1x)^{m+1}, for a = la..lab
// i = _x/_y/_z is the reduction direction (x,y,z);
// fPmQi = PmQ[i] * rho/zeta_c; i2e = .5/eta_ABC = .5/(zeta_a + zeta_b + zeta_c)
// note: the la is not for performance reasons (would likely make the actual recurrence
//       slower for small lab) but for preventing overwriting lower components.
static void OsrrB_KerM_x(double *IR_RP pOut, double const *IR_RP pIn, double fPmQx, int la, unsigned lab, double i2e)
{
   if (la < 0)
      la = 0;
   switch(la) {
   case 0:
      pOut[0] = fPmQx*pIn[0];
      if (lab == 0) return;
   case 1:
      pOut[1] = fPmQx*pIn[1] + i2e*pIn[0];
      pOut[2] = fPmQx*pIn[2];
      pOut[3] = fPmQx*pIn[3];
      if (lab == 1) return;
   case 2:
      pOut[4] = fPmQx*pIn[4] + 2*i2e*pIn[1];
      pOut[5] = fPmQx*pIn[5];
      pOut[6] = fPmQx*pIn[6];
      pOut[7] = fPmQx*pIn[7] + i2e*pIn[2];
      pOut[8] = fPmQx*pIn[8] + i2e*pIn[3];
      pOut[9] = fPmQx*pIn[9];
      if (lab == 2) return;
   case 3:
      pOut[10] = fPmQx*pIn[10] + 3*i2e*pIn[4];
      pOut[11] = fPmQx*pIn[11];
      pOut[12] = fPmQx*pIn[12];
      pOut[13] = fPmQx*pIn[13] + i2e*pIn[5];
      pOut[14] = fPmQx*pIn[14] + i2e*pIn[6];
      pOut[15] = fPmQx*pIn[15] + 2*i2e*pIn[7];
      pOut[16] = fPmQx*pIn[16];
      pOut[17] = fPmQx*pIn[17] + 2*i2e*pIn[8];
      pOut[18] = fPmQx*pIn[18];
      pOut[19] = fPmQx*pIn[19] + i2e*pIn[9];
      if (lab == 3) return;
   case 4:
      pOut[20] = fPmQx*pIn[20] + 4*i2e*pIn[10];
      pOut[21] = fPmQx*pIn[21];
      pOut[22] = fPmQx*pIn[22];
      pOut[23] = fPmQx*pIn[23] + 3*i2e*pIn[15];
      pOut[24] = fPmQx*pIn[24] + i2e*pIn[11];
      pOut[25] = fPmQx*pIn[25] + 3*i2e*pIn[17];
      pOut[26] = fPmQx*pIn[26] + i2e*pIn[12];
      pOut[27] = fPmQx*pIn[27];
      pOut[28] = fPmQx*pIn[28];
      pOut[29] = fPmQx*pIn[29] + 2*i2e*pIn[13];
      pOut[30] = fPmQx*pIn[30] + 2*i2e*pIn[14];
      pOut[31] = fPmQx*pIn[31];
      pOut[32] = fPmQx*pIn[32] + i2e*pIn[16];
      pOut[33] = fPmQx*pIn[33] + i2e*pIn[18];
      pOut[34] = fPmQx*pIn[34] + 2*i2e*pIn[19];
      if (lab == 4) return;
   case 5:
      pOut[35] = fPmQx*pIn[35] + 5*i2e*pIn[20];
      pOut[36] = fPmQx*pIn[36];
      pOut[37] = fPmQx*pIn[37];
      pOut[38] = fPmQx*pIn[38] + i2e*pIn[21];
      pOut[39] = fPmQx*pIn[39] + i2e*pIn[22];
      pOut[40] = fPmQx*pIn[40] + 4*i2e*pIn[23];
      pOut[41] = fPmQx*pIn[41];
      pOut[42] = fPmQx*pIn[42] + 4*i2e*pIn[25];
      pOut[43] = fPmQx*pIn[43];
      pOut[44] = fPmQx*pIn[44] + 3*i2e*pIn[29];
      pOut[45] = fPmQx*pIn[45] + 3*i2e*pIn[30];
      pOut[46] = fPmQx*pIn[46] + 2*i2e*pIn[24];
      pOut[47] = fPmQx*pIn[47];
      pOut[48] = fPmQx*pIn[48] + 2*i2e*pIn[26];
      pOut[49] = fPmQx*pIn[49];
      pOut[50] = fPmQx*pIn[50] + 3*i2e*pIn[34];
      pOut[51] = fPmQx*pIn[51] + i2e*pIn[27];
      pOut[52] = fPmQx*pIn[52] + i2e*pIn[28];
      pOut[53] = fPmQx*pIn[53] + i2e*pIn[31];
      pOut[54] = fPmQx*pIn[54] + 2*i2e*pIn[32];
      pOut[55] = fPmQx*pIn[55] + 2*i2e*pIn[33];
      if (lab == 5) return;
   case 6:
      pOut[56] = fPmQx*pIn[56] + 6*i2e*pIn[35];
      pOut[57] = fPmQx*pIn[57];
      pOut[58] = fPmQx*pIn[58];
      pOut[59] = fPmQx*pIn[59] + 5*i2e*pIn[40];
      pOut[60] = fPmQx*pIn[60] + i2e*pIn[36];
      pOut[61] = fPmQx*pIn[61] + 5*i2e*pIn[42];
      pOut[62] = fPmQx*pIn[62] + i2e*pIn[37];
      pOut[63] = fPmQx*pIn[63];
      pOut[64] = fPmQx*pIn[64];
      pOut[65] = fPmQx*pIn[65] + 4*i2e*pIn[44];
      pOut[66] = fPmQx*pIn[66] + 4*i2e*pIn[45];
      pOut[67] = fPmQx*pIn[67] + 2*i2e*pIn[38];
      pOut[68] = fPmQx*pIn[68] + 2*i2e*pIn[39];
      pOut[69] = fPmQx*pIn[69];
      pOut[70] = fPmQx*pIn[70];
      pOut[71] = fPmQx*pIn[71] + 3*i2e*pIn[46];
      pOut[72] = fPmQx*pIn[72] + 3*i2e*pIn[48];
      pOut[73] = fPmQx*pIn[73];
      pOut[74] = fPmQx*pIn[74] + i2e*pIn[41];
      pOut[75] = fPmQx*pIn[75] + i2e*pIn[43];
      pOut[76] = fPmQx*pIn[76] + 4*i2e*pIn[50];
      pOut[77] = fPmQx*pIn[77] + 3*i2e*pIn[54];
      pOut[78] = fPmQx*pIn[78] + i2e*pIn[47];
      pOut[79] = fPmQx*pIn[79] + 3*i2e*pIn[55];
      pOut[80] = fPmQx*pIn[80] + i2e*pIn[49];
      pOut[81] = fPmQx*pIn[81] + 2*i2e*pIn[51];
      pOut[82] = fPmQx*pIn[82] + 2*i2e*pIn[52];
      pOut[83] = fPmQx*pIn[83] + 2*i2e*pIn[53];
      if (lab == 6) return;
   case 7:
      pOut[84] = fPmQx*pIn[84] + 7*i2e*pIn[56];
      pOut[85] = fPmQx*pIn[85];
      pOut[86] = fPmQx*pIn[86];
      pOut[87] = fPmQx*pIn[87] + i2e*pIn[57];
      pOut[88] = fPmQx*pIn[88] + i2e*pIn[58];
      pOut[89] = fPmQx*pIn[89] + 6*i2e*pIn[59];
      pOut[90] = fPmQx*pIn[90];
      pOut[91] = fPmQx*pIn[91] + 6*i2e*pIn[61];
      pOut[92] = fPmQx*pIn[92];
      pOut[93] = fPmQx*pIn[93] + 5*i2e*pIn[65];
      pOut[94] = fPmQx*pIn[94] + 5*i2e*pIn[66];
      pOut[95] = fPmQx*pIn[95] + 2*i2e*pIn[60];
      pOut[96] = fPmQx*pIn[96];
      pOut[97] = fPmQx*pIn[97] + 2*i2e*pIn[62];
      pOut[98] = fPmQx*pIn[98];
      pOut[99] = fPmQx*pIn[99] + 3*i2e*pIn[67];
      pOut[100] = fPmQx*pIn[100] + 3*i2e*pIn[68];
      pOut[101] = fPmQx*pIn[101] + 4*i2e*pIn[71];
      pOut[102] = fPmQx*pIn[102];
      pOut[103] = fPmQx*pIn[103] + 4*i2e*pIn[72];
      pOut[104] = fPmQx*pIn[104];
      pOut[105] = fPmQx*pIn[105] + 5*i2e*pIn[76];
      pOut[106] = fPmQx*pIn[106] + i2e*pIn[63];
      pOut[107] = fPmQx*pIn[107] + i2e*pIn[64];
      pOut[108] = fPmQx*pIn[108] + i2e*pIn[69];
      pOut[109] = fPmQx*pIn[109] + i2e*pIn[70];
      pOut[110] = fPmQx*pIn[110] + 4*i2e*pIn[77];
      pOut[111] = fPmQx*pIn[111] + 2*i2e*pIn[74];
      pOut[112] = fPmQx*pIn[112] + 4*i2e*pIn[79];
      pOut[113] = fPmQx*pIn[113] + 2*i2e*pIn[75];
      pOut[114] = fPmQx*pIn[114] + 3*i2e*pIn[81];
      pOut[115] = fPmQx*pIn[115] + 3*i2e*pIn[82];
      pOut[116] = fPmQx*pIn[116] + i2e*pIn[73];
      pOut[117] = fPmQx*pIn[117] + 3*i2e*pIn[83];
      pOut[118] = fPmQx*pIn[118] + 2*i2e*pIn[78];
      pOut[119] = fPmQx*pIn[119] + 2*i2e*pIn[80];
      if (lab == 7) return;
   case 8:
      pOut[120] = fPmQx*pIn[120] + 8*i2e*pIn[84];
      pOut[121] = fPmQx*pIn[121];
      pOut[122] = fPmQx*pIn[122];
      pOut[123] = fPmQx*pIn[123] + 7*i2e*pIn[89];
      pOut[124] = fPmQx*pIn[124] + i2e*pIn[85];
      pOut[125] = fPmQx*pIn[125] + 7*i2e*pIn[91];
      pOut[126] = fPmQx*pIn[126] + i2e*pIn[86];
      pOut[127] = fPmQx*pIn[127];
      pOut[128] = fPmQx*pIn[128];
      pOut[129] = fPmQx*pIn[129] + 6*i2e*pIn[93];
      pOut[130] = fPmQx*pIn[130] + 6*i2e*pIn[94];
      pOut[131] = fPmQx*pIn[131] + 2*i2e*pIn[87];
      pOut[132] = fPmQx*pIn[132] + 2*i2e*pIn[88];
      pOut[133] = fPmQx*pIn[133];
      pOut[134] = fPmQx*pIn[134];
      pOut[135] = fPmQx*pIn[135] + 5*i2e*pIn[101];
      pOut[136] = fPmQx*pIn[136] + 3*i2e*pIn[95];
      pOut[137] = fPmQx*pIn[137] + 5*i2e*pIn[103];
      pOut[138] = fPmQx*pIn[138] + 3*i2e*pIn[97];
      pOut[139] = fPmQx*pIn[139];
      pOut[140] = fPmQx*pIn[140];
      pOut[141] = fPmQx*pIn[141] + 4*i2e*pIn[99];
      pOut[142] = fPmQx*pIn[142] + 4*i2e*pIn[100];
      pOut[143] = fPmQx*pIn[143];
      pOut[144] = fPmQx*pIn[144] + i2e*pIn[90];
      pOut[145] = fPmQx*pIn[145] + i2e*pIn[92];
      pOut[146] = fPmQx*pIn[146] + 6*i2e*pIn[105];
      pOut[147] = fPmQx*pIn[147] + 5*i2e*pIn[110];
      pOut[148] = fPmQx*pIn[148] + i2e*pIn[96];
      pOut[149] = fPmQx*pIn[149] + 5*i2e*pIn[112];
      pOut[150] = fPmQx*pIn[150] + i2e*pIn[98];
      pOut[151] = fPmQx*pIn[151] + 2*i2e*pIn[106];
      pOut[152] = fPmQx*pIn[152] + 2*i2e*pIn[107];
      pOut[153] = fPmQx*pIn[153] + 3*i2e*pIn[111];
      pOut[154] = fPmQx*pIn[154] + i2e*pIn[102];
      pOut[155] = fPmQx*pIn[155] + 3*i2e*pIn[113];
      pOut[156] = fPmQx*pIn[156] + i2e*pIn[104];
      pOut[157] = fPmQx*pIn[157] + 4*i2e*pIn[114];
      pOut[158] = fPmQx*pIn[158] + 4*i2e*pIn[115];
      pOut[159] = fPmQx*pIn[159] + 4*i2e*pIn[117];
      pOut[160] = fPmQx*pIn[160] + 2*i2e*pIn[108];
      pOut[161] = fPmQx*pIn[161] + 2*i2e*pIn[109];
      pOut[162] = fPmQx*pIn[162] + 3*i2e*pIn[118];
      pOut[163] = fPmQx*pIn[163] + 3*i2e*pIn[119];
      pOut[164] = fPmQx*pIn[164] + 2*i2e*pIn[116];
      if (lab == 8) return;
   }
   assert(0);
}

// Calculate (a0|c)^m from (a0|c-1y)^{m+1}, for a = la..lab
// i = _x/_y/_z is the reduction direction (x,y,z);
// fPmQi = PmQ[i] * rho/zeta_c; i2e = .5/eta_ABC = .5/(zeta_a + zeta_b + zeta_c)
// note: the la is not for performance reasons (would likely make the actual recurrence
//       slower for small lab) but for preventing overwriting lower components.
static void OsrrB_KerM_y(double *IR_RP pOut, double const *IR_RP pIn, double fPmQy, int la, unsigned lab, double i2e)
{
   if (la < 0)
      la = 0;
   switch(la) {
   case 0:
      pOut[0] = fPmQy*pIn[0];
      if (lab == 0) return;
   case 1:
      pOut[1] = fPmQy*pIn[1];
      pOut[2] = fPmQy*pIn[2] + i2e*pIn[0];
      pOut[3] = fPmQy*pIn[3];
      if (lab == 1) return;
   case 2:
      pOut[4] = fPmQy*pIn[4];
      pOut[5] = fPmQy*pIn[5] + 2*i2e*pIn[2];
      pOut[6] = fPmQy*pIn[6];
      pOut[7] = fPmQy*pIn[7] + i2e*pIn[1];
      pOut[8] = fPmQy*pIn[8];
      pOut[9] = fPmQy*pIn[9] + i2e*pIn[3];
      if (lab == 2) return;
   case 3:
      pOut[10] = fPmQy*pIn[10];
      pOut[11] = fPmQy*pIn[11] + 3*i2e*pIn[5];
      pOut[12] = fPmQy*pIn[12];
      pOut[13] = fPmQy*pIn[13] + 2*i2e*pIn[7];
      pOut[14] = fPmQy*pIn[14];
      pOut[15] = fPmQy*pIn[15] + i2e*pIn[4];
      pOut[16] = fPmQy*pIn[16] + i2e*pIn[6];
      pOut[17] = fPmQy*pIn[17];
      pOut[18] = fPmQy*pIn[18] + 2*i2e*pIn[9];
      pOut[19] = fPmQy*pIn[19] + i2e*pIn[8];
      if (lab == 3) return;
   case 4:
      pOut[20] = fPmQy*pIn[20];
      pOut[21] = fPmQy*pIn[21] + 4*i2e*pIn[11];
      pOut[22] = fPmQy*pIn[22];
      pOut[23] = fPmQy*pIn[23] + i2e*pIn[10];
      pOut[24] = fPmQy*pIn[24] + 3*i2e*pIn[13];
      pOut[25] = fPmQy*pIn[25];
      pOut[26] = fPmQy*pIn[26];
      pOut[27] = fPmQy*pIn[27] + 3*i2e*pIn[18];
      pOut[28] = fPmQy*pIn[28] + i2e*pIn[12];
      pOut[29] = fPmQy*pIn[29] + 2*i2e*pIn[15];
      pOut[30] = fPmQy*pIn[30];
      pOut[31] = fPmQy*pIn[31] + 2*i2e*pIn[16];
      pOut[32] = fPmQy*pIn[32] + i2e*pIn[14];
      pOut[33] = fPmQy*pIn[33] + 2*i2e*pIn[19];
      pOut[34] = fPmQy*pIn[34] + i2e*pIn[17];
      if (lab == 4) return;
   case 5:
      pOut[35] = fPmQy*pIn[35];
      pOut[36] = fPmQy*pIn[36] + 5*i2e*pIn[21];
      pOut[37] = fPmQy*pIn[37];
      pOut[38] = fPmQy*pIn[38] + 4*i2e*pIn[24];
      pOut[39] = fPmQy*pIn[39];
      pOut[40] = fPmQy*pIn[40] + i2e*pIn[20];
      pOut[41] = fPmQy*pIn[41] + i2e*pIn[22];
      pOut[42] = fPmQy*pIn[42];
      pOut[43] = fPmQy*pIn[43] + 4*i2e*pIn[27];
      pOut[44] = fPmQy*pIn[44] + 2*i2e*pIn[23];
      pOut[45] = fPmQy*pIn[45];
      pOut[46] = fPmQy*pIn[46] + 3*i2e*pIn[29];
      pOut[47] = fPmQy*pIn[47] + 3*i2e*pIn[31];
      pOut[48] = fPmQy*pIn[48];
      pOut[49] = fPmQy*pIn[49] + 2*i2e*pIn[28];
      pOut[50] = fPmQy*pIn[50] + i2e*pIn[25];
      pOut[51] = fPmQy*pIn[51] + 3*i2e*pIn[33];
      pOut[52] = fPmQy*pIn[52] + i2e*pIn[26];
      pOut[53] = fPmQy*pIn[53] + 2*i2e*pIn[32];
      pOut[54] = fPmQy*pIn[54] + i2e*pIn[30];
      pOut[55] = fPmQy*pIn[55] + 2*i2e*pIn[34];
      if (lab == 5) return;
   case 6:
      pOut[56] = fPmQy*pIn[56];
      pOut[57] = fPmQy*pIn[57] + 6*i2e*pIn[36];
      pOut[58] = fPmQy*pIn[58];
      pOut[59] = fPmQy*pIn[59] + i2e*pIn[35];
      pOut[60] = fPmQy*pIn[60] + 5*i2e*pIn[38];
      pOut[61] = fPmQy*pIn[61];
      pOut[62] = fPmQy*pIn[62];
      pOut[63] = fPmQy*pIn[63] + 5*i2e*pIn[43];
      pOut[64] = fPmQy*pIn[64] + i2e*pIn[37];
      pOut[65] = fPmQy*pIn[65] + 2*i2e*pIn[40];
      pOut[66] = fPmQy*pIn[66];
      pOut[67] = fPmQy*pIn[67] + 4*i2e*pIn[46];
      pOut[68] = fPmQy*pIn[68];
      pOut[69] = fPmQy*pIn[69] + 4*i2e*pIn[47];
      pOut[70] = fPmQy*pIn[70] + 2*i2e*pIn[41];
      pOut[71] = fPmQy*pIn[71] + 3*i2e*pIn[44];
      pOut[72] = fPmQy*pIn[72];
      pOut[73] = fPmQy*pIn[73] + 3*i2e*pIn[49];
      pOut[74] = fPmQy*pIn[74] + i2e*pIn[39];
      pOut[75] = fPmQy*pIn[75] + 4*i2e*pIn[51];
      pOut[76] = fPmQy*pIn[76] + i2e*pIn[42];
      pOut[77] = fPmQy*pIn[77] + i2e*pIn[45];
      pOut[78] = fPmQy*pIn[78] + 3*i2e*pIn[53];
      pOut[79] = fPmQy*pIn[79] + 2*i2e*pIn[50];
      pOut[80] = fPmQy*pIn[80] + 2*i2e*pIn[52];
      pOut[81] = fPmQy*pIn[81] + 3*i2e*pIn[55];
      pOut[82] = fPmQy*pIn[82] + i2e*pIn[48];
      pOut[83] = fPmQy*pIn[83] + 2*i2e*pIn[54];
      if (lab == 6) return;
   case 7:
      pOut[84] = fPmQy*pIn[84];
      pOut[85] = fPmQy*pIn[85] + 7*i2e*pIn[57];
      pOut[86] = fPmQy*pIn[86];
      pOut[87] = fPmQy*pIn[87] + 6*i2e*pIn[60];
      pOut[88] = fPmQy*pIn[88];
      pOut[89] = fPmQy*pIn[89] + i2e*pIn[56];
      pOut[90] = fPmQy*pIn[90] + i2e*pIn[58];
      pOut[91] = fPmQy*pIn[91];
      pOut[92] = fPmQy*pIn[92] + 6*i2e*pIn[63];
      pOut[93] = fPmQy*pIn[93] + 2*i2e*pIn[59];
      pOut[94] = fPmQy*pIn[94];
      pOut[95] = fPmQy*pIn[95] + 5*i2e*pIn[67];
      pOut[96] = fPmQy*pIn[96] + 5*i2e*pIn[69];
      pOut[97] = fPmQy*pIn[97];
      pOut[98] = fPmQy*pIn[98] + 2*i2e*pIn[64];
      pOut[99] = fPmQy*pIn[99] + 4*i2e*pIn[71];
      pOut[100] = fPmQy*pIn[100];
      pOut[101] = fPmQy*pIn[101] + 3*i2e*pIn[65];
      pOut[102] = fPmQy*pIn[102] + 3*i2e*pIn[70];
      pOut[103] = fPmQy*pIn[103];
      pOut[104] = fPmQy*pIn[104] + 4*i2e*pIn[73];
      pOut[105] = fPmQy*pIn[105] + i2e*pIn[61];
      pOut[106] = fPmQy*pIn[106] + 5*i2e*pIn[75];
      pOut[107] = fPmQy*pIn[107] + i2e*pIn[62];
      pOut[108] = fPmQy*pIn[108] + 4*i2e*pIn[78];
      pOut[109] = fPmQy*pIn[109] + 2*i2e*pIn[74];
      pOut[110] = fPmQy*pIn[110] + i2e*pIn[66];
      pOut[111] = fPmQy*pIn[111] + i2e*pIn[68];
      pOut[112] = fPmQy*pIn[112] + 2*i2e*pIn[76];
      pOut[113] = fPmQy*pIn[113] + 4*i2e*pIn[81];
      pOut[114] = fPmQy*pIn[114] + 3*i2e*pIn[79];
      pOut[115] = fPmQy*pIn[115] + i2e*pIn[72];
      pOut[116] = fPmQy*pIn[116] + 3*i2e*pIn[80];
      pOut[117] = fPmQy*pIn[117] + 2*i2e*pIn[77];
      pOut[118] = fPmQy*pIn[118] + 3*i2e*pIn[83];
      pOut[119] = fPmQy*pIn[119] + 2*i2e*pIn[82];
      if (lab == 7) return;
   case 8:
      pOut[120] = fPmQy*pIn[120];
      pOut[121] = fPmQy*pIn[121] + 8*i2e*pIn[85];
      pOut[122] = fPmQy*pIn[122];
      pOut[123] = fPmQy*pIn[123] + i2e*pIn[84];
      pOut[124] = fPmQy*pIn[124] + 7*i2e*pIn[87];
      pOut[125] = fPmQy*pIn[125];
      pOut[126] = fPmQy*pIn[126];
      pOut[127] = fPmQy*pIn[127] + 7*i2e*pIn[92];
      pOut[128] = fPmQy*pIn[128] + i2e*pIn[86];
      pOut[129] = fPmQy*pIn[129] + 2*i2e*pIn[89];
      pOut[130] = fPmQy*pIn[130];
      pOut[131] = fPmQy*pIn[131] + 6*i2e*pIn[95];
      pOut[132] = fPmQy*pIn[132];
      pOut[133] = fPmQy*pIn[133] + 6*i2e*pIn[96];
      pOut[134] = fPmQy*pIn[134] + 2*i2e*pIn[90];
      pOut[135] = fPmQy*pIn[135] + 3*i2e*pIn[93];
      pOut[136] = fPmQy*pIn[136] + 5*i2e*pIn[99];
      pOut[137] = fPmQy*pIn[137];
      pOut[138] = fPmQy*pIn[138];
      pOut[139] = fPmQy*pIn[139] + 5*i2e*pIn[104];
      pOut[140] = fPmQy*pIn[140] + 3*i2e*pIn[98];
      pOut[141] = fPmQy*pIn[141] + 4*i2e*pIn[101];
      pOut[142] = fPmQy*pIn[142];
      pOut[143] = fPmQy*pIn[143] + 4*i2e*pIn[102];
      pOut[144] = fPmQy*pIn[144] + i2e*pIn[88];
      pOut[145] = fPmQy*pIn[145] + 6*i2e*pIn[106];
      pOut[146] = fPmQy*pIn[146] + i2e*pIn[91];
      pOut[147] = fPmQy*pIn[147] + i2e*pIn[94];
      pOut[148] = fPmQy*pIn[148] + 5*i2e*pIn[108];
      pOut[149] = fPmQy*pIn[149] + 2*i2e*pIn[105];
      pOut[150] = fPmQy*pIn[150] + 2*i2e*pIn[107];
      pOut[151] = fPmQy*pIn[151] + 5*i2e*pIn[113];
      pOut[152] = fPmQy*pIn[152] + i2e*pIn[97];
      pOut[153] = fPmQy*pIn[153] + i2e*pIn[100];
      pOut[154] = fPmQy*pIn[154] + 3*i2e*pIn[109];
      pOut[155] = fPmQy*pIn[155] + 4*i2e*pIn[114];
      pOut[156] = fPmQy*pIn[156] + 4*i2e*pIn[116];
      pOut[157] = fPmQy*pIn[157] + 3*i2e*pIn[112];
      pOut[158] = fPmQy*pIn[158] + i2e*pIn[103];
      pOut[159] = fPmQy*pIn[159] + 2*i2e*pIn[110];
      pOut[160] = fPmQy*pIn[160] + 4*i2e*pIn[118];
      pOut[161] = fPmQy*pIn[161] + 2*i2e*pIn[111];
      pOut[162] = fPmQy*pIn[162] + 3*i2e*pIn[117];
      pOut[163] = fPmQy*pIn[163] + 2*i2e*pIn[115];
      pOut[164] = fPmQy*pIn[164] + 3*i2e*pIn[119];
      if (lab == 8) return;
   }
   assert(0);
}

// Calculate (a0|c)^m from (a0|c-1z)^{m+1}, for a = la..lab
// i = _x/_y/_z is the reduction direction (x,y,z);
// fPmQi = PmQ[i] * rho/zeta_c; i2e = .5/eta_ABC = .5/(zeta_a + zeta_b + zeta_c)
// note: the la is not for performance reasons (would likely make the actual recurrence
//       slower for small lab) but for preventing overwriting lower components.
static void OsrrB_KerM_z(double *IR_RP pOut, double const *IR_RP pIn, double fPmQz, int la, unsigned lab, double i2e)
{
   if (la < 0)
      la = 0;
   switch(la) {
   case 0:
      pOut[0] = fPmQz*pIn[0];
      if (lab == 0) return;
   case 1:
      pOut[1] = fPmQz*pIn[1];
      pOut[2] = fPmQz*pIn[2];
      pOut[3] = fPmQz*pIn[3] + i2e*pIn[0];
      if (lab == 1) return;
   case 2:
      pOut[4] = fPmQz*pIn[4];
      pOut[5] = fPmQz*pIn[5];
      pOut[6] = fPmQz*pIn[6] + 2*i2e*pIn[3];
      pOut[7] = fPmQz*pIn[7];
      pOut[8] = fPmQz*pIn[8] + i2e*pIn[1];
      pOut[9] = fPmQz*pIn[9] + i2e*pIn[2];
      if (lab == 2) return;
   case 3:
      pOut[10] = fPmQz*pIn[10];
      pOut[11] = fPmQz*pIn[11];
      pOut[12] = fPmQz*pIn[12] + 3*i2e*pIn[6];
      pOut[13] = fPmQz*pIn[13];
      pOut[14] = fPmQz*pIn[14] + 2*i2e*pIn[8];
      pOut[15] = fPmQz*pIn[15];
      pOut[16] = fPmQz*pIn[16] + 2*i2e*pIn[9];
      pOut[17] = fPmQz*pIn[17] + i2e*pIn[4];
      pOut[18] = fPmQz*pIn[18] + i2e*pIn[5];
      pOut[19] = fPmQz*pIn[19] + i2e*pIn[7];
      if (lab == 3) return;
   case 4:
      pOut[20] = fPmQz*pIn[20];
      pOut[21] = fPmQz*pIn[21];
      pOut[22] = fPmQz*pIn[22] + 4*i2e*pIn[12];
      pOut[23] = fPmQz*pIn[23];
      pOut[24] = fPmQz*pIn[24];
      pOut[25] = fPmQz*pIn[25] + i2e*pIn[10];
      pOut[26] = fPmQz*pIn[26] + 3*i2e*pIn[14];
      pOut[27] = fPmQz*pIn[27] + i2e*pIn[11];
      pOut[28] = fPmQz*pIn[28] + 3*i2e*pIn[16];
      pOut[29] = fPmQz*pIn[29];
      pOut[30] = fPmQz*pIn[30] + 2*i2e*pIn[17];
      pOut[31] = fPmQz*pIn[31] + 2*i2e*pIn[18];
      pOut[32] = fPmQz*pIn[32] + 2*i2e*pIn[19];
      pOut[33] = fPmQz*pIn[33] + i2e*pIn[13];
      pOut[34] = fPmQz*pIn[34] + i2e*pIn[15];
      if (lab == 4) return;
   case 5:
      pOut[35] = fPmQz*pIn[35];
      pOut[36] = fPmQz*pIn[36];
      pOut[37] = fPmQz*pIn[37] + 5*i2e*pIn[22];
      pOut[38] = fPmQz*pIn[38];
      pOut[39] = fPmQz*pIn[39] + 4*i2e*pIn[26];
      pOut[40] = fPmQz*pIn[40];
      pOut[41] = fPmQz*pIn[41] + 4*i2e*pIn[28];
      pOut[42] = fPmQz*pIn[42] + i2e*pIn[20];
      pOut[43] = fPmQz*pIn[43] + i2e*pIn[21];
      pOut[44] = fPmQz*pIn[44];
      pOut[45] = fPmQz*pIn[45] + 2*i2e*pIn[25];
      pOut[46] = fPmQz*pIn[46];
      pOut[47] = fPmQz*pIn[47] + 2*i2e*pIn[27];
      pOut[48] = fPmQz*pIn[48] + 3*i2e*pIn[30];
      pOut[49] = fPmQz*pIn[49] + 3*i2e*pIn[31];
      pOut[50] = fPmQz*pIn[50] + i2e*pIn[23];
      pOut[51] = fPmQz*pIn[51] + i2e*pIn[24];
      pOut[52] = fPmQz*pIn[52] + 3*i2e*pIn[32];
      pOut[53] = fPmQz*pIn[53] + 2*i2e*pIn[33];
      pOut[54] = fPmQz*pIn[54] + 2*i2e*pIn[34];
      pOut[55] = fPmQz*pIn[55] + i2e*pIn[29];
      if (lab == 5) return;
   case 6:
      pOut[56] = fPmQz*pIn[56];
      pOut[57] = fPmQz*pIn[57];
      pOut[58] = fPmQz*pIn[58] + 6*i2e*pIn[37];
      pOut[59] = fPmQz*pIn[59];
      pOut[60] = fPmQz*pIn[60];
      pOut[61] = fPmQz*pIn[61] + i2e*pIn[35];
      pOut[62] = fPmQz*pIn[62] + 5*i2e*pIn[39];
      pOut[63] = fPmQz*pIn[63] + i2e*pIn[36];
      pOut[64] = fPmQz*pIn[64] + 5*i2e*pIn[41];
      pOut[65] = fPmQz*pIn[65];
      pOut[66] = fPmQz*pIn[66] + 2*i2e*pIn[42];
      pOut[67] = fPmQz*pIn[67];
      pOut[68] = fPmQz*pIn[68] + 4*i2e*pIn[48];
      pOut[69] = fPmQz*pIn[69] + 2*i2e*pIn[43];
      pOut[70] = fPmQz*pIn[70] + 4*i2e*pIn[49];
      pOut[71] = fPmQz*pIn[71];
      pOut[72] = fPmQz*pIn[72] + 3*i2e*pIn[45];
      pOut[73] = fPmQz*pIn[73] + 3*i2e*pIn[47];
      pOut[74] = fPmQz*pIn[74] + 4*i2e*pIn[52];
      pOut[75] = fPmQz*pIn[75] + i2e*pIn[38];
      pOut[76] = fPmQz*pIn[76] + i2e*pIn[40];
      pOut[77] = fPmQz*pIn[77] + 2*i2e*pIn[50];
      pOut[78] = fPmQz*pIn[78] + 2*i2e*pIn[51];
      pOut[79] = fPmQz*pIn[79] + i2e*pIn[44];
      pOut[80] = fPmQz*pIn[80] + 3*i2e*pIn[53];
      pOut[81] = fPmQz*pIn[81] + i2e*pIn[46];
      pOut[82] = fPmQz*pIn[82] + 3*i2e*pIn[54];
      pOut[83] = fPmQz*pIn[83] + 2*i2e*pIn[55];
      if (lab == 6) return;
   case 7:
      pOut[84] = fPmQz*pIn[84];
      pOut[85] = fPmQz*pIn[85];
      pOut[86] = fPmQz*pIn[86] + 7*i2e*pIn[58];
      pOut[87] = fPmQz*pIn[87];
      pOut[88] = fPmQz*pIn[88] + 6*i2e*pIn[62];
      pOut[89] = fPmQz*pIn[89];
      pOut[90] = fPmQz*pIn[90] + 6*i2e*pIn[64];
      pOut[91] = fPmQz*pIn[91] + i2e*pIn[56];
      pOut[92] = fPmQz*pIn[92] + i2e*pIn[57];
      pOut[93] = fPmQz*pIn[93];
      pOut[94] = fPmQz*pIn[94] + 2*i2e*pIn[61];
      pOut[95] = fPmQz*pIn[95];
      pOut[96] = fPmQz*pIn[96] + 2*i2e*pIn[63];
      pOut[97] = fPmQz*pIn[97] + 5*i2e*pIn[68];
      pOut[98] = fPmQz*pIn[98] + 5*i2e*pIn[70];
      pOut[99] = fPmQz*pIn[99];
      pOut[100] = fPmQz*pIn[100] + 4*i2e*pIn[72];
      pOut[101] = fPmQz*pIn[101];
      pOut[102] = fPmQz*pIn[102] + 4*i2e*pIn[73];
      pOut[103] = fPmQz*pIn[103] + 3*i2e*pIn[66];
      pOut[104] = fPmQz*pIn[104] + 3*i2e*pIn[69];
      pOut[105] = fPmQz*pIn[105] + i2e*pIn[59];
      pOut[106] = fPmQz*pIn[106] + i2e*pIn[60];
      pOut[107] = fPmQz*pIn[107] + 5*i2e*pIn[74];
      pOut[108] = fPmQz*pIn[108] + 2*i2e*pIn[75];
      pOut[109] = fPmQz*pIn[109] + 4*i2e*pIn[80];
      pOut[110] = fPmQz*pIn[110] + 2*i2e*pIn[76];
      pOut[111] = fPmQz*pIn[111] + 4*i2e*pIn[82];
      pOut[112] = fPmQz*pIn[112] + i2e*pIn[65];
      pOut[113] = fPmQz*pIn[113] + i2e*pIn[67];
      pOut[114] = fPmQz*pIn[114] + i2e*pIn[71];
      pOut[115] = fPmQz*pIn[115] + 3*i2e*pIn[77];
      pOut[116] = fPmQz*pIn[116] + 3*i2e*pIn[78];
      pOut[117] = fPmQz*pIn[117] + 2*i2e*pIn[79];
      pOut[118] = fPmQz*pIn[118] + 2*i2e*pIn[81];
      pOut[119] = fPmQz*pIn[119] + 3*i2e*pIn[83];
      if (lab == 7) return;
   case 8:
      pOut[120] = fPmQz*pIn[120];
      pOut[121] = fPmQz*pIn[121];
      pOut[122] = fPmQz*pIn[122] + 8*i2e*pIn[86];
      pOut[123] = fPmQz*pIn[123];
      pOut[124] = fPmQz*pIn[124];
      pOut[125] = fPmQz*pIn[125] + i2e*pIn[84];
      pOut[126] = fPmQz*pIn[126] + 7*i2e*pIn[88];
      pOut[127] = fPmQz*pIn[127] + i2e*pIn[85];
      pOut[128] = fPmQz*pIn[128] + 7*i2e*pIn[90];
      pOut[129] = fPmQz*pIn[129];
      pOut[130] = fPmQz*pIn[130] + 2*i2e*pIn[91];
      pOut[131] = fPmQz*pIn[131];
      pOut[132] = fPmQz*pIn[132] + 6*i2e*pIn[97];
      pOut[133] = fPmQz*pIn[133] + 2*i2e*pIn[92];
      pOut[134] = fPmQz*pIn[134] + 6*i2e*pIn[98];
      pOut[135] = fPmQz*pIn[135];
      pOut[136] = fPmQz*pIn[136];
      pOut[137] = fPmQz*pIn[137] + 3*i2e*pIn[94];
      pOut[138] = fPmQz*pIn[138] + 5*i2e*pIn[100];
      pOut[139] = fPmQz*pIn[139] + 3*i2e*pIn[96];
      pOut[140] = fPmQz*pIn[140] + 5*i2e*pIn[102];
      pOut[141] = fPmQz*pIn[141];
      pOut[142] = fPmQz*pIn[142] + 4*i2e*pIn[103];
      pOut[143] = fPmQz*pIn[143] + 4*i2e*pIn[104];
      pOut[144] = fPmQz*pIn[144] + 6*i2e*pIn[107];
      pOut[145] = fPmQz*pIn[145] + i2e*pIn[87];
      pOut[146] = fPmQz*pIn[146] + i2e*pIn[89];
      pOut[147] = fPmQz*pIn[147] + 2*i2e*pIn[105];
      pOut[148] = fPmQz*pIn[148] + 2*i2e*pIn[106];
      pOut[149] = fPmQz*pIn[149] + i2e*pIn[93];
      pOut[150] = fPmQz*pIn[150] + 5*i2e*pIn[109];
      pOut[151] = fPmQz*pIn[151] + i2e*pIn[95];
      pOut[152] = fPmQz*pIn[152] + 5*i2e*pIn[111];
      pOut[153] = fPmQz*pIn[153] + 4*i2e*pIn[115];
      pOut[154] = fPmQz*pIn[154] + 4*i2e*pIn[116];
      pOut[155] = fPmQz*pIn[155] + i2e*pIn[99];
      pOut[156] = fPmQz*pIn[156] + 3*i2e*pIn[108];
      pOut[157] = fPmQz*pIn[157] + i2e*pIn[101];
      pOut[158] = fPmQz*pIn[158] + 3*i2e*pIn[110];
      pOut[159] = fPmQz*pIn[159] + 2*i2e*pIn[112];
      pOut[160] = fPmQz*pIn[160] + 2*i2e*pIn[113];
      pOut[161] = fPmQz*pIn[161] + 4*i2e*pIn[119];
      pOut[162] = fPmQz*pIn[162] + 2*i2e*pIn[114];
      pOut[163] = fPmQz*pIn[163] + 3*i2e*pIn[117];
      pOut[164] = fPmQz*pIn[164] + 3*i2e*pIn[118];
      if (lab == 8) return;
   }
   assert(0);
}

// Calculate (a0|c)^0 from (a0|0)^{m}, for a = la..lab
// Output is [nCartX(lab) - nCartX(la-1)] x (2*lc+1). Input is nCartX(lab).
// pMem must hold memory for nCartX(lab) x nCartX(lc) doubles.
// InvEtaABC = 1/(ZetaA+ZetaB+ZetaC); PmQ = (P-Q); riz = rho/zeta_c
void OsrrB_3c_shc(double *IR_RP pOut, double const *IR_RP pIn, double *IR_RP pMem, int la, unsigned lab, unsigned lc, double fPmQx, double fPmQy, double fPmQz, double InvEtaABC, double riz)
{
   unsigned nCompA = nCartX(lab);
   unsigned iCompA = nCartX(la-1);
   unsigned nCompA_Out = nCompA - iCompA;
   double i2e = .5 * InvEtaABC;
   fPmQx *= riz; fPmQy *= riz; fPmQz *= riz;
   switch(lc) {
      case 0: {
         for (unsigned i = 0; i < nCompA_Out; ++ i)
            pOut[i] = pIn[i+iCompA];
         return;
      }
      case 1: {
         double *IR_RP pOut_ = pOut - iCompA;
         OsrrB_KerM_x(&pOut_[ 0*nCompA_Out], &pIn[0], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 1*nCompA_Out], &pIn[0], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 2*nCompA_Out], &pIn[0], fPmQz, la, lab, i2e);
         return;
      }
      case 2: {
         OsrrB_KerM_x(&pMem[ 0*nCompA], &pIn[0], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[ 1*nCompA], &pIn[0], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[ 2*nCompA], &pIn[0], fPmQz, la-1, lab, i2e);
         double *IR_RP pOut_ = pMem + 3 * nCompA;
         OsrrB_KerM_x(&pOut_[ 0*nCompA_Out], &pMem[ 0*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 1*nCompA_Out], &pMem[ 1*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 2*nCompA_Out], &pMem[ 2*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 3*nCompA_Out], &pMem[ 1*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 4*nCompA_Out], &pMem[ 2*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 5*nCompA_Out], &pMem[ 2*nCompA], fPmQy, la, lab, i2e);
         ShTrC2(pOut, pOut_ + iCompA, nCompA_Out);
         return;
      }
      case 3: {
         OsrrB_KerM_x(&pMem[ 0*nCompA], &pIn[0], fPmQx, la-2, lab, i2e);
         OsrrB_KerM_y(&pMem[ 1*nCompA], &pIn[0], fPmQy, la-2, lab, i2e);
         OsrrB_KerM_z(&pMem[ 2*nCompA], &pIn[0], fPmQz, la-2, lab, i2e);
         OsrrB_KerM_x(&pMem[ 3*nCompA], &pMem[ 0*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[ 4*nCompA], &pMem[ 1*nCompA], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[ 5*nCompA], &pMem[ 2*nCompA], fPmQz, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[ 6*nCompA], &pMem[ 2*nCompA], fPmQy, la-1, lab, i2e);
         double *IR_RP pOut_ = pMem + 7 * nCompA;
         OsrrB_KerM_x(&pOut_[ 0*nCompA_Out], &pMem[ 3*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 1*nCompA_Out], &pMem[ 4*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 2*nCompA_Out], &pMem[ 5*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 3*nCompA_Out], &pMem[ 4*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 4*nCompA_Out], &pMem[ 5*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 5*nCompA_Out], &pMem[ 3*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 6*nCompA_Out], &pMem[ 5*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 7*nCompA_Out], &pMem[ 3*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 8*nCompA_Out], &pMem[ 4*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 9*nCompA_Out], &pMem[ 6*nCompA], fPmQx, la, lab, i2e);
         ShTrC3(pOut, pOut_ + iCompA, nCompA_Out);
         return;
      }
      case 4: {
         OsrrB_KerM_x(&pMem[ 0*nCompA], &pIn[0], fPmQx, la-3, lab, i2e);
         OsrrB_KerM_y(&pMem[ 1*nCompA], &pIn[0], fPmQy, la-3, lab, i2e);
         OsrrB_KerM_z(&pMem[ 2*nCompA], &pIn[0], fPmQz, la-3, lab, i2e);
         OsrrB_KerM_x(&pMem[ 3*nCompA], &pMem[ 0*nCompA], fPmQx, la-2, lab, i2e);
         OsrrB_KerM_y(&pMem[ 4*nCompA], &pMem[ 1*nCompA], fPmQy, la-2, lab, i2e);
         OsrrB_KerM_z(&pMem[ 5*nCompA], &pMem[ 2*nCompA], fPmQz, la-2, lab, i2e);
         OsrrB_KerM_x(&pMem[ 6*nCompA], &pMem[ 3*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[ 7*nCompA], &pMem[ 4*nCompA], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[ 8*nCompA], &pMem[ 5*nCompA], fPmQz, la-1, lab, i2e);
         OsrrB_KerM_x(&pMem[ 9*nCompA], &pMem[ 4*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[10*nCompA], &pMem[ 5*nCompA], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[11*nCompA], &pMem[ 3*nCompA], fPmQz, la-1, lab, i2e);
         double *IR_RP pOut_ = pMem + 12 * nCompA;
         OsrrB_KerM_x(&pOut_[ 0*nCompA_Out], &pMem[ 6*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 1*nCompA_Out], &pMem[ 7*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 2*nCompA_Out], &pMem[ 8*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 3*nCompA_Out], &pMem[ 6*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 4*nCompA_Out], &pMem[ 7*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 5*nCompA_Out], &pMem[ 6*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 6*nCompA_Out], &pMem[ 8*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 7*nCompA_Out], &pMem[ 7*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 8*nCompA_Out], &pMem[ 8*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 9*nCompA_Out], &pMem[ 9*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[10*nCompA_Out], &pMem[11*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[11*nCompA_Out], &pMem[10*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[12*nCompA_Out], &pMem[10*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[13*nCompA_Out], &pMem[ 9*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[14*nCompA_Out], &pMem[11*nCompA], fPmQy, la, lab, i2e);
         ShTrC4(pOut, pOut_ + iCompA, nCompA_Out);
         return;
      }
      case 5: {
         OsrrB_KerM_x(&pMem[ 0*nCompA], &pIn[0], fPmQx, la-4, lab, i2e);
         OsrrB_KerM_y(&pMem[ 1*nCompA], &pIn[0], fPmQy, la-4, lab, i2e);
         OsrrB_KerM_z(&pMem[ 2*nCompA], &pIn[0], fPmQz, la-4, lab, i2e);
         OsrrB_KerM_x(&pMem[ 3*nCompA], &pMem[ 0*nCompA], fPmQx, la-3, lab, i2e);
         OsrrB_KerM_y(&pMem[ 4*nCompA], &pMem[ 1*nCompA], fPmQy, la-3, lab, i2e);
         OsrrB_KerM_z(&pMem[ 5*nCompA], &pMem[ 2*nCompA], fPmQz, la-3, lab, i2e);
         OsrrB_KerM_x(&pMem[ 6*nCompA], &pMem[ 3*nCompA], fPmQx, la-2, lab, i2e);
         OsrrB_KerM_y(&pMem[ 7*nCompA], &pMem[ 4*nCompA], fPmQy, la-2, lab, i2e);
         OsrrB_KerM_z(&pMem[ 8*nCompA], &pMem[ 5*nCompA], fPmQz, la-2, lab, i2e);
         OsrrB_KerM_x(&pMem[ 9*nCompA], &pMem[ 4*nCompA], fPmQx, la-2, lab, i2e);
         OsrrB_KerM_y(&pMem[10*nCompA], &pMem[ 5*nCompA], fPmQy, la-2, lab, i2e);
         OsrrB_KerM_x(&pMem[11*nCompA], &pMem[ 6*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[12*nCompA], &pMem[ 7*nCompA], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[13*nCompA], &pMem[ 8*nCompA], fPmQz, la-1, lab, i2e);
         OsrrB_KerM_x(&pMem[14*nCompA], &pMem[ 7*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_z(&pMem[15*nCompA], &pMem[ 6*nCompA], fPmQz, la-1, lab, i2e);
         OsrrB_KerM_x(&pMem[16*nCompA], &pMem[ 8*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_x(&pMem[17*nCompA], &pMem[ 9*nCompA], fPmQx, la-1, lab, i2e);
         OsrrB_KerM_y(&pMem[18*nCompA], &pMem[10*nCompA], fPmQy, la-1, lab, i2e);
         OsrrB_KerM_x(&pMem[19*nCompA], &pMem[10*nCompA], fPmQx, la-1, lab, i2e);
         double *IR_RP pOut_ = pMem + 20 * nCompA;
         OsrrB_KerM_x(&pOut_[ 0*nCompA_Out], &pMem[11*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 1*nCompA_Out], &pMem[12*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 2*nCompA_Out], &pMem[13*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 3*nCompA_Out], &pMem[12*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 4*nCompA_Out], &pMem[13*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 5*nCompA_Out], &pMem[11*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[ 6*nCompA_Out], &pMem[13*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 7*nCompA_Out], &pMem[15*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[ 8*nCompA_Out], &pMem[12*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[ 9*nCompA_Out], &pMem[17*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[10*nCompA_Out], &pMem[15*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[11*nCompA_Out], &pMem[14*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[12*nCompA_Out], &pMem[18*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[13*nCompA_Out], &pMem[16*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[14*nCompA_Out], &pMem[18*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_y(&pOut_[15*nCompA_Out], &pMem[15*nCompA], fPmQy, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[16*nCompA_Out], &pMem[14*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[17*nCompA_Out], &pMem[19*nCompA], fPmQz, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[18*nCompA_Out], &pMem[18*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_x(&pOut_[19*nCompA_Out], &pMem[19*nCompA], fPmQx, la, lab, i2e);
         OsrrB_KerM_z(&pOut_[20*nCompA_Out], &pMem[17*nCompA], fPmQz, la, lab, i2e);
         ShTrC5(pOut, pOut_ + iCompA, nCompA_Out);
         return;
      }
   }
   assert(0);
}

// indices for factoring [CartX(la+lb) - nCartX(la)] into [CartX(lb)] x [CartY(la)]
// format: nCartY(la) x CartX(lb)
// Component order CartX(lab=0): s
static const cart_index_t iCartXY_ab0_a0[1] = {
   0
}; // 0.00 kb
// Component order CartX(lab=1): s x y z
static const cart_index_t iCartXY_ab1_a1[3] = {
   0, 1, 2
}; // 0.01 kb
// Component order CartX(lab=2): s x y z xx yy zz xy xz yz
static const cart_index_t iCartXY_ab2_a1[12] = {
   0, 1, 2, 3, 6, 7, 6, 4, 8, 7, 8, 5
}; // 0.02 kb
static const cart_index_t iCartXY_ab2_a2[6] = {
   0, 1, 2, 3, 4, 5
}; // 0.01 kb
// Component order CartX(lab=3): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz
static const cart_index_t iCartXY_ab3_a2[24] = {
   0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 13, 15, 11, 7, 12, 9, 15, 14, 13, 14, 8, 15, 10, 12
}; // 0.05 kb
static const cart_index_t iCartXY_ab3_a3[10] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9
}; // 0.02 kb
// Component order CartX(lab=4): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz xxxx yyyy zzzz xxxy xyyy xxxz xzzz yyyz yzzz xxyy xxzz yyzz xyzz xyyz xxyz
static const cart_index_t iCartXY_ab4_a2[60] = {
   0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 13, 15, 11, 7, 12, 9, 15, 14, 13, 14, 8, 15, 10, 12, 16, 25,
   26, 19, 21, 30, 25, 17, 27, 20, 29, 23, 26, 27, 18, 28, 22, 24, 19, 20, 28, 25, 30, 29, 21, 29, 22, 30,
   26, 28, 30, 23, 24, 29, 28, 27
}; // 0.12 kb
static const cart_index_t iCartXY_ab4_a3[40] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 16, 19, 20, 13, 22, 15, 23, 24, 13, 11, 18, 14, 22, 19,
   21, 24, 17, 23, 15, 17, 12, 23, 16, 24, 18, 20, 21, 22
}; // 0.08 kb
static const cart_index_t iCartXY_ab4_a4[15] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
}; // 0.03 kb
// Component order CartX(lab=5): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz xxxx yyyy zzzz xxxy xyyy xxxz xzzz yyyz yzzz xxyy xxzz yyzz xyzz xyyz xxyz xxxxx yyyyy zzzzz xyyyy xzzzz xxxxy yzzzz xxxxz yyyyz xxxyy xxxzz xxyyy yyyzz xxzzz yyzzz xxxyz xyyyz xyzzz xyyzz xxyzz xxyyz
static const cart_index_t iCartXY_ab5_a3[100] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 16, 19, 20, 13, 22, 15, 23, 24, 13, 11, 18, 14, 22, 19,
   21, 24, 17, 23, 15, 17, 12, 23, 16, 24, 18, 20, 21, 22, 25, 36, 38, 34, 35, 30, 44, 32, 45, 40, 34, 26,
   39, 28, 43, 36, 37, 45, 33, 41, 35, 37, 27, 43, 29, 44, 31, 38, 39, 42, 30, 28, 42, 36, 44, 34, 43, 40,
   41, 45, 32, 41, 29, 45, 38, 40, 42, 35, 43, 44, 40, 33, 31, 41, 42, 45, 39, 44, 37, 43
}; // 0.20 kb
static const cart_index_t iCartXY_ab5_a4[60] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 26, 22, 28, 31, 32, 24, 25,
   33, 34, 35, 30, 20, 16, 21, 24, 18, 30, 32, 23, 29, 26, 34, 27, 33, 31, 35, 22, 23, 17, 30, 31, 25, 19,
   27, 21, 35, 28, 29, 32, 33, 34
}; // 0.12 kb
static const cart_index_t iCartXY_ab5_a5[21] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
}; // 0.04 kb
// Component order CartX(lab=6): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz xxxx yyyy zzzz xxxy xyyy xxxz xzzz yyyz yzzz xxyy xxzz yyzz xyzz xyyz xxyz xxxxx yyyyy zzzzz xyyyy xzzzz xxxxy yzzzz xxxxz yyyyz xxxyy xxxzz xxyyy yyyzz xxzzz yyzzz xxxyz xyyyz xyzzz xyyzz xxyzz xxyyz xxxxxx yyyyyy zzzzzz xxxxxy xyyyyy xxxxxz xzzzzz yyyyyz yzzzzz xxxxyy xxxxzz xxyyyy xxzzzz yyyyzz yyzzzz xxxyyy xxxzzz yyyzzz xyzzzz xyyyyz xxxxyz xxxyzz xyyyzz xxxyyz xyyzzz xxyyyz xxyzzz xxyyzz
static const cart_index_t iCartXY_ab6_a3[200] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 16, 19, 20, 13, 22, 15, 23, 24, 13, 11, 18, 14, 22, 19,
   21, 24, 17, 23, 15, 17, 12, 23, 16, 24, 18, 20, 21, 22, 25, 36, 38, 34, 35, 30, 44, 32, 45, 40, 34, 26,
   39, 28, 43, 36, 37, 45, 33, 41, 35, 37, 27, 43, 29, 44, 31, 38, 39, 42, 30, 28, 42, 36, 44, 34, 43, 40,
   41, 45, 32, 41, 29, 45, 38, 40, 42, 35, 43, 44, 40, 33, 31, 41, 42, 45, 39, 44, 37, 43, 46, 61, 62, 55,
   56, 49, 67, 51, 69, 66, 61, 47, 63, 50, 68, 57, 59, 71, 53, 65, 62, 63, 48, 70, 52, 72, 54, 58, 60, 64,
   55, 50, 70, 57, 73, 61, 68, 69, 65, 71, 56, 68, 52, 73, 58, 67, 64, 62, 70, 72, 49, 57, 72, 61, 67, 55,
   73, 66, 71, 69, 67, 59, 54, 68, 64, 73, 60, 72, 63, 70, 51, 71, 58, 69, 62, 66, 72, 56, 73, 67, 69, 53,
   60, 65, 70, 71, 63, 73, 59, 68, 66, 65, 64, 71, 72, 69, 70, 67, 68, 73
}; // 0.39 kb
static const cart_index_t iCartXY_ab6_a4[150] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 26, 22, 28, 31, 32, 24, 25,
   33, 34, 35, 30, 20, 16, 21, 24, 18, 30, 32, 23, 29, 26, 34, 27, 33, 31, 35, 22, 23, 17, 30, 31, 25, 19,
   27, 21, 35, 28, 29, 32, 33, 34, 36, 47, 48, 39, 51, 41, 52, 61, 62, 45, 46, 63, 57, 59, 56, 45, 37, 50,
   51, 40, 59, 60, 43, 53, 47, 63, 49, 58, 55, 61, 46, 49, 38, 57, 58, 52, 42, 53, 44, 63, 48, 50, 54, 60,
   62, 39, 40, 54, 45, 47, 56, 62, 55, 60, 51, 57, 58, 63, 61, 59, 41, 55, 42, 56, 61, 46, 48, 58, 54, 59,
   52, 60, 62, 63, 57, 56, 43, 44, 59, 55, 57, 54, 49, 50, 61, 62, 53, 60, 58, 63
}; // 0.29 kb
static const cart_index_t iCartXY_ab6_a5[84] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 27, 32, 33,
   24, 39, 26, 40, 30, 31, 36, 43, 37, 45, 41, 46, 47, 48, 42, 44, 24, 22, 29, 25, 39, 30, 35, 41, 28, 36,
   42, 32, 34, 47, 38, 44, 40, 45, 43, 48, 46, 26, 28, 23, 40, 27, 41, 29, 31, 34, 44, 37, 46, 38, 33, 35,
   42, 43, 39, 45, 47, 48
}; // 0.16 kb
static const cart_index_t iCartXY_ab6_a6[28] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27
}; // 0.05 kb
// Component order CartX(lab=7): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz xxxx yyyy zzzz xxxy xyyy xxxz xzzz yyyz yzzz xxyy xxzz yyzz xyzz xyyz xxyz xxxxx yyyyy zzzzz xyyyy xzzzz xxxxy yzzzz xxxxz yyyyz xxxyy xxxzz xxyyy yyyzz xxzzz yyzzz xxxyz xyyyz xyzzz xyyzz xxyzz xxyyz xxxxxx yyyyyy zzzzzz xxxxxy xyyyyy xxxxxz xzzzzz yyyyyz yzzzzz xxxxyy xxxxzz xxyyyy xxzzzz yyyyzz yyzzzz xxxyyy xxxzzz yyyzzz xyzzzz xyyyyz xxxxyz xxxyzz xyyyzz xxxyyz xyyzzz xxyyyz xxyzzz xxyyzz xxxxxxx yyyyyyy zzzzzzz xyyyyyy xzzzzzz xxxxxxy yzzzzzz xxxxxxz yyyyyyz xxxxxyy xxxxxzz xxyyyyy yyyyyzz xxzzzzz yyzzzzz xxxyyyy xxxzzzz xxxxyyy yyyzzzz xxxxzzz yyyyzzz xxxxxyz xyyyyyz xyzzzzz xyyyyzz xyyzzzz xxxxyzz xxyzzzz xxxxyyz xxyyyyz xxxyyyz xxxyzzz xyyyzzz xxxyyzz xxyyyzz xxyyzzz
static const cart_index_t iCartXY_ab7_a4[300] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 26, 22, 28, 31, 32, 24, 25,
   33, 34, 35, 30, 20, 16, 21, 24, 18, 30, 32, 23, 29, 26, 34, 27, 33, 31, 35, 22, 23, 17, 30, 31, 25, 19,
   27, 21, 35, 28, 29, 32, 33, 34, 36, 47, 48, 39, 51, 41, 52, 61, 62, 45, 46, 63, 57, 59, 56, 45, 37, 50,
   51, 40, 59, 60, 43, 53, 47, 63, 49, 58, 55, 61, 46, 49, 38, 57, 58, 52, 42, 53, 44, 63, 48, 50, 54, 60,
   62, 39, 40, 54, 45, 47, 56, 62, 55, 60, 51, 57, 58, 63, 61, 59, 41, 55, 42, 56, 61, 46, 48, 58, 54, 59,
   52, 60, 62, 63, 57, 56, 43, 44, 59, 55, 57, 54, 49, 50, 61, 62, 53, 60, 58, 63, 64, 79, 80, 69, 81, 71,
   83, 94, 95, 73, 74, 97, 90, 92, 85, 81, 65, 82, 79, 67, 94, 96, 72, 84, 75, 98, 76, 88, 86, 93, 83, 84,
   66, 95, 96, 80, 68, 82, 70, 99, 77, 78, 87, 89, 91, 73, 67, 89, 81, 75, 92, 99, 86, 96, 79, 97, 88, 98,
   93, 94, 74, 88, 68, 90, 98, 83, 77, 96, 87, 97, 80, 89, 91, 99, 95, 69, 75, 91, 73, 79, 85, 95, 93, 99,
   81, 90, 98, 97, 94, 92, 90, 76, 70, 97, 88, 95, 87, 84, 78, 98, 91, 82, 89, 96, 99, 71, 93, 77, 85, 94,
   74, 80, 98, 91, 92, 83, 99, 95, 97, 90, 92, 72, 78, 94, 86, 97, 89, 76, 82, 93, 99, 84, 96, 88, 98, 85,
   86, 87, 92, 93, 90, 91, 88, 89, 94, 95, 96, 99, 98, 97
}; // 0.59 kb
static const cart_index_t iCartXY_ab7_a5[210] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 27, 32, 33,
   24, 39, 26, 40, 30, 31, 36, 43, 37, 45, 41, 46, 47, 48, 42, 44, 24, 22, 29, 25, 39, 30, 35, 41, 28, 36,
   42, 32, 34, 47, 38, 44, 40, 45, 43, 48, 46, 26, 28, 23, 40, 27, 41, 29, 31, 34, 44, 37, 46, 38, 33, 35,
   42, 43, 39, 45, 47, 48, 49, 60, 62, 64, 65, 54, 76, 56, 78, 58, 59, 66, 83, 68, 84, 70, 79, 80, 82, 75,
   77, 58, 50, 63, 52, 74, 66, 67, 77, 57, 64, 82, 60, 61, 84, 69, 79, 71, 81, 73, 83, 78, 59, 61, 51, 73,
   53, 75, 55, 68, 69, 82, 65, 83, 67, 62, 63, 80, 81, 72, 74, 76, 84, 54, 52, 72, 60, 76, 58, 74, 70, 71,
   66, 75, 64, 73, 80, 81, 77, 78, 84, 83, 82, 79, 56, 71, 53, 78, 62, 70, 72, 59, 73, 77, 68, 79, 81, 65,
   74, 75, 83, 76, 84, 80, 82, 70, 57, 55, 71, 72, 77, 63, 75, 61, 79, 80, 78, 69, 76, 67, 82, 73, 74, 81,
   84, 83
}; // 0.41 kb
static const cart_index_t iCartXY_ab7_a6[112] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 31, 32, 33, 39, 35, 41, 50, 51, 37, 38, 43, 44, 52, 53, 45, 47, 60, 55, 57, 49, 54, 62, 56,
   63, 58, 59, 61, 33, 29, 34, 37, 31, 49, 51, 36, 42, 45, 54, 39, 55, 40, 46, 43, 59, 48, 53, 50, 56, 61,
   52, 58, 60, 57, 63, 62, 35, 36, 30, 49, 50, 38, 32, 40, 34, 56, 47, 57, 41, 48, 42, 58, 44, 46, 51, 52,
   54, 59, 60, 61, 53, 62, 55, 63
}; // 0.22 kb
static const cart_index_t iCartXY_ab7_a7[36] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35
}; // 0.07 kb
// Component order CartX(lab=8): s x y z xx yy zz xy xz yz xxx yyy zzz xyy xzz xxy yzz xxz yyz xyz xxxx yyyy zzzz xxxy xyyy xxxz xzzz yyyz yzzz xxyy xxzz yyzz xyzz xyyz xxyz xxxxx yyyyy zzzzz xyyyy xzzzz xxxxy yzzzz xxxxz yyyyz xxxyy xxxzz xxyyy yyyzz xxzzz yyzzz xxxyz xyyyz xyzzz xyyzz xxyzz xxyyz xxxxxx yyyyyy zzzzzz xxxxxy xyyyyy xxxxxz xzzzzz yyyyyz yzzzzz xxxxyy xxxxzz xxyyyy xxzzzz yyyyzz yyzzzz xxxyyy xxxzzz yyyzzz xyzzzz xyyyyz xxxxyz xxxyzz xyyyzz xxxyyz xyyzzz xxyyyz xxyzzz xxyyzz xxxxxxx yyyyyyy zzzzzzz xyyyyyy xzzzzzz xxxxxxy yzzzzzz xxxxxxz yyyyyyz xxxxxyy xxxxxzz xxyyyyy yyyyyzz xxzzzzz yyzzzzz xxxyyyy xxxzzzz xxxxyyy yyyzzzz xxxxzzz yyyyzzz xxxxxyz xyyyyyz xyzzzzz xyyyyzz xyyzzzz xxxxyzz xxyzzzz xxxxyyz xxyyyyz xxxyyyz xxxyzzz xyyyzzz xxxyyzz xxyyyzz xxyyzzz xxxxxxxx yyyyyyyy zzzzzzzz xxxxxxxy xyyyyyyy xxxxxxxz xzzzzzzz yyyyyyyz yzzzzzzz xxxxxxyy xxxxxxzz xxyyyyyy xxzzzzzz yyyyyyzz yyzzzzzz xxxxxyyy xxxyyyyy xxxxxzzz xxxzzzzz yyyyyzzz yyyzzzzz xxxxyyyy xxxxzzzz yyyyzzzz xyzzzzzz xyyyyyyz xxxxxxyz xxxxxyzz xyyyyyzz xxxxxyyz xyyzzzzz xxyyyyyz xxyzzzzz xxxyzzzz xyyyzzzz xxxyyyyz xyyyyzzz xxxxyyyz xxxxyzzz xxxxyyzz xxyyyyzz xxyyzzzz xxxyyyzz xxxyyzzz xxyyyzzz
static const cart_index_t iCartXY_ab8_a4[525] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 26, 22, 28, 31, 32, 24, 25,
   33, 34, 35, 30, 20, 16, 21, 24, 18, 30, 32, 23, 29, 26, 34, 27, 33, 31, 35, 22, 23, 17, 30, 31, 25, 19,
   27, 21, 35, 28, 29, 32, 33, 34, 36, 47, 48, 39, 51, 41, 52, 61, 62, 45, 46, 63, 57, 59, 56, 45, 37, 50,
   51, 40, 59, 60, 43, 53, 47, 63, 49, 58, 55, 61, 46, 49, 38, 57, 58, 52, 42, 53, 44, 63, 48, 50, 54, 60,
   62, 39, 40, 54, 45, 47, 56, 62, 55, 60, 51, 57, 58, 63, 61, 59, 41, 55, 42, 56, 61, 46, 48, 58, 54, 59,
   52, 60, 62, 63, 57, 56, 43, 44, 59, 55, 57, 54, 49, 50, 61, 62, 53, 60, 58, 63, 64, 79, 80, 69, 81, 71,
   83, 94, 95, 73, 74, 97, 90, 92, 85, 81, 65, 82, 79, 67, 94, 96, 72, 84, 75, 98, 76, 88, 86, 93, 83, 84,
   66, 95, 96, 80, 68, 82, 70, 99, 77, 78, 87, 89, 91, 73, 67, 89, 81, 75, 92, 99, 86, 96, 79, 97, 88, 98,
   93, 94, 74, 88, 68, 90, 98, 83, 77, 96, 87, 97, 80, 89, 91, 99, 95, 69, 75, 91, 73, 79, 85, 95, 93, 99,
   81, 90, 98, 97, 94, 92, 90, 76, 70, 97, 88, 95, 87, 84, 78, 98, 91, 82, 89, 96, 99, 71, 93, 77, 85, 94,
   74, 80, 98, 91, 92, 83, 99, 95, 97, 90, 92, 72, 78, 94, 86, 97, 89, 76, 82, 93, 99, 84, 96, 88, 98, 85,
   86, 87, 92, 93, 90, 91, 88, 89, 94, 95, 96, 99, 98, 97, 100, 121, 122, 103, 115, 105, 117, 137, 138, 109, 110, 139,
   127, 129, 126, 121, 101, 123, 116, 104, 135, 136, 107, 119, 111, 140, 113, 128, 125, 131, 122, 123, 102, 133, 134, 118, 106, 120,
   108, 141, 112, 114, 124, 130, 132, 103, 116, 133, 109, 121, 126, 138, 135, 143, 115, 127, 142, 139, 137, 129, 115, 104, 134, 121,
   111, 137, 144, 125, 136, 116, 142, 128, 140, 131, 135, 105, 135, 118, 126, 137, 110, 122, 142, 133, 129, 117, 143, 138, 139, 127,
   117, 136, 106, 138, 144, 122, 112, 134, 124, 143, 118, 130, 132, 141, 133, 137, 107, 120, 135, 125, 142, 134, 113, 123, 131, 144,
   119, 136, 128, 140, 138, 119, 108, 143, 136, 133, 124, 123, 114, 144, 132, 120, 130, 134, 141, 109, 111, 141, 115, 116, 129, 143,
   131, 144, 121, 139, 140, 142, 135, 137, 110, 140, 112, 127, 142, 117, 118, 144, 132, 139, 122, 141, 133, 143, 138, 139, 113, 114,
   142, 128, 143, 130, 119, 120, 140, 141, 123, 134, 136, 144, 127, 128, 124, 139, 140, 138, 132, 136, 130, 142, 133, 134, 141, 144,
   143, 129, 125, 130, 137, 131, 139, 141, 128, 134, 135, 143, 136, 144, 140, 142, 126, 131, 132, 129, 135, 127, 133, 140, 141, 137,
   138, 144, 143, 142, 139
}; // 1.03 kb
static const cart_index_t iCartXY_ab8_a5[420] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 27, 32, 33,
   24, 39, 26, 40, 30, 31, 36, 43, 37, 45, 41, 46, 47, 48, 42, 44, 24, 22, 29, 25, 39, 30, 35, 41, 28, 36,
   42, 32, 34, 47, 38, 44, 40, 45, 43, 48, 46, 26, 28, 23, 40, 27, 41, 29, 31, 34, 44, 37, 46, 38, 33, 35,
   42, 43, 39, 45, 47, 48, 49, 60, 62, 64, 65, 54, 76, 56, 78, 58, 59, 66, 83, 68, 84, 70, 79, 80, 82, 75,
   77, 58, 50, 63, 52, 74, 66, 67, 77, 57, 64, 82, 60, 61, 84, 69, 79, 71, 81, 73, 83, 78, 59, 61, 51, 73,
   53, 75, 55, 68, 69, 82, 65, 83, 67, 62, 63, 80, 81, 72, 74, 76, 84, 54, 52, 72, 60, 76, 58, 74, 70, 71,
   66, 75, 64, 73, 80, 81, 77, 78, 84, 83, 82, 79, 56, 71, 53, 78, 62, 70, 72, 59, 73, 77, 68, 79, 81, 65,
   74, 75, 83, 76, 84, 80, 82, 70, 57, 55, 71, 72, 77, 63, 75, 61, 79, 80, 78, 69, 76, 67, 82, 73, 74, 81,
   84, 83, 85, 101, 103, 106, 107, 88, 118, 90, 120, 94, 95, 100, 127, 102, 128, 111, 122, 123, 124, 112, 114, 100, 86, 105,
   89, 119, 106, 108, 122, 92, 101, 127, 96, 98, 129, 104, 120, 110, 121, 113, 125, 116, 102, 104, 87, 121, 91, 123, 93, 107,
   108, 128, 103, 129, 105, 97, 99, 118, 119, 109, 115, 117, 126, 94, 89, 115, 96, 126, 100, 119, 114, 110, 106, 124, 101, 113,
   128, 121, 122, 116, 129, 125, 127, 120, 95, 113, 91, 125, 97, 112, 109, 102, 121, 124, 107, 127, 119, 103, 115, 123, 129, 117,
   126, 118, 128, 88, 96, 117, 101, 118, 94, 126, 111, 116, 100, 112, 106, 125, 123, 129, 114, 120, 128, 127, 124, 122, 112, 98,
   93, 113, 109, 124, 99, 123, 104, 127, 118, 125, 108, 117, 105, 128, 121, 115, 119, 126, 129, 90, 116, 97, 120, 103, 111, 117,
   95, 125, 114, 102, 122, 129, 107, 126, 112, 127, 118, 128, 123, 124, 114, 92, 99, 110, 115, 122, 105, 124, 98, 120, 128, 116,
   104, 126, 108, 127, 113, 119, 121, 129, 125, 111, 110, 109, 116, 117, 114, 115, 112, 113, 122, 123, 120, 121, 118, 119, 124, 125,
   126, 129, 128, 127
}; // 0.82 kb
static const cart_index_t iCartXY_ab8_a6[280] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 31, 32, 33, 39, 35, 41, 50, 51, 37, 38, 43, 44, 52, 53, 45, 47, 60, 55, 57, 49, 54, 62, 56,
   63, 58, 59, 61, 33, 29, 34, 37, 31, 49, 51, 36, 42, 45, 54, 39, 55, 40, 46, 43, 59, 48, 53, 50, 56, 61,
   52, 58, 60, 57, 63, 62, 35, 36, 30, 49, 50, 38, 32, 40, 34, 56, 47, 57, 41, 48, 42, 58, 44, 46, 51, 52,
   54, 59, 60, 61, 53, 62, 55, 63, 64, 75, 76, 67, 80, 69, 82, 95, 96, 73, 74, 85, 86, 104, 105, 79, 81, 108,
   97, 99, 90, 91, 106, 93, 107, 101, 102, 103, 73, 65, 78, 79, 68, 93, 94, 71, 84, 85, 103, 75, 105, 77, 87, 80,
   107, 83, 98, 89, 101, 106, 92, 99, 100, 95, 108, 104, 74, 77, 66, 91, 92, 81, 70, 83, 72, 103, 86, 104, 76, 87,
   78, 106, 82, 84, 88, 100, 102, 97, 98, 107, 94, 108, 96, 105, 67, 68, 88, 73, 75, 90, 96, 89, 94, 79, 91, 80,
   97, 92, 98, 85, 102, 100, 105, 95, 93, 103, 104, 101, 108, 99, 107, 106, 69, 89, 70, 90, 95, 74, 76, 92, 88, 93,
   81, 99, 82, 100, 94, 101, 86, 98, 96, 104, 91, 102, 108, 103, 105, 106, 97, 107, 90, 71, 72, 93, 89, 91, 88, 77,
   78, 101, 102, 95, 96, 83, 84, 99, 97, 87, 94, 92, 103, 107, 100, 106, 98, 104, 105, 108
}; // 0.55 kb
static const cart_index_t iCartXY_ab8_a7[144] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 40, 42, 47, 48, 39, 60, 41, 61, 45, 46, 52, 64, 54, 66, 57,
   58, 51, 70, 53, 72, 62, 67, 68, 76, 77, 63, 69, 65, 71, 73, 74, 80, 75, 78, 79, 39, 37, 44, 40, 60, 45,
   50, 62, 43, 51, 63, 47, 49, 68, 56, 52, 69, 57, 59, 74, 55, 65, 61, 66, 64, 70, 75, 77, 73, 67, 71, 79,
   72, 78, 76, 80, 41, 43, 38, 61, 42, 62, 44, 46, 49, 65, 53, 67, 55, 48, 50, 71, 54, 73, 56, 58, 59, 63,
   64, 60, 72, 66, 74, 68, 75, 76, 78, 69, 70, 79, 80, 77
}; // 0.28 kb
static const cart_index_t iCartXY_ab8_a8[45] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44
}; // 0.09 kb

// indices for factoring [CartY(la+lb)] into [CartY(lb)] x [CartY(la)]
// format: nCartY(la) x CartY(lb)
static const cart_index_t iCartYY_ab0_a0[1] = {
   0
}; // 0.00 kb
static const cart_index_t iCartYY_ab1_a0[3] = {
   0, 1, 2
}; // 0.01 kb
static const cart_index_t iCartYY_ab1_a1[3] = {
   0, 1, 2
}; // 0.01 kb
static const cart_index_t iCartYY_ab2_a0[6] = {
   0, 1, 2, 3, 4, 5
}; // 0.01 kb
static const cart_index_t iCartYY_ab2_a1[9] = {
   0, 3, 4, 3, 1, 5, 4, 5, 2
}; // 0.02 kb
static const cart_index_t iCartYY_ab2_a2[6] = {
   0, 1, 2, 3, 4, 5
}; // 0.01 kb
static const cart_index_t iCartYY_ab3_a0[10] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9
}; // 0.02 kb
static const cart_index_t iCartYY_ab3_a1[18] = {
   0, 5, 7, 3, 1, 8, 4, 6, 2, 5, 3, 9, 7, 9, 4, 9, 8, 6
}; // 0.04 kb
static const cart_index_t iCartYY_ab3_a2[18] = {
   0, 3, 4, 5, 7, 9, 5, 1, 6, 3, 9, 8, 7, 8, 2, 9, 4, 6
}; // 0.04 kb
static const cart_index_t iCartYY_ab3_a3[10] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9
}; // 0.02 kb
static const cart_index_t iCartYY_ab4_a0[15] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
}; // 0.03 kb
static const cart_index_t iCartYY_ab4_a1[30] = {
   0, 3, 5, 4, 1, 7, 6, 8, 2, 9, 4, 13, 10, 12, 6, 3, 9, 14, 12, 11, 8, 5, 14, 10, 13, 7,
   11, 14, 13, 12
}; // 0.06 kb
static const cart_index_t iCartYY_ab4_a2[36] = {
   0, 9, 10, 3, 5, 14, 9, 1, 11, 4, 13, 7, 10, 11, 2, 12, 6, 8, 3, 4, 12, 9, 14, 13, 5, 13,
   6, 14, 10, 12, 14, 7, 8, 13, 12, 11
}; // 0.07 kb
static const cart_index_t iCartYY_ab4_a3[30] = {
   0, 4, 6, 9, 10, 3, 12, 5, 13, 14, 3, 1, 8, 4, 12, 9, 11, 14, 7, 13, 5, 7, 2, 13, 6, 14,
   8, 10, 11, 12
}; // 0.06 kb
static const cart_index_t iCartYY_ab4_a4[15] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
}; // 0.03 kb
static const cart_index_t iCartYY_ab5_a0[21] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
}; // 0.04 kb
static const cart_index_t iCartYY_ab5_a1[45] = {
   0, 5, 7, 3, 1, 8, 4, 6, 2, 5, 9, 15, 11, 3, 16, 7, 15, 10, 13, 17, 4, 16, 8, 12, 17, 14,
   6, 9, 11, 20, 10, 19, 13, 18, 12, 14, 19, 18, 17, 20, 16, 18, 15, 20, 19
}; // 0.09 kb
static const cart_index_t iCartYY_ab5_a2[60] = {
   0, 9, 10, 5, 7, 15, 11, 1, 12, 3, 16, 8, 13, 14, 2, 17, 4, 6, 9, 3, 18, 11, 20, 16, 10, 18,
   4, 19, 13, 17, 5, 11, 19, 9, 15, 20, 19, 12, 6, 18, 17, 14, 7, 20, 13, 15, 10, 19, 20, 8, 14, 16,
   18, 12, 15, 16, 17, 20, 19, 18
}; // 0.12 kb
static const cart_index_t iCartYY_ab5_a3[60] = {
   0, 11, 13, 9, 10, 5, 19, 7, 20, 15, 9, 1, 14, 3, 18, 11, 12, 20, 8, 16, 10, 12, 2, 18, 4, 19,
   6, 13, 14, 17, 5, 3, 17, 11, 19, 9, 18, 15, 16, 20, 7, 16, 4, 20, 13, 15, 17, 10, 18, 19, 15, 8,
   6, 16, 17, 20, 14, 19, 12, 18
}; // 0.12 kb
static const cart_index_t iCartYY_ab5_a4[45] = {
   0, 3, 4, 5, 11, 7, 13, 16, 17, 9, 10, 18, 19, 20, 15, 5, 1, 6, 9, 3, 15, 17, 8, 14, 11, 19,
   12, 18, 16, 20, 7, 8, 2, 15, 16, 10, 4, 12, 6, 20, 13, 14, 17, 18, 19
}; // 0.09 kb
static const cart_index_t iCartYY_ab5_a5[21] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
}; // 0.04 kb
static const cart_index_t iCartYY_ab6_a0[28] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27
}; // 0.05 kb
static const cart_index_t iCartYY_ab6_a1[63] = {
   0, 3, 5, 4, 1, 7, 6, 8, 2, 11, 4, 19, 12, 18, 6, 3, 9, 20, 18, 14, 8, 5, 20, 10, 19, 7,
   13, 9, 15, 23, 10, 21, 16, 15, 11, 25, 22, 13, 17, 16, 26, 12, 24, 17, 14, 20, 23, 21, 25, 19, 22, 26,
   24, 18, 27, 22, 24, 21, 27, 26, 23, 25, 27
}; // 0.12 kb
static const cart_index_t iCartYY_ab6_a2[90] = {
   0, 9, 10, 3, 5, 20, 11, 1, 13, 4, 19, 7, 12, 14, 2, 18, 6, 8, 3, 15, 21, 9, 20, 23, 15, 4,
   22, 11, 25, 19, 5, 23, 16, 20, 10, 21, 16, 24, 6, 26, 12, 18, 25, 7, 17, 19, 22, 13, 26, 17, 8, 24,
   18, 14, 9, 11, 27, 15, 23, 25, 10, 27, 12, 21, 16, 26, 27, 13, 14, 22, 24, 17, 21, 22, 18, 27, 26, 24,
   23, 19, 24, 25, 27, 22, 20, 25, 26, 23, 21, 27
}; // 0.18 kb
static const cart_index_t iCartYY_ab6_a3[100] = {
   0, 15, 16, 9, 10, 3, 21, 5, 23, 20, 15, 1, 17, 4, 22, 11, 13, 25, 7, 19, 16, 17, 2, 24, 6, 26,
   8, 12, 14, 18, 9, 4, 24, 11, 27, 15, 22, 23, 19, 25, 10, 22, 6, 27, 12, 21, 18, 16, 24, 26, 3, 11,
   26, 15, 21, 9, 27, 20, 25, 23, 21, 13, 8, 22, 18, 27, 14, 26, 17, 24, 5, 25, 12, 23, 16, 20, 26, 10,
   27, 21, 23, 7, 14, 19, 24, 25, 17, 27, 13, 22, 20, 19, 18, 25, 26, 23, 24, 21, 22, 27
}; // 0.20 kb
static const cart_index_t iCartYY_ab6_a4[90] = {
   0, 11, 12, 3, 15, 5, 16, 25, 26, 9, 10, 27, 21, 23, 20, 9, 1, 14, 15, 4, 23, 24, 7, 17, 11, 27,
   13, 22, 19, 25, 10, 13, 2, 21, 22, 16, 6, 17, 8, 27, 12, 14, 18, 24, 26, 3, 4, 18, 9, 11, 20, 26,
   19, 24, 15, 21, 22, 27, 25, 23, 5, 19, 6, 20, 25, 10, 12, 22, 18, 23, 16, 24, 26, 27, 21, 20, 7, 8,
   23, 19, 21, 18, 13, 14, 25, 26, 17, 24, 22, 27
}; // 0.18 kb
static const cart_index_t iCartYY_ab6_a5[63] = {
   0, 4, 6, 11, 12, 3, 18, 5, 19, 9, 10, 15, 22, 16, 24, 20, 25, 26, 27, 21, 23, 3, 1, 8, 4, 18,
   9, 14, 20, 7, 15, 21, 11, 13, 26, 17, 23, 19, 24, 22, 27, 25, 5, 7, 2, 19, 6, 20, 8, 10, 13, 23,
   16, 25, 17, 12, 14, 21, 22, 18, 24, 26, 27
}; // 0.12 kb
static const cart_index_t iCartYY_ab6_a6[28] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27
}; // 0.05 kb
static const cart_index_t iCartYY_ab7_a0[36] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35
}; // 0.07 kb
static const cart_index_t iCartYY_ab7_a1[84] = {
   0, 5, 7, 3, 1, 8, 4, 6, 2, 5, 9, 21, 11, 3, 22, 7, 21, 10, 13, 23, 4, 22, 8, 12, 23, 14,
   6, 9, 17, 28, 10, 26, 19, 15, 11, 29, 16, 27, 13, 24, 12, 20, 25, 18, 14, 17, 15, 30, 19, 31, 16, 32,
   20, 18, 27, 25, 23, 29, 22, 24, 21, 28, 26, 26, 33, 31, 34, 24, 32, 28, 30, 33, 35, 32, 25, 30, 29, 34,
   31, 35, 27, 33, 34, 35
}; // 0.16 kb
static const cart_index_t iCartYY_ab7_a2[126] = {
   0, 9, 10, 5, 7, 21, 11, 1, 12, 3, 22, 8, 13, 14, 2, 23, 4, 6, 15, 3, 24, 11, 29, 22, 16, 25,
   4, 27, 13, 23, 5, 17, 26, 9, 21, 28, 27, 18, 6, 25, 23, 14, 7, 28, 19, 21, 10, 26, 29, 8, 20, 22,
   24, 12, 9, 15, 33, 17, 28, 30, 10, 33, 16, 26, 19, 31, 17, 11, 34, 15, 30, 29, 34, 12, 18, 24, 32, 20,
   19, 35, 13, 31, 16, 27, 35, 20, 14, 32, 25, 18, 21, 30, 31, 28, 26, 33, 30, 22, 32, 29, 34, 24, 31, 32,
   23, 35, 27, 25, 33, 24, 25, 34, 35, 32, 26, 34, 27, 33, 31, 35, 28, 29, 35, 30, 33, 34
}; // 0.25 kb
static const cart_index_t iCartYY_ab7_a3[150] = {
   0, 17, 19, 9, 10, 5, 26, 7, 28, 21, 15, 1, 20, 3, 24, 11, 12, 29, 8, 22, 16, 18, 2, 25, 4, 27,
   6, 13, 14, 23, 5, 15, 31, 17, 26, 9, 33, 21, 30, 28, 17, 3, 32, 11, 34, 15, 24, 30, 22, 29, 7, 30,
   16, 28, 19, 21, 31, 10, 33, 26, 19, 32, 4, 35, 13, 31, 23, 16, 25, 27, 30, 8, 18, 22, 32, 29, 20, 34,
   12, 24, 31, 20, 6, 32, 23, 35, 14, 27, 18, 25, 9, 11, 35, 15, 33, 17, 34, 28, 29, 30, 10, 34, 13, 33,
   16, 26, 27, 19, 35, 31, 33, 12, 14, 24, 25, 34, 18, 35, 20, 32, 26, 24, 23, 34, 27, 33, 25, 31, 32, 35,
   28, 22, 25, 29, 35, 30, 32, 33, 24, 34, 21, 29, 27, 30, 31, 28, 35, 26, 34, 33
}; // 0.29 kb
static const cart_index_t iCartYY_ab7_a4[150] = {
   0, 15, 16, 5, 17, 7, 19, 30, 31, 9, 10, 33, 26, 28, 21, 17, 1, 18, 15, 3, 30, 32, 8, 20, 11, 34,
   12, 24, 22, 29, 19, 20, 2, 31, 32, 16, 4, 18, 6, 35, 13, 14, 23, 25, 27, 9, 3, 25, 17, 11, 28, 35,
   22, 32, 15, 33, 24, 34, 29, 30, 10, 24, 4, 26, 34, 19, 13, 32, 23, 33, 16, 25, 27, 35, 31, 5, 11, 27,
   9, 15, 21, 31, 29, 35, 17, 26, 34, 33, 30, 28, 26, 12, 6, 33, 24, 31, 23, 20, 14, 34, 27, 18, 25, 32,
   35, 7, 29, 13, 21, 30, 10, 16, 34, 27, 28, 19, 35, 31, 33, 26, 28, 8, 14, 30, 22, 33, 25, 12, 18, 29,
   35, 20, 32, 24, 34, 21, 22, 23, 28, 29, 26, 27, 24, 25, 30, 31, 32, 35, 34, 33
}; // 0.29 kb
static const cart_index_t iCartYY_ab7_a5[126] = {
   0, 11, 13, 15, 16, 5, 27, 7, 29, 9, 10, 17, 34, 19, 35, 21, 30, 31, 33, 26, 28, 9, 1, 14, 3, 25,
   17, 18, 28, 8, 15, 33, 11, 12, 35, 20, 30, 22, 32, 24, 34, 29, 10, 12, 2, 24, 4, 26, 6, 19, 20, 33,
   16, 34, 18, 13, 14, 31, 32, 23, 25, 27, 35, 5, 3, 23, 11, 27, 9, 25, 21, 22, 17, 26, 15, 24, 31, 32,
   28, 29, 35, 34, 33, 30, 7, 22, 4, 29, 13, 21, 23, 10, 24, 28, 19, 30, 32, 16, 25, 26, 34, 27, 35, 31,
   33, 21, 8, 6, 22, 23, 28, 14, 26, 12, 30, 31, 29, 20, 27, 18, 33, 24, 25, 32, 35, 34
}; // 0.25 kb
static const cart_index_t iCartYY_ab7_a6[84] = {
   0, 3, 4, 5, 11, 7, 13, 22, 23, 9, 10, 15, 16, 24, 25, 17, 19, 32, 27, 29, 21, 26, 34, 28, 35, 30,
   31, 33, 5, 1, 6, 9, 3, 21, 23, 8, 14, 17, 26, 11, 27, 12, 18, 15, 31, 20, 25, 22, 28, 33, 24, 30,
   32, 29, 35, 34, 7, 8, 2, 21, 22, 10, 4, 12, 6, 28, 19, 29, 13, 20, 14, 30, 16, 18, 23, 24, 26, 31,
   32, 33, 25, 34, 27, 35
}; // 0.16 kb
static const cart_index_t iCartYY_ab7_a7[36] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35
}; // 0.07 kb
static const cart_index_t iCartYY_ab8_a0[45] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44
}; // 0.09 kb
static const cart_index_t iCartYY_ab8_a1[108] = {
   0, 3, 5, 4, 1, 7, 6, 8, 2, 11, 4, 25, 12, 24, 6, 3, 9, 26, 24, 14, 8, 5, 26, 10, 25, 7,
   13, 9, 15, 29, 10, 27, 17, 16, 11, 31, 28, 13, 19, 18, 32, 12, 30, 20, 14, 21, 16, 35, 22, 33, 18, 15,
   21, 37, 34, 23, 20, 17, 38, 22, 36, 19, 23, 26, 29, 27, 31, 25, 28, 32, 30, 24, 40, 28, 36, 41, 34, 30,
   27, 39, 38, 33, 41, 32, 29, 37, 39, 35, 31, 40, 37, 35, 42, 38, 43, 33, 44, 36, 34, 39, 42, 43, 42, 40,
   44, 43, 44, 41
}; // 0.21 kb
static const cart_index_t iCartYY_ab8_a2[168] = {
   0, 9, 10, 3, 5, 26, 11, 1, 13, 4, 25, 7, 12, 14, 2, 24, 6, 8, 3, 15, 27, 9, 26, 29, 16, 4,
   28, 11, 31, 25, 5, 29, 17, 26, 10, 27, 18, 30, 6, 32, 12, 24, 31, 7, 19, 25, 28, 13, 32, 20, 8, 30,
   24, 14, 9, 21, 39, 15, 29, 37, 10, 39, 22, 27, 17, 38, 21, 11, 40, 16, 35, 31, 22, 41, 12, 33, 18, 32,
   40, 13, 23, 28, 36, 19, 41, 23, 14, 34, 30, 20, 15, 16, 42, 21, 37, 35, 17, 43, 18, 38, 22, 33, 44, 19,
   20, 36, 34, 23, 33, 34, 24, 41, 32, 30, 35, 25, 36, 31, 40, 28, 26, 37, 38, 29, 27, 39, 27, 42, 33, 39,
   38, 43, 42, 28, 34, 40, 44, 36, 29, 35, 43, 37, 39, 42, 43, 36, 30, 44, 41, 34, 37, 31, 44, 35, 42, 40,
   38, 44, 32, 43, 33, 41, 39, 40, 41, 42, 43, 44
}; // 0.33 kb
static const cart_index_t iCartYY_ab8_a3[210] = {
   0, 15, 17, 9, 10, 3, 27, 5, 29, 26, 16, 1, 19, 4, 28, 11, 13, 31, 7, 25, 18, 20, 2, 30, 6, 32,
   8, 12, 14, 24, 21, 4, 36, 11, 40, 16, 28, 35, 25, 31, 22, 34, 6, 41, 12, 33, 24, 18, 30, 32, 3, 21,
   38, 15, 27, 9, 39, 26, 37, 29, 33, 23, 8, 34, 24, 41, 14, 32, 20, 30, 5, 37, 22, 29, 17, 26, 38, 10,
   39, 27, 35, 7, 23, 25, 36, 31, 19, 40, 13, 28, 9, 16, 43, 21, 39, 15, 42, 29, 35, 37, 10, 42, 18, 39,
   22, 27, 33, 17, 43, 38, 15, 11, 44, 16, 42, 21, 40, 37, 31, 35, 42, 13, 20, 28, 34, 40, 23, 44, 19, 36,
   17, 44, 12, 43, 18, 38, 32, 22, 41, 33, 43, 19, 14, 36, 30, 44, 20, 41, 23, 34, 26, 35, 33, 37, 38, 29,
   43, 27, 42, 39, 37, 25, 34, 31, 44, 35, 36, 42, 28, 40, 38, 36, 24, 44, 32, 43, 30, 33, 34, 41, 39, 28,
   30, 40, 41, 42, 34, 43, 36, 44, 27, 40, 32, 42, 33, 39, 41, 38, 44, 43, 29, 31, 41, 35, 43, 37, 44, 39,
   40, 42
}; // 0.41 kb
static const cart_index_t iCartYY_ab8_a4[225] = {
   0, 21, 22, 3, 15, 5, 17, 37, 38, 9, 10, 39, 27, 29, 26, 21, 1, 23, 16, 4, 35, 36, 7, 19, 11, 40,
   13, 28, 25, 31, 22, 23, 2, 33, 34, 18, 6, 20, 8, 41, 12, 14, 24, 30, 32, 3, 16, 33, 9, 21, 26, 38,
   35, 43, 15, 27, 42, 39, 37, 29, 15, 4, 34, 21, 11, 37, 44, 25, 36, 16, 42, 28, 40, 31, 35, 5, 35, 18,
   26, 37, 10, 22, 42, 33, 29, 17, 43, 38, 39, 27, 17, 36, 6, 38, 44, 22, 12, 34, 24, 43, 18, 30, 32, 41,
   33, 37, 7, 20, 35, 25, 42, 34, 13, 23, 31, 44, 19, 36, 28, 40, 38, 19, 8, 43, 36, 33, 24, 23, 14, 44,
   32, 20, 30, 34, 41, 9, 11, 41, 15, 16, 29, 43, 31, 44, 21, 39, 40, 42, 35, 37, 10, 40, 12, 27, 42, 17,
   18, 44, 32, 39, 22, 41, 33, 43, 38, 39, 13, 14, 42, 28, 43, 30, 19, 20, 40, 41, 23, 34, 36, 44, 27, 28,
   24, 39, 40, 38, 32, 36, 30, 42, 33, 34, 41, 44, 43, 29, 25, 30, 37, 31, 39, 41, 28, 34, 35, 43, 36, 44,
   40, 42, 26, 31, 32, 29, 35, 27, 33, 40, 41, 37, 38, 44, 43, 42, 39
}; // 0.44 kb
static const cart_index_t iCartYY_ab8_a5[210] = {
   0, 16, 18, 21, 22, 3, 33, 5, 35, 9, 10, 15, 42, 17, 43, 26, 37, 38, 39, 27, 29, 15, 1, 20, 4, 34,
   21, 23, 37, 7, 16, 42, 11, 13, 44, 19, 35, 25, 36, 28, 40, 31, 17, 19, 2, 36, 6, 38, 8, 22, 23, 43,
   18, 44, 20, 12, 14, 33, 34, 24, 30, 32, 41, 9, 4, 30, 11, 41, 15, 34, 29, 25, 21, 39, 16, 28, 43, 36,
   37, 31, 44, 40, 42, 35, 10, 28, 6, 40, 12, 27, 24, 17, 36, 39, 22, 42, 34, 18, 30, 38, 44, 32, 41, 33,
   43, 3, 11, 32, 16, 33, 9, 41, 26, 31, 15, 27, 21, 40, 38, 44, 29, 35, 43, 42, 39, 37, 27, 13, 8, 28,
   24, 39, 14, 38, 19, 42, 33, 40, 23, 32, 20, 43, 36, 30, 34, 41, 44, 5, 31, 12, 35, 18, 26, 32, 10, 40,
   29, 17, 37, 44, 22, 41, 27, 42, 33, 43, 38, 39, 29, 7, 14, 25, 30, 37, 20, 39, 13, 35, 43, 31, 19, 41,
   23, 42, 28, 34, 36, 44, 40, 26, 25, 24, 31, 32, 29, 30, 27, 28, 37, 38, 35, 36, 33, 34, 39, 40, 41, 44,
   43, 42
}; // 0.41 kb
static const cart_index_t iCartYY_ab8_a6[168] = {
   0, 11, 12, 3, 16, 5, 18, 31, 32, 9, 10, 21, 22, 40, 41, 15, 17, 44, 33, 35, 26, 27, 42, 29, 43, 37,
   38, 39, 9, 1, 14, 15, 4, 29, 30, 7, 20, 21, 39, 11, 41, 13, 23, 16, 43, 19, 34, 25, 37, 42, 28, 35,
   36, 31, 44, 40, 10, 13, 2, 27, 28, 17, 6, 19, 8, 39, 22, 40, 12, 23, 14, 42, 18, 20, 24, 36, 38, 33,
   34, 43, 30, 44, 32, 41, 3, 4, 24, 9, 11, 26, 32, 25, 30, 15, 27, 16, 33, 28, 34, 21, 38, 36, 41, 31,
   29, 39, 40, 37, 44, 35, 43, 42, 5, 25, 6, 26, 31, 10, 12, 28, 24, 29, 17, 35, 18, 36, 30, 37, 22, 34,
   32, 40, 27, 38, 44, 39, 41, 42, 33, 43, 26, 7, 8, 29, 25, 27, 24, 13, 14, 37, 38, 31, 32, 19, 20, 35,
   33, 23, 30, 28, 39, 43, 36, 42, 34, 40, 41, 44
}; // 0.33 kb
static const cart_index_t iCartYY_ab8_a7[108] = {
   0, 4, 6, 11, 12, 3, 24, 5, 25, 9, 10, 16, 28, 18, 30, 21, 22, 15, 34, 17, 36, 26, 31, 32, 40, 41,
   27, 33, 29, 35, 37, 38, 44, 39, 42, 43, 3, 1, 8, 4, 24, 9, 14, 26, 7, 15, 27, 11, 13, 32, 20, 16,
   33, 21, 23, 38, 19, 29, 25, 30, 28, 34, 39, 41, 37, 31, 35, 43, 36, 42, 40, 44, 5, 7, 2, 25, 6, 26,
   8, 10, 13, 29, 17, 31, 19, 12, 14, 35, 18, 37, 20, 22, 23, 27, 28, 24, 36, 30, 38, 32, 39, 40, 42, 33,
   34, 43, 44, 41
}; // 0.21 kb
static const cart_index_t iCartYY_ab8_a8[45] = {
   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
   26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44
}; // 0.09 kb

// What this does: [R] x M -> [r] x S(la) x M
// where R = lab and r == la. In order to do that, the cartesian component
// indices for nCartX(lab-la) set of nCartY(la) monomials need to be input manually (via ii).
// The actual parameters are thus: N == nSets (e.g., nCartX(lab-la)), ii: nCartY(la) x N integer array.
// M: number of sets to transform (M could in principle be done as an outer loop over the function)
void ShTrN_Indirect(double *IR_RP pOut, unsigned so, double const *IR_RP pIn, unsigned si, unsigned la, cart_index_t const *ii, unsigned N, unsigned M)
{
   switch(la) {
      case 0: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+1*n];
            for (unsigned m = 0; m != M; ++ m) {
               pOut[m*so+N*0+n] = pIn[m*si+j0];
            }
         }
         return;
      }
      case 1: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+3*n], j1 = ii[1+3*n], j2 = ii[2+3*n];
            for (unsigned m = 0; m != M; ++ m) {
               pOut[m*so+N*0+n] = pIn[m*si+j0];
               pOut[m*so+N*1+n] = pIn[m*si+j1];
               pOut[m*so+N*2+n] = pIn[m*si+j2];
            }
         }
         return;
      }
      case 2: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+6*n], j1 = ii[1+6*n], j2 = ii[2+6*n], j3 = ii[3+6*n], j4 = ii[4+6*n], j5 = ii[5+6*n];
            for (unsigned m = 0; m != M; ++ m) {
               double z0 = pIn[m*si+j0], z1 = pIn[m*si+j1];
               pOut[m*so+N*0+n] = pIn[m*si+j2] - sd0*z0 - sd0*z1;
               pOut[m*so+N*1+n] = pIn[m*si+j3]*sd1;
               pOut[m*so+N*2+n] = pIn[m*si+j4]*sd1;
               pOut[m*so+N*3+n] = sd2*z0 - sd2*z1;
               pOut[m*so+N*4+n] = pIn[m*si+j5]*sd1;
            }
         }
         return;
      }
      case 3: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+10*n], j1 = ii[1+10*n], j2 = ii[2+10*n], j3 = ii[3+10*n], j4 = ii[4+10*n], j5 = ii[5+10*n], j6 = ii[6+10*n], j7 = ii[7+10*n], j8 = ii[8+10*n], j9 = ii[9+10*n];
            for (unsigned m = 0; m != M; ++ m) {
               double z0 = pIn[m*si+j0], z1 = pIn[m*si+j1], z3 = pIn[m*si+j3], z5 = pIn[m*si+j5], z7 = pIn[m*si+j7], z8 = pIn[m*si+j8];
               pOut[m*so+N*0+n] = pIn[m*si+j4]*sd4 - sd3*z0 - sd3*z3;
               pOut[m*so+N*1+n] = pIn[m*si+j6]*sd4 - sd3*z1 - sd3*z5;
               pOut[m*so+N*2+n] = pIn[m*si+j2] - sd5*z7 - sd5*z8;
               pOut[m*so+N*3+n] = sd6*z0 - sd7*z3;
               pOut[m*so+N*4+n] = pIn[m*si+j9]*sd8;
               pOut[m*so+N*5+n] = -sd6*z1 + sd7*z5;
               pOut[m*so+N*6+n] = sd9*z7 - sd9*z8;
            }
         }
         return;
      }
      case 4: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+15*n], j1 = ii[1+15*n], j2 = ii[2+15*n], j3 = ii[3+15*n], j4 = ii[4+15*n], j5 = ii[5+15*n], j6 = ii[6+15*n], j7 = ii[7+15*n], j8 = ii[8+15*n], j9 = ii[9+15*n], ja = ii[10+15*n], jb = ii[11+15*n], jc = ii[12+15*n], jd = ii[13+15*n], je = ii[14+15*n];
            for (unsigned m = 0; m != M; ++ m) {
               double z0 = pIn[m*si+j0], z1 = pIn[m*si+j1], z3 = pIn[m*si+j3], z4 = pIn[m*si+j4], z5 = pIn[m*si+j5], z7 = pIn[m*si+j7], z9 = pIn[m*si+j9], za = pIn[m*si+ja], zb = pIn[m*si+jb], zd = pIn[m*si+jd], ze = pIn[m*si+je];
               pOut[m*so+N*0+n] = pIn[m*si+j2] + sda*z0 + sda*z1 + sdb*z9 - sdc*za - sdc*zb;
               pOut[m*so+N*1+n] = pIn[m*si+jc]*sde - sdd*z3 - sdd*z4;
               pOut[m*so+N*2+n] = pIn[m*si+j6]*sdf - sd7*z5 - sd7*zd;
               pOut[m*so+N*3+n] = sd10*z0 + sd10*z1 - sd11*z9;
               pOut[m*so+N*4+n] = pIn[m*si+j8]*sdf - sd7*z7 - sd7*ze;
               pOut[m*so+N*5+n] = -sd12*z0 + sd12*z1 + sd13*za - sd13*zb;
               pOut[m*so+N*6+n] = sd14*z3 - sd14*z4;
               pOut[m*so+N*7+n] = sd15*z5 - sd16*zd;
               pOut[m*so+N*8+n] = -sd15*z7 + sd16*ze;
            }
         }
         return;
      }
      case 5: {
         for (unsigned n = 0; n != N; ++ n) {
            unsigned j0 = ii[0+21*n], j1 = ii[1+21*n], j2 = ii[2+21*n], j3 = ii[3+21*n], j4 = ii[4+21*n], j5 = ii[5+21*n], j6 = ii[6+21*n], j7 = ii[7+21*n], j8 = ii[8+21*n], j9 = ii[9+21*n], ja = ii[10+21*n], jb = ii[11+21*n], jc = ii[12+21*n], jd = ii[13+21*n], je = ii[14+21*n], jf = ii[15+21*n], j10 = ii[16+21*n], j11 = ii[17+21*n], j12 = ii[18+21*n], j13 = ii[19+21*n], j14 = ii[20+21*n];
            for (unsigned m = 0; m != M; ++ m) {
               double z0 = pIn[m*si+j0], z1 = pIn[m*si+j1], z3 = pIn[m*si+j3], z5 = pIn[m*si+j5], z7 = pIn[m*si+j7], z8 = pIn[m*si+j8], z9 = pIn[m*si+j9], za = pIn[m*si+ja], zb = pIn[m*si+jb], zc = pIn[m*si+jc], zd = pIn[m*si+jd], ze = pIn[m*si+je], zf = pIn[m*si+jf], z10 = pIn[m*si+j10], z12 = pIn[m*si+j12], z13 = pIn[m*si+j13], z14 = pIn[m*si+j14];
               pOut[m*so+N*0+n] = pIn[m*si+j4]*sd8 + sd17*z0 + sd17*z3 + sd18*z9 - sd19*z12 - sd19*za;
               pOut[m*so+N*1+n] = pIn[m*si+j6]*sd8 + sd17*z1 + sd17*z5 + sd18*zb - sd19*z13 - sd19*zc;
               pOut[m*so+N*2+n] = -sd1a*z7 + sd1a*z8 + sd1b*zd - sd1b*ze;
               pOut[m*so+N*3+n] = -sd1c*z0 + sd1d*z9 + sd1e*za + sd1f*z3 - sd20*z12;
               pOut[m*so+N*4+n] = -sd21*z10 + sd21*zf;
               pOut[m*so+N*5+n] = sd1c*z1 - sd1d*zb - sd1e*zc - sd1f*z5 + sd20*z13;
               pOut[m*so+N*6+n] = sd22*z7 + sd22*z8 - sd23*z14;
               pOut[m*so+N*7+n] = sd24*z5 - sd25*zb + sd26*z1;
               pOut[m*so+N*8+n] = pIn[m*si+j2] + sd27*z7 + sd27*z8 + sd28*z14 - sd29*zd - sd29*ze;
               pOut[m*so+N*9+n] = sd24*z3 - sd25*z9 + sd26*z0;
               pOut[m*so+N*10+n] = pIn[m*si+j11]*sd2a - sd1b*z10 - sd1b*zf;
            }
         }
         return;
      }
   }
   assert(0);
}


// Factorize a matrix (nCartX(lab)-nCartX(a-1)) x M into nCartX(lab-la) x (2*la+1) x M. Asserts la >= lb (= lab-la)
void ShTrA_XY(double *IR_RP pOut, double const *IR_RP pIn, unsigned la, unsigned lab, unsigned M)
{
   assert(lab >= la && lab-la <= la);
   #define C_AB_A(lab,la) ((lab*(lab-1))/2 + la)
   switch(C_AB_A(lab,la)) {
      case C_AB_A( 0, 0): return ShTrN_Indirect(pOut,  1, pIn,   1, 0, &iCartXY_ab0_a0[0],  1, M); // case: 0
      case C_AB_A( 1, 1): return ShTrN_Indirect(pOut,  3, pIn,   3, 1, &iCartXY_ab1_a1[0],  1, M); // case: 1
      case C_AB_A( 2, 1): return ShTrN_Indirect(pOut, 12, pIn,   9, 1, &iCartXY_ab2_a1[0],  4, M); // case: 2
      case C_AB_A( 2, 2): return ShTrN_Indirect(pOut,  5, pIn,   6, 2, &iCartXY_ab2_a2[0],  1, M); // case: 3
      case C_AB_A( 3, 2): return ShTrN_Indirect(pOut, 20, pIn,  16, 2, &iCartXY_ab3_a2[0],  4, M); // case: 5
      case C_AB_A( 3, 3): return ShTrN_Indirect(pOut,  7, pIn,  10, 3, &iCartXY_ab3_a3[0],  1, M); // case: 6
      case C_AB_A( 4, 2): return ShTrN_Indirect(pOut, 50, pIn,  31, 2, &iCartXY_ab4_a2[0], 10, M); // case: 8
      case C_AB_A( 4, 3): return ShTrN_Indirect(pOut, 28, pIn,  25, 3, &iCartXY_ab4_a3[0],  4, M); // case: 9
      case C_AB_A( 4, 4): return ShTrN_Indirect(pOut,  9, pIn,  15, 4, &iCartXY_ab4_a4[0],  1, M); // case: 10
      case C_AB_A( 5, 3): return ShTrN_Indirect(pOut, 70, pIn,  46, 3, &iCartXY_ab5_a3[0], 10, M); // case: 13
      case C_AB_A( 5, 4): return ShTrN_Indirect(pOut, 36, pIn,  36, 4, &iCartXY_ab5_a4[0],  4, M); // case: 14
      case C_AB_A( 5, 5): return ShTrN_Indirect(pOut, 11, pIn,  21, 5, &iCartXY_ab5_a5[0],  1, M); // case: 15
      case C_AB_A( 6, 3): return ShTrN_Indirect(pOut,140, pIn,  74, 3, &iCartXY_ab6_a3[0], 20, M); // case: 18
      case C_AB_A( 6, 4): return ShTrN_Indirect(pOut, 90, pIn,  64, 4, &iCartXY_ab6_a4[0], 10, M); // case: 19
      case C_AB_A( 6, 5): return ShTrN_Indirect(pOut, 44, pIn,  49, 5, &iCartXY_ab6_a5[0],  4, M); // case: 20
      case C_AB_A( 6, 6): return ShTrN_Indirect(pOut, 13, pIn,  28, 6, &iCartXY_ab6_a6[0],  1, M); // case: 21
      case C_AB_A( 7, 4): return ShTrN_Indirect(pOut,180, pIn, 100, 4, &iCartXY_ab7_a4[0], 20, M); // case: 25
      case C_AB_A( 7, 5): return ShTrN_Indirect(pOut,110, pIn,  85, 5, &iCartXY_ab7_a5[0], 10, M); // case: 26
      case C_AB_A( 7, 6): return ShTrN_Indirect(pOut, 52, pIn,  64, 6, &iCartXY_ab7_a6[0],  4, M); // case: 27
      case C_AB_A( 7, 7): return ShTrN_Indirect(pOut, 15, pIn,  36, 7, &iCartXY_ab7_a7[0],  1, M); // case: 28
      case C_AB_A( 8, 4): return ShTrN_Indirect(pOut,315, pIn, 145, 4, &iCartXY_ab8_a4[0], 35, M); // case: 32
      case C_AB_A( 8, 5): return ShTrN_Indirect(pOut,220, pIn, 130, 5, &iCartXY_ab8_a5[0], 20, M); // case: 33
      case C_AB_A( 8, 6): return ShTrN_Indirect(pOut,130, pIn, 109, 6, &iCartXY_ab8_a6[0], 10, M); // case: 34
      case C_AB_A( 8, 7): return ShTrN_Indirect(pOut, 60, pIn,  81, 7, &iCartXY_ab8_a7[0],  4, M); // case: 35
      case C_AB_A( 8, 8): return ShTrN_Indirect(pOut, 17, pIn,  45, 8, &iCartXY_ab8_a8[0],  1, M); // case: 36
   }
   #undef C_AB_A
   assert(0);
}


// Factorize a matrix nCartY(lab) x M into nCartY(lab-la) x (2*la+1) x M. Asserts la >= lb (= lab-la)
void ShTrA_YY(double *IR_RP pOut, double const *IR_RP pIn, unsigned la, unsigned lab, unsigned M)
{
   assert(lab >= la && lab-la <= la);
   #define C_AB_A(lab,la) ((lab*(lab-1))/2 + la)
   switch(C_AB_A(lab,la)) {
      case C_AB_A( 0, 0): return ShTrN_Indirect(pOut,  1, pIn,   1, 0, &iCartYY_ab0_a0[0],  1, M); // case: 0
      case C_AB_A( 1, 1): return ShTrN_Indirect(pOut,  3, pIn,   3, 1, &iCartYY_ab1_a1[0],  1, M); // case: 1
      case C_AB_A( 2, 1): return ShTrN_Indirect(pOut,  9, pIn,   6, 1, &iCartYY_ab2_a1[0],  3, M); // case: 2
      case C_AB_A( 2, 2): return ShTrN_Indirect(pOut,  5, pIn,   6, 2, &iCartYY_ab2_a2[0],  1, M); // case: 3
      case C_AB_A( 3, 2): return ShTrN_Indirect(pOut, 15, pIn,  10, 2, &iCartYY_ab3_a2[0],  3, M); // case: 5
      case C_AB_A( 3, 3): return ShTrN_Indirect(pOut,  7, pIn,  10, 3, &iCartYY_ab3_a3[0],  1, M); // case: 6
      case C_AB_A( 4, 2): return ShTrN_Indirect(pOut, 30, pIn,  15, 2, &iCartYY_ab4_a2[0],  6, M); // case: 8
      case C_AB_A( 4, 3): return ShTrN_Indirect(pOut, 21, pIn,  15, 3, &iCartYY_ab4_a3[0],  3, M); // case: 9
      case C_AB_A( 4, 4): return ShTrN_Indirect(pOut,  9, pIn,  15, 4, &iCartYY_ab4_a4[0],  1, M); // case: 10
      case C_AB_A( 5, 3): return ShTrN_Indirect(pOut, 42, pIn,  21, 3, &iCartYY_ab5_a3[0],  6, M); // case: 13
      case C_AB_A( 5, 4): return ShTrN_Indirect(pOut, 27, pIn,  21, 4, &iCartYY_ab5_a4[0],  3, M); // case: 14
      case C_AB_A( 5, 5): return ShTrN_Indirect(pOut, 11, pIn,  21, 5, &iCartYY_ab5_a5[0],  1, M); // case: 15
      case C_AB_A( 6, 3): return ShTrN_Indirect(pOut, 70, pIn,  28, 3, &iCartYY_ab6_a3[0], 10, M); // case: 18
      case C_AB_A( 6, 4): return ShTrN_Indirect(pOut, 54, pIn,  28, 4, &iCartYY_ab6_a4[0],  6, M); // case: 19
      case C_AB_A( 6, 5): return ShTrN_Indirect(pOut, 33, pIn,  28, 5, &iCartYY_ab6_a5[0],  3, M); // case: 20
      case C_AB_A( 6, 6): return ShTrN_Indirect(pOut, 13, pIn,  28, 6, &iCartYY_ab6_a6[0],  1, M); // case: 21
      case C_AB_A( 7, 4): return ShTrN_Indirect(pOut, 90, pIn,  36, 4, &iCartYY_ab7_a4[0], 10, M); // case: 25
      case C_AB_A( 7, 5): return ShTrN_Indirect(pOut, 66, pIn,  36, 5, &iCartYY_ab7_a5[0],  6, M); // case: 26
      case C_AB_A( 7, 6): return ShTrN_Indirect(pOut, 39, pIn,  36, 6, &iCartYY_ab7_a6[0],  3, M); // case: 27
      case C_AB_A( 7, 7): return ShTrN_Indirect(pOut, 15, pIn,  36, 7, &iCartYY_ab7_a7[0],  1, M); // case: 28
      case C_AB_A( 8, 4): return ShTrN_Indirect(pOut,135, pIn,  45, 4, &iCartYY_ab8_a4[0], 15, M); // case: 32
      case C_AB_A( 8, 5): return ShTrN_Indirect(pOut,110, pIn,  45, 5, &iCartYY_ab8_a5[0], 10, M); // case: 33
      case C_AB_A( 8, 6): return ShTrN_Indirect(pOut, 78, pIn,  45, 6, &iCartYY_ab8_a6[0],  6, M); // case: 34
      case C_AB_A( 8, 7): return ShTrN_Indirect(pOut, 45, pIn,  45, 7, &iCartYY_ab8_a7[0],  3, M); // case: 35
      case C_AB_A( 8, 8): return ShTrN_Indirect(pOut, 17, pIn,  45, 8, &iCartYY_ab8_a8[0],  1, M); // case: 36
   }
   #undef C_AB_A
   assert(0);
}

// transform cartesians CartX[(x-A)^lb], centered at A, to solid harmoincs Slm(x-B), centered at B.
// Input is nCartX(lb) x nCount, output is (2l+1); sa indexes Count, sb indexes Slm(x-B).
void OsrrC(double *IR_RP pOut, unsigned sa, unsigned sb, double const *IR_RP pIn, double AmBx, double AmBy, double AmBz, unsigned lb, unsigned nCount)
{
   switch(lb) {
      case 0: {
         for (unsigned ia = 0; ia < nCount; ++ ia) {
            double const *IR_RP pAx = pIn + 1 * ia;
            pOut[sa*ia + sb*0] = pAx[0];
         }
         return;
      }
      case 1: {
         for (unsigned ia = 0; ia < nCount; ++ ia) {
            double const *IR_RP pAx = pIn + 4 * ia;
            double b_100 = AmBx*pAx[0] + pAx[1];
            double b_010 = AmBy*pAx[0] + pAx[2];
            double b_001 = AmBz*pAx[0] + pAx[3];
            pOut[sa*ia + sb*0] = b_100;
            pOut[sa*ia + sb*1] = b_010;
            pOut[sa*ia + sb*2] = b_001;
         }
         return;
      }
      case 2: {
         for (unsigned ia = 0; ia < nCount; ++ ia) {
            double const *IR_RP pAx = pIn + 10 * ia;
            double r_000_100 = AmBx*pAx[0] + pAx[1];
            double r_100_100 = AmBx*pAx[1] + pAx[4];
            double r_010_100 = AmBx*pAx[2] + pAx[7];
            double r_001_100 = AmBx*pAx[3] + pAx[8];
            double r_000_010 = AmBy*pAx[0] + pAx[2];
            double r_010_010 = AmBy*pAx[2] + pAx[5];
            double r_001_010 = AmBy*pAx[3] + pAx[9];
            double r_000_001 = AmBz*pAx[0] + pAx[3];
            double r_001_001 = AmBz*pAx[3] + pAx[6];
            double b_200 = AmBx*r_000_100 + r_100_100;
            double b_020 = AmBy*r_000_010 + r_010_010;
            double b_002 = AmBz*r_000_001 + r_001_001;
            double b_110 = AmBy*r_000_100 + r_010_100;
            double b_101 = AmBz*r_000_100 + r_001_100;
            double b_011 = AmBz*r_000_010 + r_001_010;
            pOut[sa*ia + sb*0] = b_002 - b_020*sd0 - b_200*sd0;
            pOut[sa*ia + sb*1] = b_110*sd1;
            pOut[sa*ia + sb*2] = b_101*sd1;
            pOut[sa*ia + sb*3] = -b_020*sd2 + b_200*sd2;
            pOut[sa*ia + sb*4] = b_011*sd1;
         }
         return;
      }
      case 3: {
         for (unsigned ia = 0; ia < nCount; ++ ia) {
            double const *IR_RP pAx = pIn + 20 * ia;
            double r_000_100 = AmBx*pAx[0] + pAx[1];
            double r_100_100 = AmBx*pAx[1] + pAx[4];
            double r_010_100 = AmBx*pAx[2] + pAx[7];
            double r_001_100 = AmBx*pAx[3] + pAx[8];
            double r_200_100 = AmBx*pAx[4] + pAx[10];
            double r_110_100 = AmBx*pAx[7] + pAx[15];
            double r_101_100 = AmBx*pAx[8] + pAx[17];
            double r_011_100 = AmBx*pAx[9] + pAx[19];
            double r_000_010 = AmBy*pAx[0] + pAx[2];
            double r_100_010 = AmBy*pAx[1] + pAx[7];
            double r_010_010 = AmBy*pAx[2] + pAx[5];
            double r_001_010 = AmBy*pAx[3] + pAx[9];
            double r_020_010 = AmBy*pAx[5] + pAx[11];
            double r_110_010 = AmBy*pAx[7] + pAx[13];
            double r_011_010 = AmBy*pAx[9] + pAx[18];
            double r_000_001 = AmBz*pAx[0] + pAx[3];
            double r_100_001 = AmBz*pAx[1] + pAx[8];
            double r_010_001 = AmBz*pAx[2] + pAx[9];
            double r_001_001 = AmBz*pAx[3] + pAx[6];
            double r_002_001 = AmBz*pAx[6] + pAx[12];
            double r_101_001 = AmBz*pAx[8] + pAx[14];
            double r_011_001 = AmBz*pAx[9] + pAx[16];
            double r_000_200 = AmBx*r_000_100 + r_100_100;
            double r_100_200 = AmBx*r_100_100 + r_200_100;
            double r_010_200 = AmBx*r_010_100 + r_110_100;
            double r_001_200 = AmBx*r_001_100 + r_101_100;
            double r_000_020 = AmBy*r_000_010 + r_010_010;
            double r_100_020 = AmBy*r_100_010 + r_110_010;
            double r_010_020 = AmBy*r_010_010 + r_020_010;
            double r_001_020 = AmBy*r_001_010 + r_011_010;
            double r_000_002 = AmBz*r_000_001 + r_001_001;
            double r_100_002 = AmBz*r_100_001 + r_101_001;
            double r_010_002 = AmBz*r_010_001 + r_011_001;
            double r_001_002 = AmBz*r_001_001 + r_002_001;
            double r_000_110 = AmBy*r_000_100 + r_010_100;
            double r_001_110 = AmBy*r_001_100 + r_011_100;
            double b_300 = AmBx*r_000_200 + r_100_200;
            double b_030 = AmBy*r_000_020 + r_010_020;
            double b_003 = AmBz*r_000_002 + r_001_002;
            double b_120 = AmBx*r_000_020 + r_100_020;
            double b_102 = AmBx*r_000_002 + r_100_002;
            double b_210 = AmBy*r_000_200 + r_010_200;
            double b_012 = AmBy*r_000_002 + r_010_002;
            double b_201 = AmBz*r_000_200 + r_001_200;
            double b_021 = AmBz*r_000_020 + r_001_020;
            double b_111 = AmBz*r_000_110 + r_001_110;
            pOut[sa*ia + sb*0] = b_102*sd4 - b_120*sd3 - b_300*sd3;
            pOut[sa*ia + sb*1] = b_012*sd4 - b_030*sd3 - b_210*sd3;
            pOut[sa*ia + sb*2] = b_003 - b_021*sd5 - b_201*sd5;
            pOut[sa*ia + sb*3] = -b_120*sd7 + b_300*sd6;
            pOut[sa*ia + sb*4] = b_111*sd8;
            pOut[sa*ia + sb*5] = -b_030*sd6 + b_210*sd7;
            pOut[sa*ia + sb*6] = -b_021*sd9 + b_201*sd9;
         }
         return;
      }
      case 4: {
         for (unsigned ia = 0; ia < nCount; ++ ia) {
            double const *IR_RP pAx = pIn + 35 * ia;
            double r_000_100 = AmBx*pAx[0] + pAx[1];
            double r_100_100 = AmBx*pAx[1] + pAx[4];
            double r_010_100 = AmBx*pAx[2] + pAx[7];
            double r_001_100 = AmBx*pAx[3] + pAx[8];
            double r_200_100 = AmBx*pAx[4] + pAx[10];
            double r_110_100 = AmBx*pAx[7] + pAx[15];
            double r_101_100 = AmBx*pAx[8] + pAx[17];
            double r_011_100 = AmBx*pAx[9] + pAx[19];
            double r_300_100 = AmBx*pAx[10] + pAx[20];
            double r_210_100 = AmBx*pAx[15] + pAx[23];
            double r_201_100 = AmBx*pAx[17] + pAx[25];
            double r_111_100 = AmBx*pAx[19] + pAx[34];
            double r_000_010 = AmBy*pAx[0] + pAx[2];
            double r_100_010 = AmBy*pAx[1] + pAx[7];
            double r_010_010 = AmBy*pAx[2] + pAx[5];
            double r_001_010 = AmBy*pAx[3] + pAx[9];
            double r_200_010 = AmBy*pAx[4] + pAx[15];
            double r_020_010 = AmBy*pAx[5] + pAx[11];
            double r_110_010 = AmBy*pAx[7] + pAx[13];
            double r_101_010 = AmBy*pAx[8] + pAx[19];
            double r_011_010 = AmBy*pAx[9] + pAx[18];
            double r_030_010 = AmBy*pAx[11] + pAx[21];
            double r_120_010 = AmBy*pAx[13] + pAx[24];
            double r_210_010 = AmBy*pAx[15] + pAx[29];
            double r_021_010 = AmBy*pAx[18] + pAx[27];
            double r_111_010 = AmBy*pAx[19] + pAx[33];
            double r_000_001 = AmBz*pAx[0] + pAx[3];
            double r_100_001 = AmBz*pAx[1] + pAx[8];
            double r_010_001 = AmBz*pAx[2] + pAx[9];
            double r_001_001 = AmBz*pAx[3] + pAx[6];
            double r_200_001 = AmBz*pAx[4] + pAx[17];
            double r_020_001 = AmBz*pAx[5] + pAx[18];
            double r_002_001 = AmBz*pAx[6] + pAx[12];
            double r_110_001 = AmBz*pAx[7] + pAx[19];
            double r_101_001 = AmBz*pAx[8] + pAx[14];
            double r_011_001 = AmBz*pAx[9] + pAx[16];
            double r_003_001 = AmBz*pAx[12] + pAx[22];
            double r_102_001 = AmBz*pAx[14] + pAx[26];
            double r_012_001 = AmBz*pAx[16] + pAx[28];
            double r_201_001 = AmBz*pAx[17] + pAx[30];
            double r_021_001 = AmBz*pAx[18] + pAx[31];
            double r_111_001 = AmBz*pAx[19] + pAx[32];
            double r_000_200 = AmBx*r_000_100 + r_100_100;
            double r_100_200 = AmBx*r_100_100 + r_200_100;
            double r_010_200 = AmBx*r_010_100 + r_110_100;
            double r_001_200 = AmBx*r_001_100 + r_101_100;
            double r_200_200 = AmBx*r_200_100 + r_300_100;
            double r_110_200 = AmBx*r_110_100 + r_210_100;
            double r_101_200 = AmBx*r_101_100 + r_201_100;
            double r_011_200 = AmBx*r_011_100 + r_111_100;
            double r_000_020 = AmBy*r_000_010 + r_010_010;
            double r_100_020 = AmBy*r_100_010 + r_110_010;
            double r_010_020 = AmBy*r_010_010 + r_020_010;
            double r_001_020 = AmBy*r_001_010 + r_011_010;
            double r_200_020 = AmBy*r_200_010 + r_210_010;
            double r_020_020 = AmBy*r_020_010 + r_030_010;
            double r_110_020 = AmBy*r_110_010 + r_120_010;
            double r_101_020 = AmBy*r_101_010 + r_111_010;
            double r_011_020 = AmBy*r_011_010 + r_021_010;
            double r_000_002 = AmBz*r_000_001 + r_001_001;
            double r_100_002 = AmBz*r_100_001 + r_101_001;
            double r_010_002 = AmBz*r_010_001 + r_011_001;
            double r_001_002 = AmBz*r_001_001 + r_002_001;
            double r_200_002 = AmBz*r_200_001 + r_201_001;
            double r_020_002 = AmBz*r_020_001 + r_021_001;
            double r_002_002 = AmBz*r_002_001 + r_003_001;
            double r_110_002 = AmBz*r_110_001 + r_111_001;
            double r_101_002 = AmBz*r_101_001 + r_102_001;
            double r_011_002 = AmBz*r_011_001 + r_012_001;
            double r_000_300 = AmBx*r_000_200 + r_100_200;
            double r_100_300 = AmBx*r_100_200 + r_200_200;
            double r_010_300 = AmBx*r_010_200 + r_110_200;
            double r_001_300 = AmBx*r_001_200 + r_101_200;
            double r_000_030 = AmBy*r_000_020 + r_010_020;
            double r_100_030 = AmBy*r_100_020 + r_110_020;
            double r_010_030 = AmBy*r_010_020 + r_020_020;
            double r_001_030 = AmBy*r_001_020 + r_011_020;
            double r_000_003 = AmBz*r_000_002 + r_001_002;
            double r_100_003 = AmBz*r_100_002 + r_101_002;
            double r_010_003 = AmBz*r_010_002 + r_011_002;
            double r_001_003 = AmBz*r_001_002 + r_002_002;
            double r_000_120 = AmBx*r_000_020 + r_100_020;
            double r_100_120 = AmBx*r_100_020 + r_200_020;
            double r_001_120 = AmBx*r_001_020 + r_101_020;
            double r_000_102 = AmBx*r_000_002 + r_100_002;
            double r_100_102 = AmBx*r_100_002 + r_200_002;
            double r_010_102 = AmBx*r_010_002 + r_110_002;
            double r_000_210 = AmBy*r_000_200 + r_010_200;
            double r_001_210 = AmBy*r_001_200 + r_011_200;
            double r_000_012 = AmBy*r_000_002 + r_010_002;
            double r_010_012 = AmBy*r_010_002 + r_020_002;
            double b_400 = AmBx*r_000_300 + r_100_300;
            double b_040 = AmBy*r_000_030 + r_010_030;
            double b_004 = AmBz*r_000_003 + r_001_003;
            double b_310 = AmBy*r_000_300 + r_010_300;
            double b_130 = AmBx*r_000_030 + r_100_030;
            double b_301 = AmBz*r_000_300 + r_001_300;
            double b_103 = AmBx*r_000_003 + r_100_003;
            double b_031 = AmBz*r_000_030 + r_001_030;
            double b_013 = AmBy*r_000_003 + r_010_003;
            double b_220 = AmBx*r_000_120 + r_100_120;
            double b_202 = AmBx*r_000_102 + r_100_102;
            double b_022 = AmBy*r_000_012 + r_010_012;
            double b_112 = AmBy*r_000_102 + r_010_102;
            double b_121 = AmBz*r_000_120 + r_001_120;
            double b_211 = AmBz*r_000_210 + r_001_210;
            pOut[sa*ia + sb*0] = b_004 - b_022*sdc + b_040*sda - b_202*sdc + b_220*sdb + b_400*sda;
            pOut[sa*ia + sb*1] = b_112*sde - b_130*sdd - b_310*sdd;
            pOut[sa*ia + sb*2] = b_103*sdf - b_121*sd7 - b_301*sd7;
            pOut[sa*ia + sb*3] = b_040*sd10 - b_220*sd11 + b_400*sd10;
            pOut[sa*ia + sb*4] = b_013*sdf - b_031*sd7 - b_211*sd7;
            pOut[sa*ia + sb*5] = -b_022*sd13 + b_040*sd12 + b_202*sd13 - b_400*sd12;
            pOut[sa*ia + sb*6] = -b_130*sd14 + b_310*sd14;
            pOut[sa*ia + sb*7] = -b_121*sd16 + b_301*sd15;
            pOut[sa*ia + sb*8] = -b_031*sd15 + b_211*sd16;
         }
         return;
      }
   }
   assert(0);
}

// calculate [CartY(lab)]^0 from [0]^m. Effectively, this function calculates"
//    D^r f(t)
// from f^[m](t) = (D/Dt)^m f(t) and R, where t = R^2, D^r means \prod_i (D/D{R_i})^{r_x},
// and f is some arbitrary scalar function (of which you supply the m'th derivatives with
// respect to t as [0]^m).
void ShellMdrr(double *IR_RP pOut, double const *IR_RP pIn, double Rx, double Ry, double Rz, unsigned lab)
{
   switch(lab) {
      case 0: {
         pOut[0] = pIn[0];
         return;
      }
      case 1: {
         pOut[0] = Rx*pIn[1];
         pOut[1] = Ry*pIn[1];
         pOut[2] = Rz*pIn[1];
         return;
      }
      case 2: {
         double r_100_1 = Rx*pIn[2];
         double r_010_1 = Ry*pIn[2];
         double r_001_1 = Rz*pIn[2];
         pOut[0] = Rx*r_100_1 - pIn[1];
         pOut[1] = Ry*r_010_1 - pIn[1];
         pOut[2] = Rz*r_001_1 - pIn[1];
         pOut[3] = Ry*r_100_1;
         pOut[4] = Rz*r_100_1;
         pOut[5] = Rz*r_010_1;
         return;
      }
      case 3: {
         double r_100_2 = Rx*pIn[3];
         double r_010_2 = Ry*pIn[3];
         double r_001_2 = Rz*pIn[3];
         double r_100_1 = Rx*pIn[2];
         double r_010_1 = Ry*pIn[2];
         double r_001_1 = Rz*pIn[2];
         double r_200_1 = Rx*r_100_2 - pIn[2];
         double r_020_1 = Ry*r_010_2 - pIn[2];
         double r_002_1 = Rz*r_001_2 - pIn[2];
         double r_110_1 = Ry*r_100_2;
         pOut[0] = Rx*r_200_1 - 2*r_100_1;
         pOut[1] = Ry*r_020_1 - 2*r_010_1;
         pOut[2] = Rz*r_002_1 - 2*r_001_1;
         pOut[3] = Rx*r_020_1;
         pOut[4] = Rx*r_002_1;
         pOut[5] = Ry*r_200_1;
         pOut[6] = Ry*r_002_1;
         pOut[7] = Rz*r_200_1;
         pOut[8] = Rz*r_020_1;
         pOut[9] = Rz*r_110_1;
         return;
      }
      case 4: {
         double r_100_3 = Rx*pIn[4];
         double r_010_3 = Ry*pIn[4];
         double r_001_3 = Rz*pIn[4];
         double r_100_2 = Rx*pIn[3];
         double r_010_2 = Ry*pIn[3];
         double r_001_2 = Rz*pIn[3];
         double r_200_2 = Rx*r_100_3 - pIn[3];
         double r_020_2 = Ry*r_010_3 - pIn[3];
         double r_002_2 = Rz*r_001_3 - pIn[3];
         double r_200_1 = Rx*r_100_2 - pIn[2];
         double r_020_1 = Ry*r_010_2 - pIn[2];
         double r_002_1 = Rz*r_001_2 - pIn[2];
         double r_300_1 = Rx*r_200_2 - 2*r_100_2;
         double r_030_1 = Ry*r_020_2 - 2*r_010_2;
         double r_003_1 = Rz*r_002_2 - 2*r_001_2;
         double r_120_1 = Rx*r_020_2;
         double r_102_1 = Rx*r_002_2;
         double r_210_1 = Ry*r_200_2;
         double r_012_1 = Ry*r_002_2;
         pOut[0] = Rx*r_300_1 - 3*r_200_1;
         pOut[1] = Ry*r_030_1 - 3*r_020_1;
         pOut[2] = Rz*r_003_1 - 3*r_002_1;
         pOut[3] = Ry*r_300_1;
         pOut[4] = Rx*r_030_1;
         pOut[5] = Rz*r_300_1;
         pOut[6] = Rx*r_003_1;
         pOut[7] = Rz*r_030_1;
         pOut[8] = Ry*r_003_1;
         pOut[9] = Rx*r_120_1 - r_020_1;
         pOut[10] = Rx*r_102_1 - r_002_1;
         pOut[11] = Ry*r_012_1 - r_002_1;
         pOut[12] = Ry*r_102_1;
         pOut[13] = Rz*r_120_1;
         pOut[14] = Rz*r_210_1;
         return;
      }
      case 5: {
         double r_100_4 = Rx*pIn[5];
         double r_010_4 = Ry*pIn[5];
         double r_001_4 = Rz*pIn[5];
         double r_100_3 = Rx*pIn[4];
         double r_010_3 = Ry*pIn[4];
         double r_001_3 = Rz*pIn[4];
         double r_200_3 = Rx*r_100_4 - pIn[4];
         double r_020_3 = Ry*r_010_4 - pIn[4];
         double r_002_3 = Rz*r_001_4 - pIn[4];
         double r_100_2 = Rx*pIn[3];
         double r_010_2 = Ry*pIn[3];
         double r_001_2 = Rz*pIn[3];
         double r_200_2 = Rx*r_100_3 - pIn[3];
         double r_020_2 = Ry*r_010_3 - pIn[3];
         double r_002_2 = Rz*r_001_3 - pIn[3];
         double r_300_2 = Rx*r_200_3 - 2*r_100_3;
         double r_030_2 = Ry*r_020_3 - 2*r_010_3;
         double r_003_2 = Rz*r_002_3 - 2*r_001_3;
         double r_120_2 = Rx*r_020_3;
         double r_102_2 = Rx*r_002_3;
         double r_012_2 = Ry*r_002_3;
         double r_300_1 = Rx*r_200_2 - 2*r_100_2;
         double r_030_1 = Ry*r_020_2 - 2*r_010_2;
         double r_003_1 = Rz*r_002_2 - 2*r_001_2;
         double r_400_1 = Rx*r_300_2 - 3*r_200_2;
         double r_040_1 = Ry*r_030_2 - 3*r_020_2;
         double r_004_1 = Rz*r_003_2 - 3*r_002_2;
         double r_310_1 = Ry*r_300_2;
         double r_130_1 = Rx*r_030_2;
         double r_301_1 = Rz*r_300_2;
         double r_103_1 = Rx*r_003_2;
         double r_031_1 = Rz*r_030_2;
         double r_013_1 = Ry*r_003_2;
         double r_220_1 = Rx*r_120_2 - r_020_2;
         double r_202_1 = Rx*r_102_2 - r_002_2;
         double r_022_1 = Ry*r_012_2 - r_002_2;
         pOut[0] = Rx*r_400_1 - 4*r_300_1;
         pOut[1] = Ry*r_040_1 - 4*r_030_1;
         pOut[2] = Rz*r_004_1 - 4*r_003_1;
         pOut[3] = Rx*r_040_1;
         pOut[4] = Rx*r_004_1;
         pOut[5] = Ry*r_400_1;
         pOut[6] = Ry*r_004_1;
         pOut[7] = Rz*r_400_1;
         pOut[8] = Rz*r_040_1;
         pOut[9] = Ry*r_310_1 - r_300_1;
         pOut[10] = Rz*r_301_1 - r_300_1;
         pOut[11] = Rx*r_130_1 - r_030_1;
         pOut[12] = Rz*r_031_1 - r_030_1;
         pOut[13] = Rx*r_103_1 - r_003_1;
         pOut[14] = Ry*r_013_1 - r_003_1;
         pOut[15] = Rz*r_310_1;
         pOut[16] = Rz*r_130_1;
         pOut[17] = Ry*r_103_1;
         pOut[18] = Rx*r_022_1;
         pOut[19] = Ry*r_202_1;
         pOut[20] = Rz*r_220_1;
         return;
      }
      case 6: {
         double r_100_5 = Rx*pIn[6];
         double r_010_5 = Ry*pIn[6];
         double r_001_5 = Rz*pIn[6];
         double r_100_4 = Rx*pIn[5];
         double r_010_4 = Ry*pIn[5];
         double r_001_4 = Rz*pIn[5];
         double r_200_4 = Rx*r_100_5 - pIn[5];
         double r_020_4 = Ry*r_010_5 - pIn[5];
         double r_002_4 = Rz*r_001_5 - pIn[5];
         double r_100_3 = Rx*pIn[4];
         double r_010_3 = Ry*pIn[4];
         double r_001_3 = Rz*pIn[4];
         double r_200_3 = Rx*r_100_4 - pIn[4];
         double r_020_3 = Ry*r_010_4 - pIn[4];
         double r_002_3 = Rz*r_001_4 - pIn[4];
         double r_300_3 = Rx*r_200_4 - 2*r_100_4;
         double r_030_3 = Ry*r_020_4 - 2*r_010_4;
         double r_003_3 = Rz*r_002_4 - 2*r_001_4;
         double r_012_3 = Ry*r_002_4;
         double r_200_2 = Rx*r_100_3 - pIn[3];
         double r_020_2 = Ry*r_010_3 - pIn[3];
         double r_002_2 = Rz*r_001_3 - pIn[3];
         double r_300_2 = Rx*r_200_3 - 2*r_100_3;
         double r_030_2 = Ry*r_020_3 - 2*r_010_3;
         double r_003_2 = Rz*r_002_3 - 2*r_001_3;
         double r_012_2 = Ry*r_002_3;
         double r_400_2 = Rx*r_300_3 - 3*r_200_3;
         double r_040_2 = Ry*r_030_3 - 3*r_020_3;
         double r_004_2 = Rz*r_003_3 - 3*r_002_3;
         double r_310_2 = Ry*r_300_3;
         double r_130_2 = Rx*r_030_3;
         double r_301_2 = Rz*r_300_3;
         double r_103_2 = Rx*r_003_3;
         double r_031_2 = Rz*r_030_3;
         double r_013_2 = Ry*r_003_3;
         double r_022_2 = Ry*r_012_3 - r_002_3;
         double r_400_1 = Rx*r_300_2 - 3*r_200_2;
         double r_040_1 = Ry*r_030_2 - 3*r_020_2;
         double r_004_1 = Rz*r_003_2 - 3*r_002_2;
         double r_310_1 = Ry*r_300_2;
         double r_301_1 = Rz*r_300_2;
         double r_031_1 = Rz*r_030_2;
         double r_022_1 = Ry*r_012_2 - r_002_2;
         double r_500_1 = Rx*r_400_2 - 4*r_300_2;
         double r_050_1 = Ry*r_040_2 - 4*r_030_2;
         double r_005_1 = Rz*r_004_2 - 4*r_003_2;
         double r_140_1 = Rx*r_040_2;
         double r_104_1 = Rx*r_004_2;
         double r_410_1 = Ry*r_400_2;
         double r_014_1 = Ry*r_004_2;
         double r_401_1 = Rz*r_400_2;
         double r_041_1 = Rz*r_040_2;
         double r_320_1 = Ry*r_310_2 - r_300_2;
         double r_302_1 = Rz*r_301_2 - r_300_2;
         double r_230_1 = Rx*r_130_2 - r_030_2;
         double r_032_1 = Rz*r_031_2 - r_030_2;
         double r_203_1 = Rx*r_103_2 - r_003_2;
         double r_023_1 = Ry*r_013_2 - r_003_2;
         double r_122_1 = Rx*r_022_2;
         pOut[0] = Rx*r_500_1 - 5*r_400_1;
         pOut[1] = Ry*r_050_1 - 5*r_040_1;
         pOut[2] = Rz*r_005_1 - 5*r_004_1;
         pOut[3] = Ry*r_500_1;
         pOut[4] = Rx*r_050_1;
         pOut[5] = Rz*r_500_1;
         pOut[6] = Rx*r_005_1;
         pOut[7] = Rz*r_050_1;
         pOut[8] = Ry*r_005_1;
         pOut[9] = Ry*r_410_1 - r_400_1;
         pOut[10] = Rz*r_401_1 - r_400_1;
         pOut[11] = Rx*r_140_1 - r_040_1;
         pOut[12] = Rx*r_104_1 - r_004_1;
         pOut[13] = Rz*r_041_1 - r_040_1;
         pOut[14] = Ry*r_014_1 - r_004_1;
         pOut[15] = Ry*r_320_1 - 2*r_310_1;
         pOut[16] = Rz*r_302_1 - 2*r_301_1;
         pOut[17] = Rz*r_032_1 - 2*r_031_1;
         pOut[18] = Ry*r_104_1;
         pOut[19] = Rz*r_140_1;
         pOut[20] = Rz*r_410_1;
         pOut[21] = Ry*r_302_1;
         pOut[22] = Rx*r_032_1;
         pOut[23] = Rz*r_320_1;
         pOut[24] = Rx*r_023_1;
         pOut[25] = Rz*r_230_1;
         pOut[26] = Ry*r_203_1;
         pOut[27] = Rx*r_122_1 - r_022_1;
         return;
      }
      case 7: {
         double r_100_6 = Rx*pIn[7];
         double r_010_6 = Ry*pIn[7];
         double r_001_6 = Rz*pIn[7];
         double r_100_5 = Rx*pIn[6];
         double r_010_5 = Ry*pIn[6];
         double r_001_5 = Rz*pIn[6];
         double r_200_5 = Rx*r_100_6 - pIn[6];
         double r_020_5 = Ry*r_010_6 - pIn[6];
         double r_002_5 = Rz*r_001_6 - pIn[6];
         double r_100_4 = Rx*pIn[5];
         double r_010_4 = Ry*pIn[5];
         double r_001_4 = Rz*pIn[5];
         double r_200_4 = Rx*r_100_5 - pIn[5];
         double r_020_4 = Ry*r_010_5 - pIn[5];
         double r_002_4 = Rz*r_001_5 - pIn[5];
         double r_300_4 = Rx*r_200_5 - 2*r_100_5;
         double r_030_4 = Ry*r_020_5 - 2*r_010_5;
         double r_003_4 = Rz*r_002_5 - 2*r_001_5;
         double r_100_3 = Rx*pIn[4];
         double r_010_3 = Ry*pIn[4];
         double r_001_3 = Rz*pIn[4];
         double r_200_3 = Rx*r_100_4 - pIn[4];
         double r_020_3 = Ry*r_010_4 - pIn[4];
         double r_002_3 = Rz*r_001_4 - pIn[4];
         double r_300_3 = Rx*r_200_4 - 2*r_100_4;
         double r_030_3 = Ry*r_020_4 - 2*r_010_4;
         double r_003_3 = Rz*r_002_4 - 2*r_001_4;
         double r_400_3 = Rx*r_300_4 - 3*r_200_4;
         double r_040_3 = Ry*r_030_4 - 3*r_020_4;
         double r_004_3 = Rz*r_003_4 - 3*r_002_4;
         double r_310_3 = Ry*r_300_4;
         double r_301_3 = Rz*r_300_4;
         double r_031_3 = Rz*r_030_4;
         double r_013_3 = Ry*r_003_4;
         double r_300_2 = Rx*r_200_3 - 2*r_100_3;
         double r_030_2 = Ry*r_020_3 - 2*r_010_3;
         double r_003_2 = Rz*r_002_3 - 2*r_001_3;
         double r_400_2 = Rx*r_300_3 - 3*r_200_3;
         double r_040_2 = Ry*r_030_3 - 3*r_020_3;
         double r_004_2 = Rz*r_003_3 - 3*r_002_3;
         double r_310_2 = Ry*r_300_3;
         double r_301_2 = Rz*r_300_3;
         double r_031_2 = Rz*r_030_3;
         double r_013_2 = Ry*r_003_3;
         double r_500_2 = Rx*r_400_3 - 4*r_300_3;
         double r_050_2 = Ry*r_040_3 - 4*r_030_3;
         double r_005_2 = Rz*r_004_3 - 4*r_003_3;
         double r_140_2 = Rx*r_040_3;
         double r_104_2 = Rx*r_004_3;
         double r_410_2 = Ry*r_400_3;
         double r_014_2 = Ry*r_004_3;
         double r_401_2 = Rz*r_400_3;
         double r_041_2 = Rz*r_040_3;
         double r_320_2 = Ry*r_310_3 - r_300_3;
         double r_302_2 = Rz*r_301_3 - r_300_3;
         double r_032_2 = Rz*r_031_3 - r_030_3;
         double r_023_2 = Ry*r_013_3 - r_003_3;
         double r_500_1 = Rx*r_400_2 - 4*r_300_2;
         double r_050_1 = Ry*r_040_2 - 4*r_030_2;
         double r_005_1 = Rz*r_004_2 - 4*r_003_2;
         double r_140_1 = Rx*r_040_2;
         double r_104_1 = Rx*r_004_2;
         double r_410_1 = Ry*r_400_2;
         double r_014_1 = Ry*r_004_2;
         double r_401_1 = Rz*r_400_2;
         double r_041_1 = Rz*r_040_2;
         double r_302_1 = Rz*r_301_2 - r_300_2;
         double r_032_1 = Rz*r_031_2 - r_030_2;
         double r_023_1 = Ry*r_013_2 - r_003_2;
         double r_600_1 = Rx*r_500_2 - 5*r_400_2;
         double r_060_1 = Ry*r_050_2 - 5*r_040_2;
         double r_006_1 = Rz*r_005_2 - 5*r_004_2;
         double r_510_1 = Ry*r_500_2;
         double r_150_1 = Rx*r_050_2;
         double r_501_1 = Rz*r_500_2;
         double r_105_1 = Rx*r_005_2;
         double r_051_1 = Rz*r_050_2;
         double r_015_1 = Ry*r_005_2;
         double r_420_1 = Ry*r_410_2 - r_400_2;
         double r_402_1 = Rz*r_401_2 - r_400_2;
         double r_240_1 = Rx*r_140_2 - r_040_2;
         double r_204_1 = Rx*r_104_2 - r_004_2;
         double r_042_1 = Rz*r_041_2 - r_040_2;
         double r_024_1 = Ry*r_014_2 - r_004_2;
         double r_330_1 = Ry*r_320_2 - 2*r_310_2;
         double r_303_1 = Rz*r_302_2 - 2*r_301_2;
         double r_033_1 = Rz*r_032_2 - 2*r_031_2;
         double r_312_1 = Ry*r_302_2;
         double r_132_1 = Rx*r_032_2;
         double r_123_1 = Rx*r_023_2;
         pOut[0] = Rx*r_600_1 - 6*r_500_1;
         pOut[1] = Ry*r_060_1 - 6*r_050_1;
         pOut[2] = Rz*r_006_1 - 6*r_005_1;
         pOut[3] = Rx*r_060_1;
         pOut[4] = Rx*r_006_1;
         pOut[5] = Ry*r_600_1;
         pOut[6] = Ry*r_006_1;
         pOut[7] = Rz*r_600_1;
         pOut[8] = Rz*r_060_1;
         pOut[9] = Ry*r_510_1 - r_500_1;
         pOut[10] = Rz*r_501_1 - r_500_1;
         pOut[11] = Rx*r_150_1 - r_050_1;
         pOut[12] = Rz*r_051_1 - r_050_1;
         pOut[13] = Rx*r_105_1 - r_005_1;
         pOut[14] = Ry*r_015_1 - r_005_1;
         pOut[15] = Rx*r_240_1 - 2*r_140_1;
         pOut[16] = Rx*r_204_1 - 2*r_104_1;
         pOut[17] = Ry*r_420_1 - 2*r_410_1;
         pOut[18] = Ry*r_024_1 - 2*r_014_1;
         pOut[19] = Rz*r_402_1 - 2*r_401_1;
         pOut[20] = Rz*r_042_1 - 2*r_041_1;
         pOut[21] = Rz*r_510_1;
         pOut[22] = Rz*r_150_1;
         pOut[23] = Ry*r_105_1;
         pOut[24] = Rx*r_042_1;
         pOut[25] = Rx*r_024_1;
         pOut[26] = Ry*r_402_1;
         pOut[27] = Ry*r_204_1;
         pOut[28] = Rz*r_420_1;
         pOut[29] = Rz*r_240_1;
         pOut[30] = Rz*r_330_1;
         pOut[31] = Ry*r_303_1;
         pOut[32] = Rx*r_033_1;
         pOut[33] = Ry*r_312_1 - r_302_1;
         pOut[34] = Rx*r_132_1 - r_032_1;
         pOut[35] = Rx*r_123_1 - r_023_1;
         return;
      }
      case 8: {
         double r_100_7 = Rx*pIn[8];
         double r_010_7 = Ry*pIn[8];
         double r_001_7 = Rz*pIn[8];
         double r_100_6 = Rx*pIn[7];
         double r_010_6 = Ry*pIn[7];
         double r_001_6 = Rz*pIn[7];
         double r_200_6 = Rx*r_100_7 - pIn[7];
         double r_020_6 = Ry*r_010_7 - pIn[7];
         double r_002_6 = Rz*r_001_7 - pIn[7];
         double r_100_5 = Rx*pIn[6];
         double r_010_5 = Ry*pIn[6];
         double r_001_5 = Rz*pIn[6];
         double r_200_5 = Rx*r_100_6 - pIn[6];
         double r_020_5 = Ry*r_010_6 - pIn[6];
         double r_002_5 = Rz*r_001_6 - pIn[6];
         double r_300_5 = Rx*r_200_6 - 2*r_100_6;
         double r_030_5 = Ry*r_020_6 - 2*r_010_6;
         double r_003_5 = Rz*r_002_6 - 2*r_001_6;
         double r_100_4 = Rx*pIn[5];
         double r_010_4 = Ry*pIn[5];
         double r_001_4 = Rz*pIn[5];
         double r_200_4 = Rx*r_100_5 - pIn[5];
         double r_020_4 = Ry*r_010_5 - pIn[5];
         double r_002_4 = Rz*r_001_5 - pIn[5];
         double r_300_4 = Rx*r_200_5 - 2*r_100_5;
         double r_030_4 = Ry*r_020_5 - 2*r_010_5;
         double r_003_4 = Rz*r_002_5 - 2*r_001_5;
         double r_400_4 = Rx*r_300_5 - 3*r_200_5;
         double r_040_4 = Ry*r_030_5 - 3*r_020_5;
         double r_004_4 = Rz*r_003_5 - 3*r_002_5;
         double r_310_4 = Ry*r_300_5;
         double r_301_4 = Rz*r_300_5;
         double r_031_4 = Rz*r_030_5;
         double r_200_3 = Rx*r_100_4 - pIn[4];
         double r_020_3 = Ry*r_010_4 - pIn[4];
         double r_002_3 = Rz*r_001_4 - pIn[4];
         double r_300_3 = Rx*r_200_4 - 2*r_100_4;
         double r_030_3 = Ry*r_020_4 - 2*r_010_4;
         double r_003_3 = Rz*r_002_4 - 2*r_001_4;
         double r_400_3 = Rx*r_300_4 - 3*r_200_4;
         double r_040_3 = Ry*r_030_4 - 3*r_020_4;
         double r_004_3 = Rz*r_003_4 - 3*r_002_4;
         double r_310_3 = Ry*r_300_4;
         double r_301_3 = Rz*r_300_4;
         double r_031_3 = Rz*r_030_4;
         double r_500_3 = Rx*r_400_4 - 4*r_300_4;
         double r_050_3 = Ry*r_040_4 - 4*r_030_4;
         double r_005_3 = Rz*r_004_4 - 4*r_003_4;
         double r_140_3 = Rx*r_040_4;
         double r_104_3 = Rx*r_004_4;
         double r_410_3 = Ry*r_400_4;
         double r_014_3 = Ry*r_004_4;
         double r_401_3 = Rz*r_400_4;
         double r_041_3 = Rz*r_040_4;
         double r_320_3 = Ry*r_310_4 - r_300_4;
         double r_302_3 = Rz*r_301_4 - r_300_4;
         double r_032_3 = Rz*r_031_4 - r_030_4;
         double r_400_2 = Rx*r_300_3 - 3*r_200_3;
         double r_040_2 = Ry*r_030_3 - 3*r_020_3;
         double r_004_2 = Rz*r_003_3 - 3*r_002_3;
         double r_310_2 = Ry*r_300_3;
         double r_301_2 = Rz*r_300_3;
         double r_031_2 = Rz*r_030_3;
         double r_500_2 = Rx*r_400_3 - 4*r_300_3;
         double r_050_2 = Ry*r_040_3 - 4*r_030_3;
         double r_005_2 = Rz*r_004_3 - 4*r_003_3;
         double r_140_2 = Rx*r_040_3;
         double r_104_2 = Rx*r_004_3;
         double r_410_2 = Ry*r_400_3;
         double r_014_2 = Ry*r_004_3;
         double r_401_2 = Rz*r_400_3;
         double r_041_2 = Rz*r_040_3;
         double r_320_2 = Ry*r_310_3 - r_300_3;
         double r_302_2 = Rz*r_301_3 - r_300_3;
         double r_032_2 = Rz*r_031_3 - r_030_3;
         double r_600_2 = Rx*r_500_3 - 5*r_400_3;
         double r_060_2 = Ry*r_050_3 - 5*r_040_3;
         double r_006_2 = Rz*r_005_3 - 5*r_004_3;
         double r_510_2 = Ry*r_500_3;
         double r_150_2 = Rx*r_050_3;
         double r_501_2 = Rz*r_500_3;
         double r_105_2 = Rx*r_005_3;
         double r_051_2 = Rz*r_050_3;
         double r_015_2 = Ry*r_005_3;
         double r_420_2 = Ry*r_410_3 - r_400_3;
         double r_402_2 = Rz*r_401_3 - r_400_3;
         double r_240_2 = Rx*r_140_3 - r_040_3;
         double r_204_2 = Rx*r_104_3 - r_004_3;
         double r_042_2 = Rz*r_041_3 - r_040_3;
         double r_024_2 = Ry*r_014_3 - r_004_3;
         double r_330_2 = Ry*r_320_3 - 2*r_310_3;
         double r_303_2 = Rz*r_302_3 - 2*r_301_3;
         double r_033_2 = Rz*r_032_3 - 2*r_031_3;
         double r_600_1 = Rx*r_500_2 - 5*r_400_2;
         double r_060_1 = Ry*r_050_2 - 5*r_040_2;
         double r_006_1 = Rz*r_005_2 - 5*r_004_2;
         double r_510_1 = Ry*r_500_2;
         double r_150_1 = Rx*r_050_2;
         double r_501_1 = Rz*r_500_2;
         double r_105_1 = Rx*r_005_2;
         double r_051_1 = Rz*r_050_2;
         double r_015_1 = Ry*r_005_2;
         double r_402_1 = Rz*r_401_2 - r_400_2;
         double r_240_1 = Rx*r_140_2 - r_040_2;
         double r_204_1 = Rx*r_104_2 - r_004_2;
         double r_042_1 = Rz*r_041_2 - r_040_2;
         double r_024_1 = Ry*r_014_2 - r_004_2;
         double r_330_1 = Ry*r_320_2 - 2*r_310_2;
         double r_303_1 = Rz*r_302_2 - 2*r_301_2;
         double r_033_1 = Rz*r_032_2 - 2*r_031_2;
         double r_700_1 = Rx*r_600_2 - 6*r_500_2;
         double r_070_1 = Ry*r_060_2 - 6*r_050_2;
         double r_007_1 = Rz*r_006_2 - 6*r_005_2;
         double r_160_1 = Rx*r_060_2;
         double r_106_1 = Rx*r_006_2;
         double r_610_1 = Ry*r_600_2;
         double r_016_1 = Ry*r_006_2;
         double r_601_1 = Rz*r_600_2;
         double r_061_1 = Rz*r_060_2;
         double r_520_1 = Ry*r_510_2 - r_500_2;
         double r_502_1 = Rz*r_501_2 - r_500_2;
         double r_250_1 = Rx*r_150_2 - r_050_2;
         double r_052_1 = Rz*r_051_2 - r_050_2;
         double r_205_1 = Rx*r_105_2 - r_005_2;
         double r_025_1 = Ry*r_015_2 - r_005_2;
         double r_340_1 = Rx*r_240_2 - 2*r_140_2;
         double r_304_1 = Rx*r_204_2 - 2*r_104_2;
         double r_430_1 = Ry*r_420_2 - 2*r_410_2;
         double r_034_1 = Ry*r_024_2 - 2*r_014_2;
         double r_403_1 = Rz*r_402_2 - 2*r_401_2;
         double r_043_1 = Rz*r_042_2 - 2*r_041_2;
         double r_142_1 = Rx*r_042_2;
         double r_124_1 = Rx*r_024_2;
         double r_412_1 = Ry*r_402_2;
         double r_331_1 = Rz*r_330_2;
         double r_313_1 = Ry*r_303_2;
         double r_133_1 = Rx*r_033_2;
         pOut[0] = Rx*r_700_1 - 7*r_600_1;
         pOut[1] = Ry*r_070_1 - 7*r_060_1;
         pOut[2] = Rz*r_007_1 - 7*r_006_1;
         pOut[3] = Ry*r_700_1;
         pOut[4] = Rx*r_070_1;
         pOut[5] = Rz*r_700_1;
         pOut[6] = Rx*r_007_1;
         pOut[7] = Rz*r_070_1;
         pOut[8] = Ry*r_007_1;
         pOut[9] = Ry*r_610_1 - r_600_1;
         pOut[10] = Rz*r_601_1 - r_600_1;
         pOut[11] = Rx*r_160_1 - r_060_1;
         pOut[12] = Rx*r_106_1 - r_006_1;
         pOut[13] = Rz*r_061_1 - r_060_1;
         pOut[14] = Ry*r_016_1 - r_006_1;
         pOut[15] = Ry*r_520_1 - 2*r_510_1;
         pOut[16] = Rx*r_250_1 - 2*r_150_1;
         pOut[17] = Rz*r_502_1 - 2*r_501_1;
         pOut[18] = Rx*r_205_1 - 2*r_105_1;
         pOut[19] = Rz*r_052_1 - 2*r_051_1;
         pOut[20] = Ry*r_025_1 - 2*r_015_1;
         pOut[21] = Rx*r_340_1 - 3*r_240_1;
         pOut[22] = Rx*r_304_1 - 3*r_204_1;
         pOut[23] = Ry*r_034_1 - 3*r_024_1;
         pOut[24] = Ry*r_106_1;
         pOut[25] = Rz*r_160_1;
         pOut[26] = Rz*r_610_1;
         pOut[27] = Ry*r_502_1;
         pOut[28] = Rx*r_052_1;
         pOut[29] = Rz*r_520_1;
         pOut[30] = Rx*r_025_1;
         pOut[31] = Rz*r_250_1;
         pOut[32] = Ry*r_205_1;
         pOut[33] = Ry*r_304_1;
         pOut[34] = Rx*r_034_1;
         pOut[35] = Rz*r_340_1;
         pOut[36] = Rx*r_043_1;
         pOut[37] = Rz*r_430_1;
         pOut[38] = Ry*r_403_1;
         pOut[39] = Ry*r_412_1 - r_402_1;
         pOut[40] = Rx*r_142_1 - r_042_1;
         pOut[41] = Rx*r_124_1 - r_024_1;
         pOut[42] = Rz*r_331_1 - r_330_1;
         pOut[43] = Ry*r_313_1 - r_303_1;
         pOut[44] = Rx*r_133_1 - r_033_1;
         return;
      }
      case 9: {
         double r_100_8 = Rx*pIn[9];
         double r_010_8 = Ry*pIn[9];
         double r_001_8 = Rz*pIn[9];
         double r_100_7 = Rx*pIn[8];
         double r_010_7 = Ry*pIn[8];
         double r_001_7 = Rz*pIn[8];
         double r_200_7 = Rx*r_100_8 - pIn[8];
         double r_020_7 = Ry*r_010_8 - pIn[8];
         double r_002_7 = Rz*r_001_8 - pIn[8];
         double r_100_6 = Rx*pIn[7];
         double r_010_6 = Ry*pIn[7];
         double r_001_6 = Rz*pIn[7];
         double r_200_6 = Rx*r_100_7 - pIn[7];
         double r_020_6 = Ry*r_010_7 - pIn[7];
         double r_002_6 = Rz*r_001_7 - pIn[7];
         double r_300_6 = Rx*r_200_7 - 2*r_100_7;
         double r_030_6 = Ry*r_020_7 - 2*r_010_7;
         double r_003_6 = Rz*r_002_7 - 2*r_001_7;
         double r_100_5 = Rx*pIn[6];
         double r_010_5 = Ry*pIn[6];
         double r_001_5 = Rz*pIn[6];
         double r_200_5 = Rx*r_100_6 - pIn[6];
         double r_020_5 = Ry*r_010_6 - pIn[6];
         double r_002_5 = Rz*r_001_6 - pIn[6];
         double r_300_5 = Rx*r_200_6 - 2*r_100_6;
         double r_030_5 = Ry*r_020_6 - 2*r_010_6;
         double r_003_5 = Rz*r_002_6 - 2*r_001_6;
         double r_400_5 = Rx*r_300_6 - 3*r_200_6;
         double r_040_5 = Ry*r_030_6 - 3*r_020_6;
         double r_004_5 = Rz*r_003_6 - 3*r_002_6;
         double r_310_5 = Ry*r_300_6;
         double r_100_4 = Rx*pIn[5];
         double r_010_4 = Ry*pIn[5];
         double r_001_4 = Rz*pIn[5];
         double r_200_4 = Rx*r_100_5 - pIn[5];
         double r_020_4 = Ry*r_010_5 - pIn[5];
         double r_002_4 = Rz*r_001_5 - pIn[5];
         double r_300_4 = Rx*r_200_5 - 2*r_100_5;
         double r_030_4 = Ry*r_020_5 - 2*r_010_5;
         double r_003_4 = Rz*r_002_5 - 2*r_001_5;
         double r_400_4 = Rx*r_300_5 - 3*r_200_5;
         double r_040_4 = Ry*r_030_5 - 3*r_020_5;
         double r_004_4 = Rz*r_003_5 - 3*r_002_5;
         double r_310_4 = Ry*r_300_5;
         double r_500_4 = Rx*r_400_5 - 4*r_300_5;
         double r_050_4 = Ry*r_040_5 - 4*r_030_5;
         double r_005_4 = Rz*r_004_5 - 4*r_003_5;
         double r_140_4 = Rx*r_040_5;
         double r_104_4 = Rx*r_004_5;
         double r_410_4 = Ry*r_400_5;
         double r_014_4 = Ry*r_004_5;
         double r_401_4 = Rz*r_400_5;
         double r_041_4 = Rz*r_040_5;
         double r_320_4 = Ry*r_310_5 - r_300_5;
         double r_300_3 = Rx*r_200_4 - 2*r_100_4;
         double r_030_3 = Ry*r_020_4 - 2*r_010_4;
         double r_003_3 = Rz*r_002_4 - 2*r_001_4;
         double r_400_3 = Rx*r_300_4 - 3*r_200_4;
         double r_040_3 = Ry*r_030_4 - 3*r_020_4;
         double r_004_3 = Rz*r_003_4 - 3*r_002_4;
         double r_310_3 = Ry*r_300_4;
         double r_500_3 = Rx*r_400_4 - 4*r_300_4;
         double r_050_3 = Ry*r_040_4 - 4*r_030_4;
         double r_005_3 = Rz*r_004_4 - 4*r_003_4;
         double r_140_3 = Rx*r_040_4;
         double r_104_3 = Rx*r_004_4;
         double r_410_3 = Ry*r_400_4;
         double r_014_3 = Ry*r_004_4;
         double r_401_3 = Rz*r_400_4;
         double r_041_3 = Rz*r_040_4;
         double r_320_3 = Ry*r_310_4 - r_300_4;
         double r_600_3 = Rx*r_500_4 - 5*r_400_4;
         double r_060_3 = Ry*r_050_4 - 5*r_040_4;
         double r_006_3 = Rz*r_005_4 - 5*r_004_4;
         double r_510_3 = Ry*r_500_4;
         double r_150_3 = Rx*r_050_4;
         double r_501_3 = Rz*r_500_4;
         double r_105_3 = Rx*r_005_4;
         double r_051_3 = Rz*r_050_4;
         double r_015_3 = Ry*r_005_4;
         double r_420_3 = Ry*r_410_4 - r_400_4;
         double r_402_3 = Rz*r_401_4 - r_400_4;
         double r_240_3 = Rx*r_140_4 - r_040_4;
         double r_204_3 = Rx*r_104_4 - r_004_4;
         double r_042_3 = Rz*r_041_4 - r_040_4;
         double r_024_3 = Ry*r_014_4 - r_004_4;
         double r_330_3 = Ry*r_320_4 - 2*r_310_4;
         double r_500_2 = Rx*r_400_3 - 4*r_300_3;
         double r_050_2 = Ry*r_040_3 - 4*r_030_3;
         double r_005_2 = Rz*r_004_3 - 4*r_003_3;
         double r_140_2 = Rx*r_040_3;
         double r_104_2 = Rx*r_004_3;
         double r_410_2 = Ry*r_400_3;
         double r_014_2 = Ry*r_004_3;
         double r_401_2 = Rz*r_400_3;
         double r_041_2 = Rz*r_040_3;
         double r_600_2 = Rx*r_500_3 - 5*r_400_3;
         double r_060_2 = Ry*r_050_3 - 5*r_040_3;
         double r_006_2 = Rz*r_005_3 - 5*r_004_3;
         double r_510_2 = Ry*r_500_3;
         double r_150_2 = Rx*r_050_3;
         double r_501_2 = Rz*r_500_3;
         double r_105_2 = Rx*r_005_3;
         double r_051_2 = Rz*r_050_3;
         double r_015_2 = Ry*r_005_3;
         double r_420_2 = Ry*r_410_3 - r_400_3;
         double r_402_2 = Rz*r_401_3 - r_400_3;
         double r_240_2 = Rx*r_140_3 - r_040_3;
         double r_204_2 = Rx*r_104_3 - r_004_3;
         double r_042_2 = Rz*r_041_3 - r_040_3;
         double r_024_2 = Ry*r_014_3 - r_004_3;
         double r_330_2 = Ry*r_320_3 - 2*r_310_3;
         double r_700_2 = Rx*r_600_3 - 6*r_500_3;
         double r_070_2 = Ry*r_060_3 - 6*r_050_3;
         double r_007_2 = Rz*r_006_3 - 6*r_005_3;
         double r_160_2 = Rx*r_060_3;
         double r_106_2 = Rx*r_006_3;
         double r_610_2 = Ry*r_600_3;
         double r_016_2 = Ry*r_006_3;
         double r_601_2 = Rz*r_600_3;
         double r_061_2 = Rz*r_060_3;
         double r_520_2 = Ry*r_510_3 - r_500_3;
         double r_502_2 = Rz*r_501_3 - r_500_3;
         double r_250_2 = Rx*r_150_3 - r_050_3;
         double r_052_2 = Rz*r_051_3 - r_050_3;
         double r_205_2 = Rx*r_105_3 - r_005_3;
         double r_025_2 = Ry*r_015_3 - r_005_3;
         double r_340_2 = Rx*r_240_3 - 2*r_140_3;
         double r_304_2 = Rx*r_204_3 - 2*r_104_3;
         double r_430_2 = Ry*r_420_3 - 2*r_410_3;
         double r_034_2 = Ry*r_024_3 - 2*r_014_3;
         double r_403_2 = Rz*r_402_3 - 2*r_401_3;
         double r_043_2 = Rz*r_042_3 - 2*r_041_3;
         double r_331_2 = Rz*r_330_3;
         double r_700_1 = Rx*r_600_2 - 6*r_500_2;
         double r_070_1 = Ry*r_060_2 - 6*r_050_2;
         double r_007_1 = Rz*r_006_2 - 6*r_005_2;
         double r_160_1 = Rx*r_060_2;
         double r_106_1 = Rx*r_006_2;
         double r_610_1 = Ry*r_600_2;
         double r_016_1 = Ry*r_006_2;
         double r_601_1 = Rz*r_600_2;
         double r_061_1 = Rz*r_060_2;
         double r_520_1 = Ry*r_510_2 - r_500_2;
         double r_502_1 = Rz*r_501_2 - r_500_2;
         double r_250_1 = Rx*r_150_2 - r_050_2;
         double r_052_1 = Rz*r_051_2 - r_050_2;
         double r_205_1 = Rx*r_105_2 - r_005_2;
         double r_025_1 = Ry*r_015_2 - r_005_2;
         double r_340_1 = Rx*r_240_2 - 2*r_140_2;
         double r_304_1 = Rx*r_204_2 - 2*r_104_2;
         double r_430_1 = Ry*r_420_2 - 2*r_410_2;
         double r_034_1 = Ry*r_024_2 - 2*r_014_2;
         double r_403_1 = Rz*r_402_2 - 2*r_401_2;
         double r_043_1 = Rz*r_042_2 - 2*r_041_2;
         double r_331_1 = Rz*r_330_2;
         double r_800_1 = Rx*r_700_2 - 7*r_600_2;
         double r_080_1 = Ry*r_070_2 - 7*r_060_2;
         double r_008_1 = Rz*r_007_2 - 7*r_006_2;
         double r_710_1 = Ry*r_700_2;
         double r_170_1 = Rx*r_070_2;
         double r_701_1 = Rz*r_700_2;
         double r_107_1 = Rx*r_007_2;
         double r_071_1 = Rz*r_070_2;
         double r_017_1 = Ry*r_007_2;
         double r_620_1 = Ry*r_610_2 - r_600_2;
         double r_602_1 = Rz*r_601_2 - r_600_2;
         double r_260_1 = Rx*r_160_2 - r_060_2;
         double r_206_1 = Rx*r_106_2 - r_006_2;
         double r_062_1 = Rz*r_061_2 - r_060_2;
         double r_026_1 = Ry*r_016_2 - r_006_2;
         double r_530_1 = Ry*r_520_2 - 2*r_510_2;
         double r_350_1 = Rx*r_250_2 - 2*r_150_2;
         double r_503_1 = Rz*r_502_2 - 2*r_501_2;
         double r_305_1 = Rx*r_205_2 - 2*r_105_2;
         double r_053_1 = Rz*r_052_2 - 2*r_051_2;
         double r_035_1 = Ry*r_025_2 - 2*r_015_2;
         double r_440_1 = Rx*r_340_2 - 3*r_240_2;
         double r_404_1 = Rx*r_304_2 - 3*r_204_2;
         double r_044_1 = Ry*r_034_2 - 3*r_024_2;
         double r_512_1 = Ry*r_502_2;
         double r_152_1 = Rx*r_052_2;
         double r_125_1 = Rx*r_025_2;
         double r_314_1 = Ry*r_304_2;
         double r_134_1 = Rx*r_034_2;
         double r_341_1 = Rz*r_340_2;
         double r_143_1 = Rx*r_043_2;
         double r_431_1 = Rz*r_430_2;
         double r_413_1 = Ry*r_403_2;
         double r_332_1 = Rz*r_331_2 - r_330_2;
         pOut[0] = Rx*r_800_1 - 8*r_700_1;
         pOut[1] = Ry*r_080_1 - 8*r_070_1;
         pOut[2] = Rz*r_008_1 - 8*r_007_1;
         pOut[3] = Rx*r_080_1;
         pOut[4] = Rx*r_008_1;
         pOut[5] = Ry*r_800_1;
         pOut[6] = Ry*r_008_1;
         pOut[7] = Rz*r_800_1;
         pOut[8] = Rz*r_080_1;
         pOut[9] = Ry*r_710_1 - r_700_1;
         pOut[10] = Rz*r_701_1 - r_700_1;
         pOut[11] = Rx*r_170_1 - r_070_1;
         pOut[12] = Rz*r_071_1 - r_070_1;
         pOut[13] = Rx*r_107_1 - r_007_1;
         pOut[14] = Ry*r_017_1 - r_007_1;
         pOut[15] = Rx*r_260_1 - 2*r_160_1;
         pOut[16] = Rx*r_206_1 - 2*r_106_1;
         pOut[17] = Ry*r_620_1 - 2*r_610_1;
         pOut[18] = Ry*r_026_1 - 2*r_016_1;
         pOut[19] = Rz*r_602_1 - 2*r_601_1;
         pOut[20] = Rz*r_062_1 - 2*r_061_1;
         pOut[21] = Ry*r_530_1 - 3*r_520_1;
         pOut[22] = Rz*r_503_1 - 3*r_502_1;
         pOut[23] = Rx*r_350_1 - 3*r_250_1;
         pOut[24] = Rz*r_053_1 - 3*r_052_1;
         pOut[25] = Rx*r_305_1 - 3*r_205_1;
         pOut[26] = Ry*r_035_1 - 3*r_025_1;
         pOut[27] = Rz*r_710_1;
         pOut[28] = Rz*r_170_1;
         pOut[29] = Ry*r_107_1;
         pOut[30] = Rx*r_062_1;
         pOut[31] = Rx*r_026_1;
         pOut[32] = Ry*r_602_1;
         pOut[33] = Ry*r_206_1;
         pOut[34] = Rz*r_620_1;
         pOut[35] = Rz*r_260_1;
         pOut[36] = Rz*r_530_1;
         pOut[37] = Ry*r_503_1;
         pOut[38] = Rz*r_350_1;
         pOut[39] = Ry*r_305_1;
         pOut[40] = Rx*r_053_1;
         pOut[41] = Rx*r_035_1;
         pOut[42] = Rx*r_044_1;
         pOut[43] = Ry*r_404_1;
         pOut[44] = Rz*r_440_1;
         pOut[45] = Ry*r_512_1 - r_502_1;
         pOut[46] = Rx*r_152_1 - r_052_1;
         pOut[47] = Rx*r_125_1 - r_025_1;
         pOut[48] = Rz*r_341_1 - r_340_1;
         pOut[49] = Ry*r_314_1 - r_304_1;
         pOut[50] = Rz*r_431_1 - r_430_1;
         pOut[51] = Rx*r_134_1 - r_034_1;
         pOut[52] = Ry*r_413_1 - r_403_1;
         pOut[53] = Rx*r_143_1 - r_043_1;
         pOut[54] = Rz*r_332_1 - 2*r_331_1;
         return;
      }
      case 10: {
         double r_100_9 = Rx*pIn[10];
         double r_010_9 = Ry*pIn[10];
         double r_001_9 = Rz*pIn[10];
         double r_100_8 = Rx*pIn[9];
         double r_010_8 = Ry*pIn[9];
         double r_001_8 = Rz*pIn[9];
         double r_200_8 = Rx*r_100_9 - pIn[9];
         double r_020_8 = Ry*r_010_9 - pIn[9];
         double r_002_8 = Rz*r_001_9 - pIn[9];
         double r_100_7 = Rx*pIn[8];
         double r_010_7 = Ry*pIn[8];
         double r_001_7 = Rz*pIn[8];
         double r_200_7 = Rx*r_100_8 - pIn[8];
         double r_020_7 = Ry*r_010_8 - pIn[8];
         double r_002_7 = Rz*r_001_8 - pIn[8];
         double r_300_7 = Rx*r_200_8 - 2*r_100_8;
         double r_030_7 = Ry*r_020_8 - 2*r_010_8;
         double r_003_7 = Rz*r_002_8 - 2*r_001_8;
         double r_100_6 = Rx*pIn[7];
         double r_010_6 = Ry*pIn[7];
         double r_001_6 = Rz*pIn[7];
         double r_200_6 = Rx*r_100_7 - pIn[7];
         double r_020_6 = Ry*r_010_7 - pIn[7];
         double r_002_6 = Rz*r_001_7 - pIn[7];
         double r_300_6 = Rx*r_200_7 - 2*r_100_7;
         double r_030_6 = Ry*r_020_7 - 2*r_010_7;
         double r_003_6 = Rz*r_002_7 - 2*r_001_7;
         double r_400_6 = Rx*r_300_7 - 3*r_200_7;
         double r_040_6 = Ry*r_030_7 - 3*r_020_7;
         double r_004_6 = Rz*r_003_7 - 3*r_002_7;
         double r_100_5 = Rx*pIn[6];
         double r_010_5 = Ry*pIn[6];
         double r_001_5 = Rz*pIn[6];
         double r_200_5 = Rx*r_100_6 - pIn[6];
         double r_020_5 = Ry*r_010_6 - pIn[6];
         double r_002_5 = Rz*r_001_6 - pIn[6];
         double r_300_5 = Rx*r_200_6 - 2*r_100_6;
         double r_030_5 = Ry*r_020_6 - 2*r_010_6;
         double r_003_5 = Rz*r_002_6 - 2*r_001_6;
         double r_400_5 = Rx*r_300_6 - 3*r_200_6;
         double r_040_5 = Ry*r_030_6 - 3*r_020_6;
         double r_004_5 = Rz*r_003_6 - 3*r_002_6;
         double r_500_5 = Rx*r_400_6 - 4*r_300_6;
         double r_050_5 = Ry*r_040_6 - 4*r_030_6;
         double r_005_5 = Rz*r_004_6 - 4*r_003_6;
         double r_140_5 = Rx*r_040_6;
         double r_104_5 = Rx*r_004_6;
         double r_410_5 = Ry*r_400_6;
         double r_014_5 = Ry*r_004_6;
         double r_200_4 = Rx*r_100_5 - pIn[5];
         double r_020_4 = Ry*r_010_5 - pIn[5];
         double r_002_4 = Rz*r_001_5 - pIn[5];
         double r_300_4 = Rx*r_200_5 - 2*r_100_5;
         double r_030_4 = Ry*r_020_5 - 2*r_010_5;
         double r_003_4 = Rz*r_002_5 - 2*r_001_5;
         double r_400_4 = Rx*r_300_5 - 3*r_200_5;
         double r_040_4 = Ry*r_030_5 - 3*r_020_5;
         double r_004_4 = Rz*r_003_5 - 3*r_002_5;
         double r_500_4 = Rx*r_400_5 - 4*r_300_5;
         double r_050_4 = Ry*r_040_5 - 4*r_030_5;
         double r_005_4 = Rz*r_004_5 - 4*r_003_5;
         double r_140_4 = Rx*r_040_5;
         double r_104_4 = Rx*r_004_5;
         double r_410_4 = Ry*r_400_5;
         double r_014_4 = Ry*r_004_5;
         double r_600_4 = Rx*r_500_5 - 5*r_400_5;
         double r_060_4 = Ry*r_050_5 - 5*r_040_5;
         double r_006_4 = Rz*r_005_5 - 5*r_004_5;
         double r_510_4 = Ry*r_500_5;
         double r_150_4 = Rx*r_050_5;
         double r_501_4 = Rz*r_500_5;
         double r_105_4 = Rx*r_005_5;
         double r_051_4 = Rz*r_050_5;
         double r_015_4 = Ry*r_005_5;
         double r_420_4 = Ry*r_410_5 - r_400_5;
         double r_240_4 = Rx*r_140_5 - r_040_5;
         double r_204_4 = Rx*r_104_5 - r_004_5;
         double r_024_4 = Ry*r_014_5 - r_004_5;
         double r_400_3 = Rx*r_300_4 - 3*r_200_4;
         double r_040_3 = Ry*r_030_4 - 3*r_020_4;
         double r_004_3 = Rz*r_003_4 - 3*r_002_4;
         double r_500_3 = Rx*r_400_4 - 4*r_300_4;
         double r_050_3 = Ry*r_040_4 - 4*r_030_4;
         double r_005_3 = Rz*r_004_4 - 4*r_003_4;
         double r_140_3 = Rx*r_040_4;
         double r_104_3 = Rx*r_004_4;
         double r_410_3 = Ry*r_400_4;
         double r_014_3 = Ry*r_004_4;
         double r_600_3 = Rx*r_500_4 - 5*r_400_4;
         double r_060_3 = Ry*r_050_4 - 5*r_040_4;
         double r_006_3 = Rz*r_005_4 - 5*r_004_4;
         double r_510_3 = Ry*r_500_4;
         double r_150_3 = Rx*r_050_4;
         double r_501_3 = Rz*r_500_4;
         double r_105_3 = Rx*r_005_4;
         double r_051_3 = Rz*r_050_4;
         double r_015_3 = Ry*r_005_4;
         double r_420_3 = Ry*r_410_4 - r_400_4;
         double r_240_3 = Rx*r_140_4 - r_040_4;
         double r_204_3 = Rx*r_104_4 - r_004_4;
         double r_024_3 = Ry*r_014_4 - r_004_4;
         double r_700_3 = Rx*r_600_4 - 6*r_500_4;
         double r_070_3 = Ry*r_060_4 - 6*r_050_4;
         double r_007_3 = Rz*r_006_4 - 6*r_005_4;
         double r_160_3 = Rx*r_060_4;
         double r_106_3 = Rx*r_006_4;
         double r_610_3 = Ry*r_600_4;
         double r_016_3 = Ry*r_006_4;
         double r_601_3 = Rz*r_600_4;
         double r_061_3 = Rz*r_060_4;
         double r_520_3 = Ry*r_510_4 - r_500_4;
         double r_502_3 = Rz*r_501_4 - r_500_4;
         double r_250_3 = Rx*r_150_4 - r_050_4;
         double r_052_3 = Rz*r_051_4 - r_050_4;
         double r_205_3 = Rx*r_105_4 - r_005_4;
         double r_025_3 = Ry*r_015_4 - r_005_4;
         double r_340_3 = Rx*r_240_4 - 2*r_140_4;
         double r_304_3 = Rx*r_204_4 - 2*r_104_4;
         double r_430_3 = Ry*r_420_4 - 2*r_410_4;
         double r_034_3 = Ry*r_024_4 - 2*r_014_4;
         double r_600_2 = Rx*r_500_3 - 5*r_400_3;
         double r_060_2 = Ry*r_050_3 - 5*r_040_3;
         double r_006_2 = Rz*r_005_3 - 5*r_004_3;
         double r_510_2 = Ry*r_500_3;
         double r_150_2 = Rx*r_050_3;
         double r_501_2 = Rz*r_500_3;
         double r_105_2 = Rx*r_005_3;
         double r_051_2 = Rz*r_050_3;
         double r_015_2 = Ry*r_005_3;
         double r_240_2 = Rx*r_140_3 - r_040_3;
         double r_204_2 = Rx*r_104_3 - r_004_3;
         double r_024_2 = Ry*r_014_3 - r_004_3;
         double r_700_2 = Rx*r_600_3 - 6*r_500_3;
         double r_070_2 = Ry*r_060_3 - 6*r_050_3;
         double r_007_2 = Rz*r_006_3 - 6*r_005_3;
         double r_160_2 = Rx*r_060_3;
         double r_106_2 = Rx*r_006_3;
         double r_610_2 = Ry*r_600_3;
         double r_016_2 = Ry*r_006_3;
         double r_601_2 = Rz*r_600_3;
         double r_061_2 = Rz*r_060_3;
         double r_520_2 = Ry*r_510_3 - r_500_3;
         double r_502_2 = Rz*r_501_3 - r_500_3;
         double r_250_2 = Rx*r_150_3 - r_050_3;
         double r_052_2 = Rz*r_051_3 - r_050_3;
         double r_205_2 = Rx*r_105_3 - r_005_3;
         double r_025_2 = Ry*r_015_3 - r_005_3;
         double r_340_2 = Rx*r_240_3 - 2*r_140_3;
         double r_304_2 = Rx*r_204_3 - 2*r_104_3;
         double r_430_2 = Ry*r_420_3 - 2*r_410_3;
         double r_034_2 = Ry*r_024_3 - 2*r_014_3;
         double r_800_2 = Rx*r_700_3 - 7*r_600_3;
         double r_080_2 = Ry*r_070_3 - 7*r_060_3;
         double r_008_2 = Rz*r_007_3 - 7*r_006_3;
         double r_710_2 = Ry*r_700_3;
         double r_170_2 = Rx*r_070_3;
         double r_701_2 = Rz*r_700_3;
         double r_107_2 = Rx*r_007_3;
         double r_071_2 = Rz*r_070_3;
         double r_017_2 = Ry*r_007_3;
         double r_620_2 = Ry*r_610_3 - r_600_3;
         double r_602_2 = Rz*r_601_3 - r_600_3;
         double r_260_2 = Rx*r_160_3 - r_060_3;
         double r_206_2 = Rx*r_106_3 - r_006_3;
         double r_062_2 = Rz*r_061_3 - r_060_3;
         double r_026_2 = Ry*r_016_3 - r_006_3;
         double r_530_2 = Ry*r_520_3 - 2*r_510_3;
         double r_350_2 = Rx*r_250_3 - 2*r_150_3;
         double r_503_2 = Rz*r_502_3 - 2*r_501_3;
         double r_305_2 = Rx*r_205_3 - 2*r_105_3;
         double r_053_2 = Rz*r_052_3 - 2*r_051_3;
         double r_035_2 = Ry*r_025_3 - 2*r_015_3;
         double r_440_2 = Rx*r_340_3 - 3*r_240_3;
         double r_404_2 = Rx*r_304_3 - 3*r_204_3;
         double r_044_2 = Ry*r_034_3 - 3*r_024_3;
         double r_314_2 = Ry*r_304_3;
         double r_341_2 = Rz*r_340_3;
         double r_431_2 = Rz*r_430_3;
         double r_800_1 = Rx*r_700_2 - 7*r_600_2;
         double r_080_1 = Ry*r_070_2 - 7*r_060_2;
         double r_008_1 = Rz*r_007_2 - 7*r_006_2;
         double r_710_1 = Ry*r_700_2;
         double r_170_1 = Rx*r_070_2;
         double r_701_1 = Rz*r_700_2;
         double r_107_1 = Rx*r_007_2;
         double r_071_1 = Rz*r_070_2;
         double r_017_1 = Ry*r_007_2;
         double r_620_1 = Ry*r_610_2 - r_600_2;
         double r_602_1 = Rz*r_601_2 - r_600_2;
         double r_260_1 = Rx*r_160_2 - r_060_2;
         double r_206_1 = Rx*r_106_2 - r_006_2;
         double r_062_1 = Rz*r_061_2 - r_060_2;
         double r_026_1 = Ry*r_016_2 - r_006_2;
         double r_530_1 = Ry*r_520_2 - 2*r_510_2;
         double r_350_1 = Rx*r_250_2 - 2*r_150_2;
         double r_503_1 = Rz*r_502_2 - 2*r_501_2;
         double r_305_1 = Rx*r_205_2 - 2*r_105_2;
         double r_053_1 = Rz*r_052_2 - 2*r_051_2;
         double r_035_1 = Ry*r_025_2 - 2*r_015_2;
         double r_440_1 = Rx*r_340_2 - 3*r_240_2;
         double r_404_1 = Rx*r_304_2 - 3*r_204_2;
         double r_044_1 = Ry*r_034_2 - 3*r_024_2;
         double r_314_1 = Ry*r_304_2;
         double r_341_1 = Rz*r_340_2;
         double r_431_1 = Rz*r_430_2;
         double r_900_1 = Rx*r_800_2 - 8*r_700_2;
         double r_090_1 = Ry*r_080_2 - 8*r_070_2;
         double r_009_1 = Rz*r_008_2 - 8*r_007_2;
         double r_180_1 = Rx*r_080_2;
         double r_108_1 = Rx*r_008_2;
         double r_810_1 = Ry*r_800_2;
         double r_018_1 = Ry*r_008_2;
         double r_801_1 = Rz*r_800_2;
         double r_081_1 = Rz*r_080_2;
         double r_720_1 = Ry*r_710_2 - r_700_2;
         double r_702_1 = Rz*r_701_2 - r_700_2;
         double r_270_1 = Rx*r_170_2 - r_070_2;
         double r_072_1 = Rz*r_071_2 - r_070_2;
         double r_207_1 = Rx*r_107_2 - r_007_2;
         double r_027_1 = Ry*r_017_2 - r_007_2;
         double r_360_1 = Rx*r_260_2 - 2*r_160_2;
         double r_306_1 = Rx*r_206_2 - 2*r_106_2;
         double r_630_1 = Ry*r_620_2 - 2*r_610_2;
         double r_036_1 = Ry*r_026_2 - 2*r_016_2;
         double r_603_1 = Rz*r_602_2 - 2*r_601_2;
         double r_063_1 = Rz*r_062_2 - 2*r_061_2;
         double r_540_1 = Ry*r_530_2 - 3*r_520_2;
         double r_504_1 = Rz*r_503_2 - 3*r_502_2;
         double r_450_1 = Rx*r_350_2 - 3*r_250_2;
         double r_054_1 = Rz*r_053_2 - 3*r_052_2;
         double r_405_1 = Rx*r_305_2 - 3*r_205_2;
         double r_045_1 = Ry*r_035_2 - 3*r_025_2;
         double r_162_1 = Rx*r_062_2;
         double r_126_1 = Rx*r_026_2;
         double r_612_1 = Ry*r_602_2;
         double r_531_1 = Rz*r_530_2;
         double r_513_1 = Ry*r_503_2;
         double r_351_1 = Rz*r_350_2;
         double r_315_1 = Ry*r_305_2;
         double r_153_1 = Rx*r_053_2;
         double r_135_1 = Rx*r_035_2;
         double r_144_1 = Rx*r_044_2;
         double r_414_1 = Ry*r_404_2;
         double r_441_1 = Rz*r_440_2;
         double r_342_1 = Rz*r_341_2 - r_340_2;
         double r_324_1 = Ry*r_314_2 - r_304_2;
         double r_432_1 = Rz*r_431_2 - r_430_2;
         pOut[0] = Rx*r_900_1 - 9*r_800_1;
         pOut[1] = Ry*r_090_1 - 9*r_080_1;
         pOut[2] = Rz*r_009_1 - 9*r_008_1;
         pOut[3] = Ry*r_900_1;
         pOut[4] = Rx*r_090_1;
         pOut[5] = Rz*r_900_1;
         pOut[6] = Rx*r_009_1;
         pOut[7] = Rz*r_090_1;
         pOut[8] = Ry*r_009_1;
         pOut[9] = Ry*r_810_1 - r_800_1;
         pOut[10] = Rz*r_801_1 - r_800_1;
         pOut[11] = Rx*r_180_1 - r_080_1;
         pOut[12] = Rx*r_108_1 - r_008_1;
         pOut[13] = Rz*r_081_1 - r_080_1;
         pOut[14] = Ry*r_018_1 - r_008_1;
         pOut[15] = Ry*r_720_1 - 2*r_710_1;
         pOut[16] = Rx*r_270_1 - 2*r_170_1;
         pOut[17] = Rz*r_702_1 - 2*r_701_1;
         pOut[18] = Rx*r_207_1 - 2*r_107_1;
         pOut[19] = Rz*r_072_1 - 2*r_071_1;
         pOut[20] = Ry*r_027_1 - 2*r_017_1;
         pOut[21] = Ry*r_630_1 - 3*r_620_1;
         pOut[22] = Rz*r_603_1 - 3*r_602_1;
         pOut[23] = Rx*r_360_1 - 3*r_260_1;
         pOut[24] = Rx*r_306_1 - 3*r_206_1;
         pOut[25] = Rz*r_063_1 - 3*r_062_1;
         pOut[26] = Ry*r_036_1 - 3*r_026_1;
         pOut[27] = Ry*r_540_1 - 4*r_530_1;
         pOut[28] = Rz*r_504_1 - 4*r_503_1;
         pOut[29] = Rz*r_054_1 - 4*r_053_1;
         pOut[30] = Ry*r_108_1;
         pOut[31] = Rz*r_180_1;
         pOut[32] = Rz*r_810_1;
         pOut[33] = Ry*r_702_1;
         pOut[34] = Rx*r_072_1;
         pOut[35] = Rz*r_720_1;
         pOut[36] = Rx*r_027_1;
         pOut[37] = Rz*r_270_1;
         pOut[38] = Ry*r_207_1;
         pOut[39] = Ry*r_306_1;
         pOut[40] = Rx*r_036_1;
         pOut[41] = Rz*r_360_1;
         pOut[42] = Rx*r_063_1;
         pOut[43] = Rz*r_630_1;
         pOut[44] = Ry*r_603_1;
         pOut[45] = Ry*r_504_1;
         pOut[46] = Rx*r_054_1;
         pOut[47] = Rz*r_540_1;
         pOut[48] = Rx*r_045_1;
         pOut[49] = Rz*r_450_1;
         pOut[50] = Ry*r_405_1;
         pOut[51] = Ry*r_612_1 - r_602_1;
         pOut[52] = Rx*r_162_1 - r_062_1;
         pOut[53] = Rx*r_126_1 - r_026_1;
         pOut[54] = Rz*r_531_1 - r_530_1;
         pOut[55] = Rz*r_351_1 - r_350_1;
         pOut[56] = Ry*r_513_1 - r_503_1;
         pOut[57] = Ry*r_315_1 - r_305_1;
         pOut[58] = Rx*r_153_1 - r_053_1;
         pOut[59] = Rx*r_135_1 - r_035_1;
         pOut[60] = Rz*r_441_1 - r_440_1;
         pOut[61] = Ry*r_414_1 - r_404_1;
         pOut[62] = Rx*r_144_1 - r_044_1;
         pOut[63] = Ry*r_324_1 - 2*r_314_1;
         pOut[64] = Rz*r_342_1 - 2*r_341_1;
         pOut[65] = Rz*r_432_1 - 2*r_431_1;
         return;
      }
   }
   assert(0);
}

// Form [R]_out := [R + 2ix]_in + [R + 2iy]_in + [R + 2iz]_in: Contract nCartY(lab+2) into nCartY(lab).
void ShellLaplace(double *IR_RP pOut, double const *IR_RP pIn, unsigned LaplaceOrder, unsigned lab)
{
   assert(LaplaceOrder == 1);
   switch(lab) {
      case 0: {
         pOut[0] = pIn[0] + pIn[1] + pIn[2];
         return;
      }
      case 1: {
         pOut[0] = pIn[0] + pIn[3] + pIn[4];
         pOut[1] = pIn[5] + pIn[1] + pIn[6];
         pOut[2] = pIn[7] + pIn[8] + pIn[2];
         return;
      }
      case 2: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[9] + pIn[1] + pIn[11];
         pOut[2] = pIn[10] + pIn[11] + pIn[2];
         pOut[3] = pIn[3] + pIn[4] + pIn[12];
         pOut[4] = pIn[5] + pIn[13] + pIn[6];
         pOut[5] = pIn[14] + pIn[7] + pIn[8];
         return;
      }
      case 3: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[12];
         pOut[2] = pIn[13] + pIn[14] + pIn[2];
         pOut[3] = pIn[9] + pIn[3] + pIn[18];
         pOut[4] = pIn[10] + pIn[18] + pIn[4];
         pOut[5] = pIn[5] + pIn[11] + pIn[19];
         pOut[6] = pIn[19] + pIn[12] + pIn[6];
         pOut[7] = pIn[7] + pIn[20] + pIn[13];
         pOut[8] = pIn[20] + pIn[8] + pIn[14];
         pOut[9] = pIn[15] + pIn[16] + pIn[17];
         return;
      }
      case 4: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[13];
         pOut[2] = pIn[12] + pIn[14] + pIn[2];
         pOut[3] = pIn[3] + pIn[15] + pIn[21];
         pOut[4] = pIn[15] + pIn[4] + pIn[22];
         pOut[5] = pIn[5] + pIn[23] + pIn[16];
         pOut[6] = pIn[16] + pIn[24] + pIn[6];
         pOut[7] = pIn[25] + pIn[7] + pIn[17];
         pOut[8] = pIn[26] + pIn[17] + pIn[8];
         pOut[9] = pIn[9] + pIn[11] + pIn[27];
         pOut[10] = pIn[10] + pIn[27] + pIn[12];
         pOut[11] = pIn[27] + pIn[13] + pIn[14];
         pOut[12] = pIn[21] + pIn[22] + pIn[18];
         pOut[13] = pIn[23] + pIn[19] + pIn[24];
         pOut[14] = pIn[20] + pIn[25] + pIn[26];
         return;
      }
      case 5: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[12];
         pOut[2] = pIn[13] + pIn[14] + pIn[2];
         pOut[3] = pIn[15] + pIn[3] + pIn[24];
         pOut[4] = pIn[16] + pIn[25] + pIn[4];
         pOut[5] = pIn[5] + pIn[17] + pIn[26];
         pOut[6] = pIn[27] + pIn[18] + pIn[6];
         pOut[7] = pIn[7] + pIn[28] + pIn[19];
         pOut[8] = pIn[29] + pIn[8] + pIn[20];
         pOut[9] = pIn[9] + pIn[15] + pIn[33];
         pOut[10] = pIn[10] + pIn[33] + pIn[16];
         pOut[11] = pIn[17] + pIn[11] + pIn[34];
         pOut[12] = pIn[34] + pIn[12] + pIn[18];
         pOut[13] = pIn[19] + pIn[35] + pIn[13];
         pOut[14] = pIn[35] + pIn[20] + pIn[14];
         pOut[15] = pIn[21] + pIn[30] + pIn[31];
         pOut[16] = pIn[30] + pIn[22] + pIn[32];
         pOut[17] = pIn[31] + pIn[32] + pIn[23];
         pOut[18] = pIn[33] + pIn[24] + pIn[25];
         pOut[19] = pIn[26] + pIn[34] + pIn[27];
         pOut[20] = pIn[28] + pIn[29] + pIn[35];
         return;
      }
      case 6: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[13];
         pOut[2] = pIn[12] + pIn[14] + pIn[2];
         pOut[3] = pIn[3] + pIn[15] + pIn[27];
         pOut[4] = pIn[16] + pIn[4] + pIn[28];
         pOut[5] = pIn[5] + pIn[29] + pIn[17];
         pOut[6] = pIn[18] + pIn[30] + pIn[6];
         pOut[7] = pIn[31] + pIn[7] + pIn[19];
         pOut[8] = pIn[32] + pIn[20] + pIn[8];
         pOut[9] = pIn[9] + pIn[21] + pIn[39];
         pOut[10] = pIn[10] + pIn[39] + pIn[22];
         pOut[11] = pIn[21] + pIn[11] + pIn[40];
         pOut[12] = pIn[22] + pIn[41] + pIn[12];
         pOut[13] = pIn[40] + pIn[13] + pIn[23];
         pOut[14] = pIn[41] + pIn[23] + pIn[14];
         pOut[15] = pIn[15] + pIn[16] + pIn[42];
         pOut[16] = pIn[17] + pIn[43] + pIn[18];
         pOut[17] = pIn[44] + pIn[19] + pIn[20];
         pOut[18] = pIn[33] + pIn[34] + pIn[24];
         pOut[19] = pIn[35] + pIn[25] + pIn[36];
         pOut[20] = pIn[26] + pIn[37] + pIn[38];
         pOut[21] = pIn[27] + pIn[42] + pIn[33];
         pOut[22] = pIn[42] + pIn[28] + pIn[34];
         pOut[23] = pIn[29] + pIn[35] + pIn[43];
         pOut[24] = pIn[43] + pIn[36] + pIn[30];
         pOut[25] = pIn[37] + pIn[31] + pIn[44];
         pOut[26] = pIn[38] + pIn[44] + pIn[32];
         pOut[27] = pIn[39] + pIn[40] + pIn[41];
         return;
      }
      case 7: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[12];
         pOut[2] = pIn[13] + pIn[14] + pIn[2];
         pOut[3] = pIn[15] + pIn[3] + pIn[30];
         pOut[4] = pIn[16] + pIn[31] + pIn[4];
         pOut[5] = pIn[5] + pIn[17] + pIn[32];
         pOut[6] = pIn[33] + pIn[18] + pIn[6];
         pOut[7] = pIn[7] + pIn[34] + pIn[19];
         pOut[8] = pIn[35] + pIn[8] + pIn[20];
         pOut[9] = pIn[9] + pIn[21] + pIn[45];
         pOut[10] = pIn[10] + pIn[45] + pIn[22];
         pOut[11] = pIn[23] + pIn[11] + pIn[46];
         pOut[12] = pIn[46] + pIn[12] + pIn[24];
         pOut[13] = pIn[25] + pIn[47] + pIn[13];
         pOut[14] = pIn[47] + pIn[26] + pIn[14];
         pOut[15] = pIn[21] + pIn[15] + pIn[48];
         pOut[16] = pIn[22] + pIn[49] + pIn[16];
         pOut[17] = pIn[17] + pIn[23] + pIn[50];
         pOut[18] = pIn[51] + pIn[24] + pIn[18];
         pOut[19] = pIn[19] + pIn[52] + pIn[25];
         pOut[20] = pIn[53] + pIn[20] + pIn[26];
         pOut[21] = pIn[27] + pIn[36] + pIn[37];
         pOut[22] = pIn[38] + pIn[28] + pIn[40];
         pOut[23] = pIn[39] + pIn[41] + pIn[29];
         pOut[24] = pIn[48] + pIn[30] + pIn[42];
         pOut[25] = pIn[49] + pIn[42] + pIn[31];
         pOut[26] = pIn[32] + pIn[50] + pIn[43];
         pOut[27] = pIn[43] + pIn[51] + pIn[33];
         pOut[28] = pIn[34] + pIn[44] + pIn[52];
         pOut[29] = pIn[44] + pIn[35] + pIn[53];
         pOut[30] = pIn[36] + pIn[38] + pIn[54];
         pOut[31] = pIn[37] + pIn[54] + pIn[39];
         pOut[32] = pIn[54] + pIn[40] + pIn[41];
         pOut[33] = pIn[45] + pIn[48] + pIn[49];
         pOut[34] = pIn[50] + pIn[46] + pIn[51];
         pOut[35] = pIn[52] + pIn[53] + pIn[47];
         return;
      }
      case 8: {
         pOut[0] = pIn[0] + pIn[9] + pIn[10];
         pOut[1] = pIn[11] + pIn[1] + pIn[13];
         pOut[2] = pIn[12] + pIn[14] + pIn[2];
         pOut[3] = pIn[3] + pIn[15] + pIn[33];
         pOut[4] = pIn[16] + pIn[4] + pIn[34];
         pOut[5] = pIn[5] + pIn[35] + pIn[17];
         pOut[6] = pIn[18] + pIn[36] + pIn[6];
         pOut[7] = pIn[37] + pIn[7] + pIn[19];
         pOut[8] = pIn[38] + pIn[20] + pIn[8];
         pOut[9] = pIn[9] + pIn[21] + pIn[51];
         pOut[10] = pIn[10] + pIn[51] + pIn[22];
         pOut[11] = pIn[23] + pIn[11] + pIn[52];
         pOut[12] = pIn[24] + pIn[53] + pIn[12];
         pOut[13] = pIn[52] + pIn[13] + pIn[25];
         pOut[14] = pIn[53] + pIn[26] + pIn[14];
         pOut[15] = pIn[15] + pIn[27] + pIn[54];
         pOut[16] = pIn[27] + pIn[16] + pIn[55];
         pOut[17] = pIn[17] + pIn[56] + pIn[28];
         pOut[18] = pIn[28] + pIn[57] + pIn[18];
         pOut[19] = pIn[58] + pIn[19] + pIn[29];
         pOut[20] = pIn[59] + pIn[29] + pIn[20];
         pOut[21] = pIn[21] + pIn[23] + pIn[60];
         pOut[22] = pIn[22] + pIn[61] + pIn[24];
         pOut[23] = pIn[62] + pIn[25] + pIn[26];
         pOut[24] = pIn[39] + pIn[40] + pIn[30];
         pOut[25] = pIn[41] + pIn[31] + pIn[42];
         pOut[26] = pIn[32] + pIn[43] + pIn[44];
         pOut[27] = pIn[33] + pIn[54] + pIn[45];
         pOut[28] = pIn[55] + pIn[34] + pIn[46];
         pOut[29] = pIn[35] + pIn[47] + pIn[56];
         pOut[30] = pIn[57] + pIn[48] + pIn[36];
         pOut[31] = pIn[49] + pIn[37] + pIn[58];
         pOut[32] = pIn[50] + pIn[59] + pIn[38];
         pOut[33] = pIn[45] + pIn[63] + pIn[39];
         pOut[34] = pIn[63] + pIn[46] + pIn[40];
         pOut[35] = pIn[47] + pIn[41] + pIn[64];
         pOut[36] = pIn[64] + pIn[42] + pIn[48];
         pOut[37] = pIn[43] + pIn[49] + pIn[65];
         pOut[38] = pIn[44] + pIn[65] + pIn[50];
         pOut[39] = pIn[51] + pIn[60] + pIn[61];
         pOut[40] = pIn[60] + pIn[52] + pIn[62];
         pOut[41] = pIn[61] + pIn[62] + pIn[53];
         pOut[42] = pIn[54] + pIn[55] + pIn[63];
         pOut[43] = pIn[56] + pIn[64] + pIn[57];
         pOut[44] = pIn[65] + pIn[58] + pIn[59];
         return;
      }
   }
   assert(0);
}

unsigned char iCartPow[56][3] = {
   {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}, {2,0,0}, {0,2,0}, {0,0,2}, {1,1,0},
   {1,0,1}, {0,1,1}, {3,0,0}, {0,3,0}, {0,0,3}, {1,2,0}, {1,0,2}, {2,1,0},
   {0,1,2}, {2,0,1}, {0,2,1}, {1,1,1}, {4,0,0}, {0,4,0}, {0,0,4}, {3,1,0},
   {1,3,0}, {3,0,1}, {1,0,3}, {0,3,1}, {0,1,3}, {2,2,0}, {2,0,2}, {0,2,2},
   {1,1,2}, {1,2,1}, {2,1,1}, {5,0,0}, {0,5,0}, {0,0,5}, {1,4,0}, {1,0,4},
   {4,1,0}, {0,1,4}, {4,0,1}, {0,4,1}, {3,2,0}, {3,0,2}, {2,3,0}, {0,3,2},
   {2,0,3}, {0,2,3}, {3,1,1}, {1,3,1}, {1,1,3}, {1,2,2}, {2,1,2}, {2,2,1}
}; // 0.16 kb

} // namespace ir
