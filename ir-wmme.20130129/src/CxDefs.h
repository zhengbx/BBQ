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

#ifndef CX_DEFS_H
#define CX_DEFS_H

#define AIC_NO_THROW throw()
#define CX_NO_THROW AIC_NO_THROW

// define macros for restricted pointer extensions. By defining a
// pointer as ``restricted'', we promise the compiler that
// the pointer is not aliased in the current scope

#ifdef __GNUC__ // g++
   #define AIC_RP __restrict__
#elif _MSC_VER // microsoft c++. intel c++ may also understand this syntax.
   #define AIC_RP __restrict
#elif __INTEL_COMPILER
   // compile with -restrict command line option (linux) or -Qrestrict (windows).
   #define AIC_RP restrict
#else
   #define AIC_RP
#endif

#define RESTRICT AIC_RP


#ifdef assert
    #undef assert
#endif
#ifdef assert_rt
   #undef assert_rt
#endif
void AicAssertFail( char const *pExpr, char const *pFile, int iLine );
#define assert_rt(x) if(x) {} else AicAssertFail(#x,__FILE__,__LINE__)

#ifdef _DEBUG
    #define assert(x) assert_rt(x)
#else
    #define assert(x) ((void) 0)
#endif



#endif // CX_DEFS_H
