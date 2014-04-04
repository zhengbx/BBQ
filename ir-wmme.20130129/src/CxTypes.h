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

#ifndef _CX_TYPES_H
#define _CX_TYPES_H

#ifndef _for_each
    #define _for_each(it,con) for ( (it) = (con).begin(); (it) != (con).end(); ++(it) )
#endif

#define __STDC_CONSTANT_MACROS
// ^- ask stdint.h to include fixed-size literal macros (e.g., UINT64_C).
#include <boost/cstdint.hpp>
using boost::uint64_t;
using boost::uint32_t;
using boost::uint16_t;
using boost::uint8_t;
using boost::int64_t;
using boost::int32_t;
using boost::int16_t;
using boost::int8_t;
using std::size_t;
using std::ptrdiff_t;

typedef unsigned int
    uint;
typedef unsigned char
    uchar;
typedef unsigned int
    uint;

#include "CxDefs.h"
#define RESTRICT AIC_RP


namespace ct {
    struct FIntrusivePtrDest;
}

void intrusive_ptr_add_ref( ct::FIntrusivePtrDest const *pExpr );
void intrusive_ptr_release( ct::FIntrusivePtrDest const *pExpr );


namespace ct {
    /// A base class for reference counted objects. Classes derived from this can
    /// be used as target for boost::intrusive_ptr.
    struct FIntrusivePtrDest
    {
        FIntrusivePtrDest() : m_RefCount(0) {};
        inline virtual ~FIntrusivePtrDest() = 0;

        mutable int m_RefCount;
        friend void ::intrusive_ptr_add_ref( FIntrusivePtrDest const *Expr );
        friend void ::intrusive_ptr_release( FIntrusivePtrDest const *Expr );
    };

    inline FIntrusivePtrDest::~FIntrusivePtrDest()
    {
    };
} // namespace ct

inline void intrusive_ptr_add_ref( ct::FIntrusivePtrDest const *pExpr ) {
    pExpr->m_RefCount += 1;
}

inline void intrusive_ptr_release( ct::FIntrusivePtrDest const *pExpr ) {
    assert( pExpr->m_RefCount > 0 );
    pExpr->m_RefCount -= 1;
    if ( pExpr->m_RefCount == 0 )
        delete pExpr;
}

#include <boost/intrusive_ptr.hpp>
// ^- that's for gcc 4.7., which otherwise refuses to instantiate intrusive_ptr,
//    since due to changes in two-phase name lookup, it cannot anymore instantiate *any*
//    templates depending on (non-argument dependent) code which was declared only
//    later in the translation unit.
//    Yes, you heard right: In the case of intrusive_ptr, it means that all
//    reference counting functions must be declared *BEFORE* intrusive_ptr.hpp
//    is #included *ANYWHERE*. Sounds like a fun can of worms. For the
//    beginning: don't include #intrusive_ptr.hpp directly, that might break
//    random code.


#endif // _CX_TYPES_H

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4; indent-width 4; mixedindent off; indent-mode normal;
