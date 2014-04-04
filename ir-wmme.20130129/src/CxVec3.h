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

#ifndef CX_VEC3_H
#define CX_VEC3_H

#define TVec3 TVector3<T>

namespace ct {

// that is only here to avoid external dependencies...
template<class T>
struct TVector3
{
   typedef T
      value_type;
   T m[3];

   TVector3() {};
   explicit TVector3(T const *xyz) {m[0] = xyz[0]; m[1] = xyz[1]; m[2] = xyz[2];}
   TVector3(T x, T y, T z) {m[0] = x; m[1] = y; m[2] = z;}

   void operator += (TVec3 const &other) {this->m[0] += other.m[0]; this->m[1] += other.m[1]; this->m[2] += other.m[2];}
   void operator -= (TVec3 const &other) {this->m[0] -= other.m[0]; this->m[1] -= other.m[1]; this->m[2] -= other.m[2];}
   void operator *= (T f) {this->m[0] *= f; this->m[1] *= f; this->m[2] *= f;}
   void operator /= (T f) {*this *= 1/f;}

   T &operator[] (unsigned i) { return this->m[i]; }
   T const &operator[] (unsigned i) const { return this->m[i]; }

   operator T* () {return &this->m[0];}
   operator T const* () const {return &this->m[0];}
};

template<class T> TVec3 operator + (TVec3 const &a, TVec3 const &b) { return TVec3(a[0]+b[0], a[1]+b[1], a[2]+b[2]); }
template<class T> TVec3 operator - (TVec3 const &a, TVec3 const &b) { return TVec3(a[0]-b[0], a[1]-b[1], a[2]-b[2]); }
template<class T> TVec3 operator * (T f, TVec3 const &b) { return TVec3(f*b[0], f*b[1], f*b[2]); }
template<class T> TVec3 operator * (TVec3 const &b, T f) { return TVec3(f*b[0], f*b[1], f*b[2]); }
template<class T> T Dot(TVec3 const &a, TVec3 const &b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
template<class T> T LengthSq(TVec3 const &a) { return Dot(a,a); }
template<class T> T DistSq(TVec3 const &a, TVec3 const &b) { return LengthSq(a-b); }

} // namespace ct

#undef TVEC3

#endif // CX_VEC3_H
