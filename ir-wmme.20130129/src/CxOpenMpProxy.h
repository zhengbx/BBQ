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

// include OpenMP header if available or define inline dummy functions
// for OMP primitives if not.
#ifndef OPENMP_PROXY_H
#define OPENMP_PROXY_H

#ifdef _OPENMP
   #include <omp.h>
#else
   inline int omp_get_thread_num() { return 0; } // current thread id
   inline void omp_set_num_threads(int) {};
   inline int omp_get_max_threads() { return 1; } // total number of threads supposed to be running.

   struct omp_lock_t {};
   inline void omp_destroy_lock(omp_lock_t *){};
   inline void omp_init_lock(omp_lock_t *){};
   inline void omp_set_lock(omp_lock_t *){};
#endif

#endif // OPENMP_PROXY_H
