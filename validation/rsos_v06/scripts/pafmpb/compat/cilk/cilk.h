/* Cilk Plus -> OpenMP shim for building PAFMPB (AEGB_v2_0) with gcc, which no longer
   supports Cilk Plus.  PAFMPB uses only cilk_for loops, whose bodies are independent
   iterations, so each becomes an OpenMP parallel loop.  Nested cilk_for loops become
   nested parallel regions, which OpenMP runs serially unless nesting is enabled. */
#ifndef PAFMPB_CILK_COMPAT_H
#define PAFMPB_CILK_COMPAT_H
#include <omp.h>
#define cilk_for _Pragma("omp parallel for schedule(guided)") for
#define cilk_spawn
#define cilk_sync
#endif
