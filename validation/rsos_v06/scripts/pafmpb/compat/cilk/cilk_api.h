/* Cilk runtime calls used by PAFMPB, mapped to OpenMP (see cilk.h). */
#ifndef PAFMPB_CILK_API_COMPAT_H
#define PAFMPB_CILK_API_COMPAT_H
#include <omp.h>
#include <stdlib.h>
#include <string.h>
static inline int __cilkrts_get_nworkers(void) { return omp_get_max_threads(); }
static inline int __cilkrts_set_param(const char *name, const char *value)
{
  if (strcmp(name, "nworkers") != 0) return 1;
  int n = atoi(value);
  if (n < 1) return 1;
  omp_set_num_threads(n);
  return 0;
}
#endif
