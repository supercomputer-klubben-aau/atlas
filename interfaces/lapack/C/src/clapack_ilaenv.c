/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#define INT
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#include "atlas_lapack.h"
#include "clapack.h"

int clapack_ilaenv(enum ATL_ISPEC ISPEC, enum ATL_LAROUT ROUT,
                   unsigned int OPTS, int N1, int N2, int N3, int N4)
{
   int ATL_ilaenv(enum ATL_ISPEC ISPEC, enum ATL_LAROUT ROUT, unsigned int OPTS,
                  int N1, int N2, int N3, int N4);
   return(ATL_ilaenv(ISPEC, ROUT, OPTS, N1, N2, N3, N4));
}
