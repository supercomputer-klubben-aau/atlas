/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2000 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_mv.h"
#include "atlas_r1.h"

int main(int nargs, char **args)
{
   int m, n, lda=4;
   TYPE a[4];

   printf("\nGEMM: NB=%d, lat=%d, mu=%d, nu=%d, ku=%d\n\n", NB, ATL_mmLAT,
          ATL_mmMU, ATL_mmNU, ATL_mmKU);

   ATL_GetPartMVN(a, lda, &m, &n);
   printf("mvN  block = %d x %d; mu=%d, nu=%d\n", m, n, ATL_mvNMU, ATL_mvNNU);
   ATL_GetPartMVT(a, lda, &m, &n);
   printf("mvT  block = %d x %d; mu=%d, nu=%d\n\n", m, n, ATL_mvTMU, ATL_mvTNU);

   ATL_GetPartSYMV(a, lda, &m, &n);
   printf("symv block = %d x %d\n\n", m, n);

   ATL_GetPartR1(a, lda, m, n);
   printf("ger  block = %d x %d; mu=%d, nu=%d\n\n", m, n, ATL_r1MU, ATL_r1NU);

   return(0);
}
