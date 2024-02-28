/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
void Mjoin(PATL,geswapT)(ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT lda,
                         TYPE *B, ATL_CINT ldb)
/*
 * This routine swaps with transpose two arrays.  In detail:
 *    C' = A**T
 *    A' = C**T
 *    A is of size MxN, B is of size NxM
 * NOTE: not TLB-optimized, so for use on small matrices.
 * NOTE: A and C cannot overlap.
 */
{
   ATL_INT i;

   for (i=0; i < M; i++)
      Mjoin(PATL,swap)(N, A+(i SHIFT), lda, B+(i SHIFT)*ldb, 1);
}
