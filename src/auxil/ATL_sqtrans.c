/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
#define NB 32
#define ATL_MulByNB(n_) ((n_)<<5)
#define ATL_DivByNB(n_) ((n_)>>5)

static void Mjoin(PATL,sqtrans0)(ATL_CINT N, TYPE *C, ATL_CINT ldc)
/*
 * Does an in-place transpose of a square matrix.
 * NOTE: this should only be used on small matrices, as it is not optimized
 *       for the TLB.
 */
{
   ATL_INT j;
/*
 * We will work by reflecting swapping columns & rows across diagonal,
 * starting from the last column, so that early cols are retained in cache
 */
   for (j=N-1; j; j--)
      Mjoin(PATL,swap)(j, C+((size_t)ldc)*(j SHIFT), 1, C+(j SHIFT), ldc);
}

void Mjoin(PATL,sqtrans)(ATL_CINT N, TYPE *C, ATL_CINT ldc)
/*
 * Does an in-place transpose of a square matrix.  This routine is blocked
 * to help with TLB
 */
{
   const size_t ldt = ldc;
   ATL_CINT Nnb = ATL_MulByNB(ATL_DivByNB(N)), Nr = N - Nnb;
   ATL_INT i, j;

   if (N < NB+NB)
   {
      Mjoin(PATL,sqtrans0)(N, C, ldc);
      return;
   }
/*
 * Loop in reverse order, so first part of matrix retained in cache
 */
   if (Nr)
   {
      for (i=0; i < Nnb; i += NB)
         Mjoin(PATL,geswapT)(NB, Nr, C+((Nnb*ldt+i)SHIFT), ldc,
                             C+((Nnb+i*ldt)SHIFT), ldc);
      Mjoin(PATL,sqtrans0)(Nr, C+((Nnb*(ldt+1))SHIFT), ldc);
   }
   for (j=Nnb-NB; j >= 0; j -= NB)
   {

      for (i=0; i < j; i += NB)
         Mjoin(PATL,geswapT)(NB, NB, C+((j*ldt+i)SHIFT), ldc,
                             C+((j+i*ldt)SHIFT), ldc);
      Mjoin(PATL,sqtrans0)(NB, C+((j*(ldt+1))SHIFT), ldc);
   }
}
