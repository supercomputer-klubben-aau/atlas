/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"

#if defined(DCPLX)
   #define iamax ATL_izamax
#elif defined(SCPLX)
   #define iamax ATL_icamax
#elif defined(SREAL)
   #define iamax ATL_isamax
#else
   #define iamax ATL_idamax
#endif
TYPE Mjoin(PATL,gemaxnrm)
   (ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT lda)
/*
 * Returns maximum absolute value in A;
 * for complex value CV returns maximum of abs(real(CV)) + abs(cplx(CV))
 */
{
   TYPE maxval=0.0, mv;
   ATL_INT j, i;
   #ifdef TCPLX
      ATL_CINT lda2 = lda+lda;
   #else
      #define lda2 lda
   #endif

   for (j=0; j < N; j++, A += lda2)
   {
      i = iamax(M, A, 1);
      #ifdef TCPLX
         i += i;
         mv = Mabs(A[i]) + Mabs(A[i+1]);
      #else
         mv = Mabs(A[i]);
      #endif
      maxval = (maxval >= mv) ? maxval : mv;
   }
   return(maxval);
}
#ifndef TCPLX
   #undef lda2
#endif
