/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include "atlas_aux.h"
/*
 * Expected defs: BETA[0,1,N,X]
 */
#ifdef BETA1
   #define BETA 1
#elif defined(BETA0)
   #define BETA 0
#elif defined(BETAN)
   #define BETA -1
#else
   #define BETA 2
#endif
#if BETA == 1
   #define BEN _b1
#elif BETA == -1
   #define BEN _bN
#elif BETA == 0
   #define BEN _b0
#else
   #define BEN _bX
#endif
#ifdef Conj_
void Mjoin(Mjoin(PATL,geput1H),BEN)
#else
void Mjoin(Mjoin(PATL,geput1T),BEN)
#endif
   (ATL_CSZT M, ATL_CSZT N, const TYPE *A, ATL_CSZT lda,
    const SCALAR beta, TYPE *C, ATL_CSZT ldc)
/*
 *  C = beta*C + A^T or C = beta*C + A^H, C is MxN  and A is NxM
 *  This algorithm is unblocked, and so should only be called with data that
 *  fits in the cache.  It reads A row-wise, so if lda is large it can have
 *  TLB problems.  A blocked operation (geput) should be made on top of this
 *  one if you need to copy matrices of unbounded size.
 */
{
   ATL_UINT j;
   #ifdef TCPLX
      register TYPE rb=*beta, ib=beta[1];
      const TYPE ONE[2] = {ATL_rone, ATL_rzero};
      ATL_CSZT ldc2 = ldc+ldc, M2=M+M, lda2=lda+lda;
      for (j=0; j < N; j++, C += ldc2)
      {
         const TYPE *a=A+j+j;
         ATL_UINT i;
         for (i=0; i < M2; i += 2, a += lda2)
         {
            #if BETA == 1
               C[i] += *a;
               #ifdef Conj_
                  C[i+1] -= a[1];
               #else
                  C[i+1] += a[1];
               #endif
            #elif BETA == 0
               C[i] = *a;
               #if Conj_
                  C[i+1] = -a[1];
               #else
                  C[i+1] = a[1];
               #endif
            #elif BETA == -1
               C[i] = *a - C[i];
               #ifdef Conj_
                  C[i+1] = -(a[1]+C[i+1]);
               #else
                  C[i+1] = a[1] - C[i+1];
               #endif
            #else
               const register TYPE rc=C[i], ic=C[i+1];
               C[i]   = *a   + rc*rb - ic*ib;
               #ifdef Conj_
                  C[i+1] = rc*ib + ic*rb - a[1];
               #else
                  C[i+1] = a[1] + rc*ib + ic*rb;
               #endif
            #endif
         }
      }
   #else
      for (j=0; j < N; j++, C += ldc)
      {
         const TYPE *a=A+j;
         ATL_UINT i;
         for (i=0; i < M; i++, a += lda)
         #if BETA == 0
            C[i] = *a;
         #elif BETA == 1
            C[i] += *a;
         #elif BETA == -1
            C[i] = *a - C[i];
         #else
            C[i] = beta*C[i] + *a;
         #endif
      }
   #endif
}
