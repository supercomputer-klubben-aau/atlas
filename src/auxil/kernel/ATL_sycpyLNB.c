/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef Conj_
   #define sycpy Mjoin(PATL,hecpy)
#else
   #define sycpy Mjoin(PATL,sycpy)
#endif
#ifdef ALPHA1
   #define SYCPY Mjoin(sycpy,LNB_a1)
#elif defined(ALPHAN)
   #define SYCPY Mjoin(sycpy,LNB_an)
#elif defined(ALPHAR)
   #define SYCPY Mjoin(sycpy,LNB_ar)
#else
   #define SYCPY Mjoin(sycpy,LNB_aX)
#endif
void SYCPY(ATL_CSZT N, const SCALAR alpha, const TYPE *A, ATL_CSZT lda,
           TYPE *C, ATL_CSZT ldc)
/*
 * C = alpha *A, A is NxN Lower symetric, C is dense scaled copy
 * Assumes C has ldc ~= N, and ldcxN matrix fits in cache.
 */
{
#ifdef TCPLX
   ATL_CSZT ldc2=ldc+ldc;
   #if defined(ALPHAR) || defined(ALPHAX)
      const register TYPE ral=(*alpha);
      #ifdef ALPHAX
         const register TYPE ial=alpha[1];
      #endif
   #endif
   ATL_SZT j;
   ATL_CSZT incC = (ldc+1)<<1, incA=(lda+1)<<1;
   for (j=0; j < N; j++, C += incC, A += incA)
   {
      TYPE *cr=C+ldc2;
      ATL_CSZT NN = (N-j)<<1;
      ATL_SZT i;
      #ifndef ALPHAX
         register TYPE iA;
      #endif
      #ifdef ALPHA1
         *C = *A;
         #ifndef Conj_
            C[1] = A[1];
         #endif
      #elif defined(ALPHAN)
         *C = -(*A);
         #ifndef Conj_
            C[1] = -A[1];
         #endif
      #elif defined(ALPHAR)
         register TYPE rA = *A;
         *C = ral*rA;
         #ifndef Conj_
            C[1] = ral*A[1];
         #endif
      #elif defined(ALPHAX)
         register TYPE rA = *A, iA;
         #ifdef Conj_
            *C = rA * ral;
            C[1] = rA*ial;
         #else
            iA = A[1];
            *C = rA * ral - iA * ial;
            C[1] = rA*ial + iA*ral;
         #endif
      #endif
      #if !defined(ALPHAX) && defined(Conj_)
         C[1] = ATL_rzero;
      #endif
      #ifdef Conj_
         #define IAC -iA
      #else
         #define IAC iA
      #endif
      for (i=2; i < NN; i += 2, cr += ldc2)
      {
         #ifdef ALPHA1
            C[i] = *cr = A[i];
            C[i+1] = iA = A[i+1];
            cr[1] = IAC;
         #elif defined(ALPHAN)
            C[i] = *cr = -A[i];
            #ifdef Conj_
               cr[1] = iA = A[i+1];
               C[i+1] = -iA;
            #else
               cr[1] = C[i+1] = -A[i+1];
            #endif
         #elif defined(ALPHAR)
            C[i] = *cr = ral*A[i];
            C[i+1] = iA = ral*A[i+1];
            cr[1] = IAC;
         #elif defined(ALPHAX)
            rA = A[i]; iA = A[i+1];
            C[i] = *cr = ral*rA - ial*iA;
            C[i+1] = iA = ral*iA + ial*rA;
            cr[1] = IAC;
         #endif
      }
   }
#else
   ATL_SZT j;
   ATL_CSZT incC = (ldc+1), incA=(lda+1);
   for (j=0; j < N; j++, C += incC, A += incA)
   {
      TYPE *cr=C+ldc;
      ATL_CSZT NN = N-j;
      ATL_SZT i;
      #ifdef ALPHA1
         *C = *A;
      #elif defined(ALPHAN)
         *C = -(*A);
      #elif defined(ALPHAX)
         *C = *A * alpha;
      #endif
      for (i=1; i < NN; i++, cr += ldc)
      #ifdef ALPHA1
         C[i] = *cr = A[i];
      #elif defined(ALPHAN)
         C[i] = *cr = -A[i];
      #elif defined(ALPHAX)
         C[i] = *cr = alpha * A[i];
      #endif
   }
#endif
}
