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
   #define SYCPY Mjoin(sycpy,UNB_a1)
#elif defined(ALPHAN)
   #define SYCPY Mjoin(sycpy,UNB_an)
#elif defined(ALPHAR)
   #define SYCPY Mjoin(sycpy,UNB_ar)
#else
   #define SYCPY Mjoin(sycpy,UNB_aX)
#endif
void SYCPY(ATL_CSZT N, const SCALAR alpha, const TYPE *A, ATL_CSZT lda,
           TYPE *C, ATL_CSZT ldc)
/*
 * C = alpha *A, A is NxN Upper symmetric/hermitian, C is dense scaled copy
 * Assumes C has ldc ~= N, and ldcxN matrix fits in cache.
 */
{
#ifdef TCPLX
   ATL_CSZT ldc2=ldc+ldc, lda2=lda<<1;
   #if defined(ALPHAR) || defined(ALPHAX)
      const register TYPE ral=(*alpha);
      #ifdef ALPHAX
         const register TYPE ial=alpha[1];
      #endif
   #endif
   ATL_SZT j;
   TYPE *CR = C;
   for (j=0; j < N; j++, C += ldc2, A += lda2, CR += 2)
   {
      TYPE *cr=CR;
      ATL_CSZT NN = j+j;
      ATL_SZT i;
      #if defined(ALPHAR) || defined(ALPHAX)
         register TYPE rA;
      #endif
      register TYPE iA;
      #ifdef Conj_
         #define IAC -iA
      #else
         #define IAC iA
      #endif
      for (i=0; i < NN; i += 2, cr += ldc2)
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
      #ifdef ALPHA1
         C[NN] = A[NN];
         #ifndef Conj_
            C[NN+1] = A[NN+1];
         #endif
      #elif defined(ALPHAN)
         C[NN] = -A[NN];
         #ifndef Conj_
            C[NN+1] = -A[NN+1];
         #endif
      #elif defined(ALPHAR)
         rA = A[NN];
         C[NN] = ral*rA;
         #ifndef Conj_
            C[NN+1] = ral*A[NN+1];
         #endif
      #elif defined(ALPHAX)
         rA = A[NN];
         #ifdef Conj_
            C[NN] = rA * ral;
            C[NN+1] = rA*ial;
         #else
            iA = A[NN+1];
            C[NN] = rA * ral - iA * ial;
            C[NN+1] = rA*ial + iA*ral;
         #endif
      #endif
      #if !defined(ALPHAX) && defined(Conj_)
         C[NN+1] = ATL_rzero;
      #endif
   }
#else
   ATL_SZT j;
   TYPE *CR = C;
   for (j=0; j < N; j++, C += ldc, A += lda)
   {
      ATL_SZT i;
      TYPE *cr = CR++;
      for (i=0; i < j; i++, cr += ldc)
         #ifdef ALPHA1
            C[i] = *cr = A[i];
         #elif defined(ALPHAN)
            C[i] = *cr = -A[i];
         #elif defined(ALPHAX)
            C[i] = *cr = alpha * A[i];
         #endif
      #ifdef ALPHA1
         C[j] = A[j];
      #elif defined(ALPHAN)
         C[j] = -A[j];
      #elif defined(ALPHAX)
         C[j] = A[j] * alpha;
      #endif
   }
#endif
}
