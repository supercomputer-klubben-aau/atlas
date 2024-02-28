/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef TREAL

void Mjoin(Mjoin(PATL,sycopyU),NM)
   (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C)
/*
 * Copies a symmetric matrix stored in Upper format to a dense matrix
 */
{
   int i, j;
   const register TYPE alpha=alpha0;
   const TYPE *Ar, *Ac = A;
   TYPE *rC, *cC=C;

   if (N > 1)
   {
      for (j=0; j != N; j++)
      {
         for (i=0; i <= j; i++) C[i] = ATL_MulByALPHA(Ac[i]);
         Ar = Ac + j + lda;
         for (i=j+1; i < N; i++, Ar += lda) C[i] = ATL_MulByALPHA(*Ar);
         C += N;
         Ac += lda;
      }
   }
   else if (N == 1) *C = ATL_MulByALPHA(*A);
}

#else

#ifdef Herm_
void Mjoin(PATL,hecopyU)
#else
void Mjoin(PATL,sycopyU)
#endif
   (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C)
/*
 * Copies a symmetric matrix stored in Upper format to a dense matrix
 */
{
   int i, j2;
   ATL_CSZT  N2=N<<1, lda2=lda<<1;
   const TYPE *a;
   #define ldc2 N2;
   #ifdef Herm_
      const TYPE zero=0.0;
   #endif

   for (j2=0; j2 != N2; j2 += 2)
   {
      a = A;
      for (i=0; i != j2; i++) C[i] = *a++;
      C[j2] = *a;
      #ifdef Herm_
         C[j2+1] = zero;
      #else
         C[j2+1] = a[1];
      #endif
      if (j2 != N2)
      {
         for (i=j2+2; i != N2; i += 2)
         {
            a += lda2;
            C[i] = *a;
            #ifdef Herm_
               C[i+1] = -a[1];
            #else
               C[i+1] = a[1];
            #endif
         }
      }
      A += lda2;
      C += ldc2;
   }
}

#endif
