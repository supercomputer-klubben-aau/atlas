/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef TREAL

Mjoin(void Mjoin(PATL,sycopyL),NM)
   (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C)
/*
 * Copies a symmetric matrix stored in Lower format to a dense matrix
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
         Ar = A + j;
         for (i=0; i <= j; i++, Ar += lda) C[i] = ATL_MulByALPHA(*Ar);
         for (i=j+1; i < N; i++) C[i] = ATL_MulByALPHA(Ac[i]);
         C += N;
         Ac += lda;
      }
   }
   else if (N == 1) *C = ATL_MulByALPHA(*A);
}

#else

#ifdef Herm_
   void Mjoin(PATL,hecopyL)
#else
   void Mjoin(PATL,sycopyL)
#endif
   (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C)
{
   int j2, i;
   ATL_CSZT  N2 = N<<1, lda2=lda<<1;
   #define ldc2 N2
   const TYPE *a;
   #ifdef Herm_
      const TYPE zero=0.0;
   #endif
   TYPE *c;

   for (j2=0; j2 != N2; j2 += 2)
   {
      a = A + j2;
      for (i=0; i != j2; i += 2)
      {
         C[i] = *a;
         #ifdef Herm_
            C[i+1] = -a[1];
         #else
            C[i+1] = a[1];
         #endif
         a += lda2;
      }
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
            a += 2;
            C[i] = *a;
            C[i+1] = a[1];
         }
      }
      C += ldc2;
   }
}

#endif
