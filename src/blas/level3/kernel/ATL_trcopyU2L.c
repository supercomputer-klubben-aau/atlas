/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef TREAL

#ifdef SDREAL
   #ifdef UnitDiag_
   void Mjoin(ATL_sdtrcopyU2L_U),NM)
   #else
   void Mjoin(PATL_sdtrcopyU2L_N,NM)
   #endif
   (ATL_CSZT  N, const double alpha0, const float *A, ATL_CSZT lda, double *C)
#else
   #ifdef UnitDiag_
   void Mjoin(Mjoin(PATL,trcopyU2L_U),NM)
   #else
   void Mjoin(Mjoin(PATL,trcopyU2L_N),NM)
   #endif
      (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Takes an Upper matrix, transposes it, and copies it to a dense matrix with
 * zeroes above the diagonal (i.e., lower form)
 */
{
   int i, j;
   ATL_CSZT ldap1 = lda+1;
   const register TYPE alpha=alpha0;
   const TYPE *Ac = A, *Ar;
   TYPE *rC, *cC=C;

   if (N > 1)
   {
      for (j=0; j != N; j++)
      {
         for (i=0; i != j; i++) C[i] = 0.0;
         #ifdef UnitDiag_
            C[j] = alpha;
         #else
            C[j] = ATL_MulByALPHA(*Ac);
         #endif
         Ar = Ac + lda;
         for (i=j+1; i < N; i++, Ar += lda) C[i] = ATL_MulByALPHA(*Ar);
         C += N;
         Ac += ldap1;
      }
   }
   else if (N == 1)
   {
      #ifdef UnitDiag_
         *C = alpha;
      #else
         *C = ATL_MulByALPHA(*A);
      #endif
   }
}

#else

#ifdef SDCPLX
#ifdef UnitDiag_
   #ifdef ConjTrans_
      void ATL_cztrcopyU2Lc_U
   #else
      void ATL_cztrcopyU2L_U
   #endif
#else
   #ifdef ConjTrans_
      void ATL_cztrcopyU2Lc_N
   #else
      void ATL_cztrcopyU2L_N
   #endif
#endif
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, TYPE *C)
#else
#ifdef UnitDiag_
   #ifdef ConjTrans_
      void Mjoin(PATL,trcopyU2Lc_U)
   #else
      void Mjoin(PATL,trcopyU2L_U)
   #endif
#else
   #ifdef ConjTrans_
      void Mjoin(PATL,trcopyU2Lc_N)
   #else
      void Mjoin(PATL,trcopyU2L_N)
   #endif
#endif
   (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Takes an Upper matrix, transposes it, and copies it to a dense matrix with
 * zeroes above the diagonal (i.e., lower form)
 */
{
   int i, j2;
   ATL_CSZT lda2=lda<<1, N2=N<<1, ldap1=lda2+2;
   #define ldc N2
   const TYPE *a;
   #ifdef UnitDiag_
      const TYPE one=1.0, zero=0.0;
   #else
      const TYPE zero=0.0;
   #endif

   for (j2=0; j2 != N2; j2 += 2)
   {
      for (i=0; i != j2; i += 2) C[i] = C[i+1] = zero;
      #ifdef UnitDiag_
         C[j2] = one;
         C[j2+1] = zero;
      #else
         C[j2] = *A;
         #ifdef ConjTrans_
            C[j2+1] = -A[1];
         #else
            C[j2+1] = A[1];
         #endif
      #endif
      a = A + lda2;
      if (j2 != N2)
      {
         for (i=j2+2; i != N2; i += 2)
         {
            C[i] = *a;
            #ifdef ConjTrans_
               C[i+1] = -a[1];
            #else
               C[i+1] = a[1];
            #endif
            a += lda2;
         }
      }
      A += ldap1;
      C += ldc;
   }
}

#endif

