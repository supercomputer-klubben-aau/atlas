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
   void Mjoin(Mjoin(ATL_sdtrcopyU2U_U,NM)
   #else
   void Mjoin(Mjoin(ATL_sdtrcopyU2U_N,NM)
   #endif
   (ATL_CSZT  N, const double alpha0, const float *A, ATL_CSZT lda, double *C)
#else
   #ifdef UnitDiag_
   void Mjoin(Mjoin(PATL,trcopyU2U_U),NM)
   #else
   void Mjoin(Mjoin(PATL,trcopyU2U_N),NM)
   #endif
   (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Copies an Upper matrix to a dense matrix with zeros below the diagonal
 */
{
   int i, j;
   const register TYPE alpha=alpha0;
   const TYPE *Ac = A;

   if (N > 1)
   {
      for (j=0; j != N; j++)
      {
         for (i=0; i != j; i++) C[i] = ATL_MulByALPHA(Ac[i]);
         #ifdef UnitDiag_
            C[j] = alpha;
         #else
            C[j] = ATL_MulByALPHA(Ac[j]);
         #endif
         for (i=j+1; i < N; i++) C[i] = 0.0;
         C += N;
         Ac += lda;
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
         void ATL_cztrcopyU2Uc_U
      #else
         void ATL_cztrcopyU2U_U
      #endif
   #else
      #ifdef ConjTrans_
         void ATL_cztrcopyU2Uc_N
      #else
         void ATL_cztrcopyU2U_N
      #endif
   #endif
      (ATL_CSZT  N, const float *A, ATL_CSZT lda, double *C)
#else
   #ifdef UnitDiag_
      #ifdef ConjTrans_
         void Mjoin(PATL,trcopyU2Uc_U)
      #else
         void Mjoin(PATL,trcopyU2U_U)
      #endif
   #else
      #ifdef ConjTrans_
         void Mjoin(PATL,trcopyU2Uc_N)
      #else
         void Mjoin(PATL,trcopyU2U_N)
      #endif
   #endif
      (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Copies an Upper matrix to a dense matrix with zeros below the diagonal
 */
{
   int i, j2;
   ATL_CSZT ldc2=N<<1, lda2=lda<<1;
   #define N2 ldc2
   #ifdef UnitDiag_
      const TYPE one=1.0, zero=0.0;
   #else
      const TYPE zero=0.0;
   #endif

   for (j2=0; j2 != N2; j2 += 2)
   {
      for (i=0; i != j2; i += 2)
      {
         C[i] = A[i];
         #ifdef ConjTrans_
            C[i+1] = -A[i+1];
         #else
            C[i+1] = A[i+1];
         #endif
      }
      #ifdef UnitDiag_
         C[j2] = one;
         C[j2+1] = zero;
      #else
         C[j2] = A[j2];
         #ifdef ConjTrans_
            C[j2+1] = -A[j2+1];
         #else
            C[j2+1] = A[j2+1];
         #endif
      #endif
      if (j2 != N2) for (i=j2+2; i != N2; i += 2) C[i] = C[i+1] = zero;
      C += ldc2;
      A += lda2;
   }
}

#endif
