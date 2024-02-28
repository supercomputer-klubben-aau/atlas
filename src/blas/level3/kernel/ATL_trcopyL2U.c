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
      void Mjoin(ATL_sdtrcopyL2U_U,NM)
   #else
      void Mjoin(ATL_sdtrcopyL2U_N,NM)
   #endif
   (ATL_CSZT  N, const double alpha0, const float *A, ATL_CSZT lda, double *C)
#else
   #ifdef UnitDiag_
      void Mjoin(Mjoin(PATL,trcopyL2U_U),NM)
   #else
      void Mjoin(Mjoin(PATL,trcopyL2U_N),NM)
   #endif
   (ATL_CSZT  N, const SCALAR alpha0, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Takes a Lower matrix, transposes it, and copies it to a dense matrix with
 * zeros below the diagonal (i.e., upper form)
 */
{
   int i, j;
   const register TYPE alpha=alpha0;
   const TYPE *Ar;

   if (N > 1)
   {
      for (j=0; j != N; j++)
      {
         Ar = A + j;
         for (i=0; i != j; i++, Ar += lda) C[i] = ATL_MulByALPHA(*Ar);
         #ifdef UnitDiag_
            C[j] = alpha;
         #else
            C[j] = ATL_MulByALPHA(*Ar);
         #endif
         for (i=j+1; i < N; i++) C[i] = 0.0;
         C += N;
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
         void ATL_cztrcopyL2Uc_U
      #else
         void ATL_cztrcopyL2U_U
      #endif
   #else
      #ifdef ConjTrans_
         void ATL_cztrcopyL2Uc_N
      #else
         void ATL_cztrcopyL2U_N
      #endif
   #endif
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, double *C)
#else
   #ifdef UnitDiag_
      #ifdef ConjTrans_
         void Mjoin(PATL,trcopyL2Uc_U)
      #else
         void Mjoin(PATL,trcopyL2U_U)
      #endif
   #else
      #ifdef ConjTrans_
         void Mjoin(PATL,trcopyL2Uc_N)
      #else
         void Mjoin(PATL,trcopyL2U_N)
      #endif
   #endif
   (ATL_CSZT  N, const TYPE *A, ATL_CSZT lda, TYPE *C)
#endif
/*
 * Takes a Lower matrix, transposes it, and copies it to a dense matrix with
 * zeros below the diagonal (i.e., upper form)
 */
{
   int i, j2;
   ATL_CSZT  N2=N<<1, lda2=lda<<1;
   #define ldc N2
   const TYPE *a;
   const TYPE zero = 0.0;
   #ifdef UnitDiag_
      const TYPE one = 1.0;
   #endif

   for (j2=0; j2 != N2; j2 += 2)
   {
      a = A;
      for (i=0; i != j2; i += 2)
      {
         C[i] = *a;
         #ifdef ConjTrans_
            C[i+1] = -a[1];
         #else
            C[i+1] = a[1];
         #endif
         a += lda2;
      }
      #ifdef UnitDiag_
         C[j2] = one;
         C[j2+1] = zero;
      #else
         C[i] = *a;
         #ifdef ConjTrans_
            C[i+1] = -a[1];
         #else
            C[i+1] = a[1];
         #endif
      #endif
      if (j2 != N2) for (i=j2+2; i != N2; i += 2) C[i] = C[i+1] = zero;
      A += 2;
      C += ldc;
   }
}

#endif
