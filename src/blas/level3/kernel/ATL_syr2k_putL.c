/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#include "atlas_kern3.h"
#ifdef USE_AMM
   #undef USE_AMM
#endif

#ifdef Herm_
   #define MPi -
#else
   #define MPi +
#endif
#ifdef TREAL

void Mjoin(Mjoin(PATL,syr2k_putL),BNM)
   (ATL_CSZT  N, const TYPE *D, const SCALAR beta0, TYPE *A, ATL_CSZT lda)

/*
 * Takes D with property (D + D') = (D + D')', and writes it to
 * lower part of symmetric matrix A
 */
{
   register int i, j;
   const TYPE *Dc=D, *Dr;
   TYPE *Ac=A;
   const register SCALAR beta=beta0;

   for (j=0; j != N; j++)
   {
      for (Dr=Dc+j, i=j; i != N; i++, Dr += N)
         #if defined(BETA1)
            Ac[i] += Dc[i] + *Dr;
         #elif defined (BETA0)
            Ac[i] =  Dc[i] + *Dr;
         #else
            Ac[i] = Ac[i]*beta + Dc[i] + *Dr;
         #endif
      Ac += lda;
      Dc += N;
   }
}

#else

/*
 * Workaround for icc errors on IA64Itan2
 */
#ifdef ATL_IntelIccBugs
#pragma optimize("", off)
#endif
#ifdef Herm_
   void Mjoin(Mjoin(PATL,her2k_putL),BNM)
#else
   void Mjoin(Mjoin(PATL,syr2k_putL),BNM)
#endif
   (ATL_CSZT  N, const TYPE *D, const SCALAR beta0, TYPE *A, ATL_CSZT lda)

/*
 * Takes D with property (D + D') = (D + D')', and writes it to
 * lower part of symmetric matrix A
 */
{
   register int i, j2;
   ATL_CSZT  N2=N<<1, lda2=lda<<1;
   #define ldD2 N2
   #ifdef Herm_
      const TYPE zero=0.0;
   #endif
   const TYPE *Dc=D, *Dr;
   #ifdef BETAXI0
      const register TYPE rbeta=*beta0;
   #elif defined(BETAX)
      register TYPE ra, ia;
      const register TYPE rbeta=*beta0, ibeta=beta0[1];
   #endif

   for (j2=0; j2 != N2; j2 += 2)
   {
      Dr = Dc + j2 + ldD2;
      #ifdef BETA1
         A[j2] += Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] += Dc[j2+1] + Dc[j2+1];
         #endif
         if (N2 != j2)
         {
            for (i=j2+2; i != N2; i += 2, Dr += ldD2)
            {
               A[i] += Dc[i] + *Dr;
               A[i+1] += Dc[i+1] MPi Dr[1];
            }
         }
      #elif defined(BETA0)
         A[j2] = Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] = Dc[j2+1] + Dc[j2+1];
         #endif
         if (N2 != j2)
         {
            for (i=j2+2; i != N2; i += 2, Dr += ldD2)
            {
               A[i] = Dc[i] + *Dr;
               A[i+1] = Dc[i+1] MPi Dr[1];
            }
         }
      #elif defined(BETAN1) || defined(BETAXI0)
         A[j2] = ATL_MulByBETA(A[j2]) + Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] = ATL_MulByBETA(A[j2+1]) + Dc[j2+1] + Dc[j2+1];
         #endif
         if (N2 != j2)
         {
            for (i=j2+2; i != N2; i += 2, Dr += ldD2)
            {
               A[i] = ATL_MulByBETA(A[i]) + Dc[i] + *Dr;
               A[i+1] = ATL_MulByBETA(A[i+1]) + Dc[i+1] MPi Dr[1];
            }
         }
      #else
         ra = A[j2];
         ia = A[j2+1];
         A[j2] = ra*rbeta - ia*ibeta + Dc[j2] + Dc[j2];
         A[j2+1] = rbeta*ia + ibeta*ra + Dc[j2+1] + Dc[j2+1];
         if (j2 != N2)
         {
            for (i=j2+2; i != N2; i += 2, Dr += ldD2)
            {
               ra = A[i];
               ia = A[i+1];
               A[i] = ra*rbeta - ia*ibeta + Dc[i] + *Dr;
               A[i+1] = rbeta*ia + ibeta*ra + Dc[i+1] + Dr[1];
            }
         }
      #endif
      A += lda2;
      Dc += ldD2;
   }
}

#endif
#undef MPi
