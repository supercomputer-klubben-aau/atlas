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
void Mjoin(Mjoin(PATL,syr2k_putU),BNM)
   (ATL_CSZT  N, const TYPE *D, const SCALAR beta0, TYPE *A, ATL_CSZT lda)

/*
 * Takes D with property (D + D') = (D + D')', and writes it to
 * upper part of symetric matrix A
 */
{
   register int i, j;
   ATL_CSZT ldap1 = lda+1;
   const TYPE *Dc=D, *Dr;
   TYPE *Ad=A, *Ar;
   const register SCALAR beta=beta0;

   for (j=0; j != N; j++)
   {
      for (Dr=Dc+j, Ar=Ad, i=j; i != N; i++, Dr += N, Ar += lda)
         #if defined(BETA1)
            *Ar += Dc[i] + *Dr;
         #elif defined (BETA0)
            *Ar =  Dc[i] + *Dr;
         #else
            *Ar = *Ar*beta + Dc[i] + *Dr;
         #endif
      Ad += ldap1;
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
   void Mjoin(Mjoin(PATL,her2k_putU),BNM)
#else
   void Mjoin(Mjoin(PATL,syr2k_putU),BNM)
#endif
   (ATL_CSZT  N, const TYPE *D, const SCALAR beta0, TYPE *A, ATL_CSZT lda)

/*
 * Takes D with property (D + D') = (D + D')', and writes it to
 * upper part of symetric matrix A
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
      Dr = D + j2;
      #ifdef BETA1
         for (i=0; i != j2; i += 2, Dr += ldD2)
         {
            A[i] += Dc[i] + *Dr;
            A[i+1] += Dc[i+1] MPi Dr[1];
         }
         A[j2] += Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] += Dc[j2+1] + Dc[j2+1];
         #endif
      #elif defined(BETA0)
         for (i=0; i != j2; i += 2, Dr += ldD2)
         {
            A[i] = Dc[i] + *Dr;
            A[i+1] = Dc[i+1] MPi Dr[1];
         }
         A[j2] = Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] = Dc[j2+1] + Dc[j2+1];
         #endif
      #elif defined(BETAN1) || defined(BETAXI0)
         for (i=0; i != j2; i += 2, Dr += ldD2)
         {
            A[i] = ATL_MulByBETA(A[i]) + Dc[i] + *Dr;
            A[i+1] = ATL_MulByBETA(A[i+1]) + Dc[i+1] MPi Dr[1];
         }
         A[j2] = ATL_MulByBETA(A[j2]) + Dc[j2] + Dc[j2];
         #ifdef Herm_
            A[j2+1] = zero;
         #else
            A[j2+1] = ATL_MulByBETA(A[j2+1]) + Dc[j2+1] + Dc[j2+1];
         #endif
      #else
         for (i=0; i != j2; i += 2, Dr += ldD2)
         {
            ra = A[i]; ia = A[i+1];
            A[i] = ra*rbeta - ia*ibeta + Dc[i] + *Dr;
            A[i+1] = rbeta*ia + ibeta*ra + Dc[i+1] + Dr[1];
         }
         #ifdef Herm_
            A[j2] = A[j2]*rbeta + 2.0*Dc[j2];
            A[j2+1] = ATL_rzero;
         #else
            ra = A[j2]; ia = A[j2+1];
            A[j2] = ra*rbeta - ia*ibeta + 2.0*Dc[j2];
            A[j2+1] = rbeta*ia + ibeta*ra + 2.0*Dc[j2+1];
         #endif
      #endif
      A += lda2;
      Dc += ldD2;
   }
}

#endif
#undef MPi
