/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_level2.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,syrk_view.h)) /* blocking & perf info */
#include Mstr(Mjoin(ATLAS_PRE,amm_sqsyrk.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_umsyrk.h))

static INLINE void ATL_syr1
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT N,
#ifdef Conj_
   const TYPE ralpha,
#else
   const SCALAR alpha,
#endif
   const TYPE *A,
   ATL_CSZT lda,
#ifdef Conj_
   const TYPE rbeta,
#else
   const SCALAR beta,
#endif
   TYPE *C,
   ATL_CSZT ldc
)
{
/*
 * First, see if we can simply call the Level-2 BLAS
 */
   #ifdef Conj_
      TYPE *c=C+1;
   #endif
   #ifndef TCPLX
      if (beta == 1.0)
      {
         Mjoin(PATL,syr)(Uplo, N, alpha, A, TA==AtlasNoTrans?1:lda, C, ldc);
         return;
      }
   #endif
   {
      size_t k;
      TYPE *X=(TYPE*)A;
      void *vp=NULL;
      #ifdef Conj_
         const TYPE beta[2] = {rbeta, ATL_rzero};
      #endif
/*
 *    Copy A if it's a row so it is contiguous for all N axpby calls
 */
      #ifdef Conj_
      if (TA == AtlasConjTrans)
      #else
      if (TA == AtlasTrans || TA == AtlasConjTrans)
      #endif
      {
         vp = malloc(ATL_MulBySize(N)+ATL_Cachelen);
         ATL_assert(vp);
         X = ATL_AlignPtr(vp);
         #ifdef Conj_
            Mjoin(PATL,copyConj)(N, A, lda, X, 1);
         #else
            Mjoin(PATL,copy)(N, A, lda, X, 1);
         #endif
      }
      #ifdef Conj_
         if (rbeta == 1.0)
         {
            Mjoin(PATL,her)(Uplo, N, ralpha, X, 1, C, ldc);
            if (vp)
               free(vp);
            return;
         }
      #endif  /* NO syr for complex, so cplx SYR has to use axpy-based code */
      #ifdef TCPLX
         if (Uplo == AtlasUpper)
         {
            size_t ldc2 = ldc+ldc;
            #ifndef Conj_
               const register TYPE ra=(*alpha), ia=alpha[1];
            #endif
            for (k=0; k < N; k++, C += ldc2)
            {
               size_t k2=k+k;
               #ifndef Conj_
                  const register TYPE rx=X[k2], ix=X[k2+1];
               #endif
               TYPE scal[2];
               #ifdef Conj_
                  scal[0] = ralpha*X[k2];
                  scal[1] = -ralpha*X[k2+1];
               #else
                  scal[0] = ra*rx - ia*ix;
                  scal[1] = ra*ix + ia*rx;
               #endif
               Mjoin(PATL,axpby)(k+1, scal, X, 1, beta, C, 1);
            }
         }
         else   /* Uplo == AtlasLower */
         {
            size_t ldcp1 = (ldc+1)SHIFT;
            #ifndef Conj_
               const register TYPE ra=(*alpha), ia=alpha[1];
            #endif
            for (k=0; k < N; k++, C += ldcp1, X += 2)
            {
               #ifndef Conj_
                  const register TYPE rx=*X, ix=X[1];
               #endif
               TYPE scal[2];
               #ifdef Conj_
                  scal[0] = ralpha* *X;
                  scal[1] = ralpha * (-X[1]);
               #else
                  scal[0] = ra*rx - ia*ix;
                  scal[1] = ra*ix + ia*rx;
               #endif
               Mjoin(PATL,axpby)(N-k, scal, X, 1, beta, C, 1);
            }
         }
      #else  /* real with beta != 1.0 */
         if (Uplo == AtlasUpper)
         {
            for (k=0; k < N; k++, C += ldc)
               Mjoin(PATL,axpby)(k+1, X[k]*alpha,  X, 1, beta, C, 1);
         }
         else   /* Uplo == AtlasLower */
         {
            for (k=0; k < N; k++, C += (ldc+1), X++)
               Mjoin(PATL,axpby)(N-k, (*X)*alpha,  X, 1, beta, C, 1);
         }
      #endif
      if (vp)
         free(vp);
   }
   #ifdef Conj_  /* must zero complex part of diagonal for HERK! */
      Mjoin(PATLU,zero)(N, c, (ldc+1)SHIFT);
   #endif
}

#include "atlas_cache.h"
#include "atlas_reflevel3.h"
/*
 * This routine assumes 1st col of A is x, 2nd y, and NoTrans case
 */
static void ATL_syr2_axpy
(
   const enum ATLAS_UPLO Uplo,
   ATL_CSZT N,
#ifdef Conj_
   const TYPE ralpha,
#else
   const SCALAR alpha,
#endif
   const TYPE *x,
   const TYPE *y,
#ifdef Conj_
   const TYPE rbeta,
#else
   const SCALAR beta,
#endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   #ifdef Conj_
      const TYPE beta[2] ={rbeta,  ATL_rzero};
   #endif
   if (Uplo == AtlasLower)
   {
      ATL_CSZT ldcp1 = (ldc+1)SHIFT;
      ATL_SZT j;
      for (j=0; j < N; j++, C += ldcp1)
      {
         ATL_CSZT n = N-j;
         #ifdef Conj_
            const TYPE alpx[2] = {ralpha*x[j+j], -ralpha*x[j+j+1]};
            const TYPE alpy[2] = {ralpha*y[j+j], -ralpha*y[j+j+1]};
         #elif defined(TCPLX)
            const TYPE ralp=alpha[0], ialp=alpha[1];
            TYPE alpx[2], alpy[2], rx=x[j+j], ix=x[j+j+1];
            if (ialp == ATL_rzero)
            {
               alpx[0] = ralp*rx;
               alpx[1] = ralp*ix;
               rx=y[j+j], ix=y[j+j+1];
               alpy[0] = ralp*rx;
               alpy[1] = ralp*ix;
            }
            else
            {
               alpx[0] = ralp*rx - ialp*ix;
               alpx[1] = ralp*ix + ialp*rx;
               rx=y[j+j], ix=y[j+j+1];
               alpy[0] = ralp*rx - ialp*ix;
               alpy[1] = ralp*ix + ialp*rx;
            }
         #else
            const TYPE alpx = alpha * x[j], alpy = alpha * y[j];
         #endif
         Mjoin(PATL,axpby)(n, alpx, x+(j SHIFT), 1, beta, C, 1);
         Mjoin(PATL,axpy)(n, alpy, y+(j SHIFT), 1, C, 1);
         #ifdef Conj_
            C[1] = ATL_rzero;
         #endif
      }
   }
   else /* Uplo == AtlasUpper */
   {
      ATL_CSZT ldc2 = ldc SHIFT;
      ATL_SZT j;
      for (j=0; j < N; j++, C += ldc2)
      {
         #ifdef Conj_
            const TYPE alpx[2] = {ralpha*x[j+j], -ralpha*x[j+j+1]};
            const TYPE alpy[2] = {ralpha*y[j+j], -ralpha*y[j+j+1]};
         #elif defined(TCPLX)
            const TYPE ralp=alpha[0], ialp=alpha[1];
            TYPE alpx[2], alpy[2], rx=x[j+j], ix=x[j+j+1];
            if (ialp == ATL_rzero)
            {
               alpx[0] = ralp*rx;
               alpx[1] = ralp*ix;
               rx=y[j+j], ix=y[j+j+1];
               alpy[0] = ralp*rx;
               alpy[1] = ralp*ix;
            }
            else
            {
               alpx[0] = ralp*rx - ialp*ix;
               alpx[1] = ralp*ix + ialp*rx;
               rx=y[j+j], ix=y[j+j+1];
               alpy[0] = ralp*rx - ialp*ix;
               alpy[1] = ralp*ix + ialp*rx;
            }
         #else
            const TYPE alpx = alpha * x[j], alpy = alpha * y[j];
         #endif
         Mjoin(PATL,axpby)(j+1, alpx, x, 1, beta, C, 1);
         Mjoin(PATL,axpy)(j+1, alpy, y, 1, C, 1);
         #ifdef Conj_
            C[j+j+1] = ATL_rzero;
         #endif
      }
   }
}

#ifdef Conj_
   #define ATL_ger2 Mjoin(PATL,ger2c)
#elif defined(TCPLX)
   #define ATL_ger2 Mjoin(PATL,ger2u)
#else
   #define ATL_ger2 Mjoin(PATL,ger2)
#endif
static void ATL_syr2
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT N,
#ifdef Conj_
   const TYPE ralpha,
#else
   const SCALAR alpha,
#endif
   const TYPE *A,
   ATL_CSZT lda,
#ifdef Conj_
   const TYPE rbeta,
#else
   const SCALAR beta,
#endif
   TYPE *C,
   ATL_CSZT ldc
)
{
   void *vp=NULL;
   TYPE *x, *y;
   #ifdef TCPLX
      ATL_CSZT lda2 = lda+lda, ldc2 = ldc+ldc;
   #else
      #define lda2 lda
      #define ldc2 ldc
   #endif
   #ifdef Conj_
      const TYPE alpha[2]={ralpha, ATL_rzero};
      const TYPE beta[2] ={rbeta,  ATL_rzero};
   #endif
/*
 * If A cols are strided, copy them to contiguous & aligned storage
 */
   if (TA != AtlasNoTrans)
   {
      vp = malloc(ATL_MulBySize(N+N) + ATL_Cachelen+ATL_Cachelen);
      ATL_assert(vp);
      x = ATL_AlignPtr(vp);
      y = x + (N SHIFT);
      y = ATL_AlignPtr(y);
      #ifdef Conj_
         Mjoin(PATL,copyConj)(N, A, lda, x, 1);
         Mjoin(PATL,copyConj)(N, A+(1 SHIFT), lda, y, 1);
      #else
         Mjoin(PATL,copy)(N, A, lda, x, 1);
         Mjoin(PATL,copy)(N, A+(1 SHIFT), lda, y, 1);
      #endif
   }
   else
   {
      x = (TYPE *) A;
      y = (TYPE *) (A + lda2);
   }
/*
 * If BETA != 1, then we must pass through C twice, so use daxpy-based code
 */
   #ifdef Conj_
   if (rbeta != ATL_rone)
   #else
   if (!SCALAR_IS_ONE(beta))
   #endif
   {
      #ifdef Conj_
         ATL_syr2_axpy(Uplo, N, ralpha, x, y, rbeta, C, ldc);
      #else
         ATL_syr2_axpy(Uplo, N, alpha, x, y, beta, C, ldc);
      #endif
   }
   else /* beta=1 can use ger2 based code, make only 1 pass thru C */
   {
      #ifdef Conj_
         const TYPE alpha[2] = {ralpha, ATL_rzero};
      #endif
      #if L1C_ELTS >= 16384
         unsigned int Np = (N > 128) ? 120 : N;
      #elif L1C_ELTS >= 8192
         unsigned int Np = (N > 90) ? 80 : N;
      #elif L1C_ELTS >= 4096
         unsigned int Np = (N > 64) ? 60 : N;
      #elif L1C_ELTS >= 2048
         unsigned int Np = (N > 48) ? 40 : N;
      #else
         unsigned int Np = (N > 40) ? 24 : N;
      #endif
      if (Uplo == AtlasLower)
      {
         ATL_CSZT incC=((ldc+1)SHIFT)*Np, incX=(Np SHIFT);
         ATL_SZT j;
         for (j=0; j < N; j += Np, C += incC, x += incX, y += incX)
         {
            const unsigned nb = Mmin(N-j, Np), nb2=nb SHIFT;
            #ifdef Conj_
               ATL_syr2_axpy(Uplo, nb, ralpha, x, y, rbeta, C, ldc);
            #else
               ATL_syr2_axpy(Uplo, nb, alpha, x, y, beta, C, ldc);
            #endif
            ATL_ger2(N-j-nb, nb, alpha, x+nb2, 1, x, 1, alpha, y+nb2, 1, y, 1,
                     C+nb2, ldc);
         }
      }
      else
      {
         const TYPE *x0 = x, *y0 = y;
         ATL_CSZT incC=((ldc)SHIFT)*Np, incX=(Np SHIFT);
         ATL_SZT j;
         for (j=0; j < N; j += Np, C += incC, x += incX, y += incX)
         {
            const unsigned nb = Mmin(N-j, Np), nb2=nb SHIFT;
            if (j)
               ATL_ger2(j, nb, alpha, x0, 1, x, 1, alpha, y0, 1, y, 1,
                        C, ldc);
            #ifdef Conj_
               ATL_syr2_axpy(Uplo, nb, ralpha, x, y, rbeta, C+(j SHIFT), ldc);
            #else
               ATL_syr2_axpy(Uplo, nb, alpha, x, y, beta, C+(j SHIFT), ldc);
            #endif
         }
      }
   }
   if (vp)
      free(vp);
}

static INLINE int N2Idx(const unsigned int N)
/*
 * Uses Recursive halving to find smallest NB >= N in VWsyrk.
 */
{
   unsigned int idxB=ATL_VWsyrk_NCASES-1, idxS=0, idxM;
   unsigned int nbB, nbS=ATL_VWsyrk_MIN_NB, nbM;
   if (N <= ATL_VWsyrk_MIN_NB)
      return(0);
   nbB = ATL_GetVWsyrkMB(idxB);
   if (N >= nbB)
      return(ATL_VWsyrk_NCASES-1);
   KEEP_ON:
      idxM = ((idxB-idxS)>>1)+idxS;
      if (idxM == idxS)
         return((nbS >= N) ? idxS:idxB);
      nbM = ATL_GetVWsyrkMB(idxM);
      if (nbM > N)
      {
         idxB = idxM;
         nbB = nbM;
      }
      else if (nbM < N)
      {
         idxS = idxM;
         nbS = nbM;
      }
      else /* if (nbM == N) */
         return(idxM);
   goto KEEP_ON;
}

#ifdef Conj_
void Mjoin(PATL,ipherk)
#else
void Mjoin(PATL,ipsyrk)
#endif
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CSZT N,
   ATL_CSZT K,
#ifdef Conj_
   const TYPE ralpha,
#else
   const SCALAR alpha,
#endif
   const TYPE *A,
   ATL_CSZT lda,
#ifdef Conj_
   const TYPE rbeta,
#else
   const SCALAR beta,
#endif
   TYPE *C,
   ATL_CSZT ldc
)
/*
 * C NxN, A NxK
 * SYRK:
 *    C = alpha * A * A^T + beta*C, if TA == AtlasNoTrans
 *    C = alpha * A^T * A + beta*C, if TA == AtlasTrans
 * HERK:
 *    C = alpha * A * A^H + beta*C, if TA == AtlasNoTrans
 *    C = alpha * A^H * A + beta*C, if TA == AtlasTrans
 */
{
   cm2am_t a2blk;
   ablk2cmat_t blk2c;
   void *vp;
   #ifdef TCPLX
      TYPE *rA, *iA, *rC, *iC, *c;
      TYPE one[2] = {ATL_rone, ATL_rzero};
      #ifdef Conj_
         TYPE *crA, *ciA;
         const TYPE alpha[2] = {ralpha, ATL_rzero};
         const TYPE beta[2] = {rbeta, ATL_rzero};
      #endif
   #else
      TYPE *pA, *pC, *c;
      #define one ATL_rone
   #endif
   size_t szA, szC, szE, incAk, nkb, k;
   ATL_CUINT nnu = (N+ATL_SYRKK_NU-1)/ATL_SYRKK_NU, NN=nnu*ATL_SYRKK_NU;
   ATL_UINT kb, kbS, kb0, KB0, idx;
/*
 * Handle degenerate cases
 */
   if (!N)                            /* no output! */
      return;
   if (SCALAR_IS_ZERO(alpha) || !K)  /* really scale of C */
   {
      if (SCALAR_IS_ONE(beta))       /* no-op */
         return;
      if (SCALAR_IS_ZERO(beta))      /* explicit zero */
      {
         if (Uplo == AtlasLower)
            Mjoin(PATL,trsetL)(N, N, beta, beta, C, ldc);
         else
            Mjoin(PATL,trsetU)(N, N, beta, beta, C, ldc);
         return;
      }
      Mjoin(PATL,trscal)(Uplo, N, N, beta, C, ldc);
      #ifdef Conj_  /* must zero complex part of diagonal for HERK! */
         Mjoin(PATLU,zero)(N, C+1, (ldc+1)SHIFT);
      #endif
      return;
   }
   if (K == 1)  /* Level-2/1 BLAS */
   {
      #ifdef Conj_
         ATL_syr1(Uplo, TA, N, ralpha, A, lda, rbeta, C, ldc);
      #else
         ATL_syr1(Uplo, TA, N, alpha, A, lda, beta, C, ldc);
      #endif
      return;
   }
   if (K == 2)  /* Level-2/1 BLAS */
   {
      #ifdef Conj_
         ATL_syr2(Uplo, TA, N, ralpha, A, lda, rbeta, C, ldc);
      #else
         ATL_syr2(Uplo, TA, N, alpha, A, lda, beta, C, ldc);
      #endif
      return;
   }
   if (N == 1)  /* dot product */
   {
      const size_t incA = (TA==AtlasNoTrans) ? lda : 1;
      #ifdef TCPLX
         TYPE dot[2];
         #ifdef Conj_
            Mjoin(PATL,dotc_sub)(K, A, incA, A, incA, dot);
            *dot   *= ralpha;
            if (rbeta != ATL_rzero)
               *dot   += rbeta * *C;
            *C = *dot;
            C[1] = ATL_rzero;
         #else
            const register TYPE ra=(*alpha), ia=alpha[1];
            register TYPE rd, id, rr;
            Mjoin(PATL,dotu_sub)(K, A, incA, A, incA, dot);
            rr = rd = dot[0];
            id = dot[1];
            rd = rr*ra - id*ia;
            id = rr*ia + id*ra;
            if (!SCALAR_IS_ZERO(beta))
            {
               const register TYPE rb=(*beta), ib=beta[1];
               const register TYPE rc=(*C), ic=C[1];
               rd += rb*rc - ib*ic;
               id += rb*ic + ib*rc;
            }
            C[0] = rd;
            C[1] = id;
         #endif
      #else
         TYPE dot;
         dot = Mjoin(PATL,dot)(K, A, incA, A, incA);
         dot *= alpha;
         if (beta != ATL_rzero)
            dot += beta * *C;
         *C = dot;
      #endif
      return;
   }
/*
 * Find NB closest to our present N, and set initial KB to its tuned version
 */
   idx = N2Idx(NN);
   kb = ATL_GetVWsyrkKB(idx);
   szC = ((ATL_SYRKK_NU*ATL_SYRKK_NU+ATL_SYRKK_VLEN-1)/ATL_SYRKK_VLEN)
         * ATL_SYRKK_VLEN;
   szC *= ((nnu+1)*nnu)>>1;  /* only need lower tri blks, not full nnu*nnu */
/*
 * Our SYRK C copy is always to Lower, so if output is Upper, will need
 * workspace to put the Lower part, before reflecting it to Upper.
 * Need max of this extra space, or preload distance.
 */
   szE = (ATL_SYRKK_NU*ATL_SYRKK_NU)<<1;
   if (K > kb)
   {
      nkb = K / kb;
      kb0 = K - nkb*kb;
      if (!kb0)
      {
         kb0 = kb;
         nkb--;
      }
      KB0 = ((kb0+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
      kbS = ((kb+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   }
   else
   {
      nkb = 0;
      kb0 = K;
      kbS = kb = KB0 = ((kb0+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   }
   if (IS_COLMAJ(TA))
   {
      incAk = lda*(kb SHIFT);
      a2blk = Mjoin(PATL,a2blk_syrkT);
   }
   else
   {
      incAk = kb SHIFT;
      a2blk = Mjoin(PATL,a2blk_syrkN);
   }
   szA = nnu*ATL_SYRKK_NU * kbS;
   vp = malloc(ATL_MulBySize(szA + szC + szE)
               + 2*ATL_Cachelen);
   ATL_assert(vp);
   #ifdef TCPLX
      iA = ATL_AlignPtr(vp);
      rA = iA + szA;
      iC = rA + szA;
      iC = ATL_AlignPtr(iC);
      rC = iC + szC;
      #ifdef Conj_
         crA = (TA == AtlasNoTrans) ? rA : iA;
         ciA = (TA == AtlasNoTrans) ? iA : rA;
      #endif
      c = rC + szC;
   #else
      pA = ATL_AlignPtr(vp);
      pC = pA + (szA SHIFT);
      pC = ATL_AlignPtr(pC);
      c = pC + szC;
   #endif
#if 0
   ipinfo_t ip;
   int i, flg;
   Mjoin(PATL,ipgenInfo)(&ip, 0, TA, TA, N, N, K, lda, lda, ldc, alpha, beta);
   flg = (Uplo == AtlasLower) ? 0 : 1;
   if (TA == AtlasNoTrans)
      flg |= 2;
   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb1);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1b1):Mjoin(PATL,SyrkIntoC_aXb1);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbN);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1bN):Mjoin(PATL,SyrkIntoC_aXbN);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb0);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1b0):Mjoin(PATL,SyrkIntoC_aXb0);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbX);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1bX):Mjoin(PATL,SyrkIntoC_aXbX);
      }
   }
   else
   {
      if (SCALAR_IS_NONE(alpha))
         blk2c = Mjoin(PATL,SyrkIntoC_aNb0);
      else
         blk2c = SCALAR_IS_ONE(alpha) ?
                 Mjoin(PATL,SyrkIntoC_a1b0):Mjoin(PATL,SyrkIntoC_aXb0);
   }
   Mjoin(PATL,syrkBlk)(&ip, flg|4, 0, 0, A, a2blk, K<=ip.kb?blk2c:NULL, beta, C,
                       pA, NULL, NULL, NULL, NULL, pC, pA);
/*   A += (kb0 SHIFT) * (IS_COLMAJ(TA) ? lda : 1); */
   nkb = K / ip.kb;
   if (nkb * ip.kb == K)
      nkb--;
   for (k=0; k < nkb; k++)
      Mjoin(PATL,syrkBlk)(&ip, flg, 0, k+1, A, a2blk, k == nkb-1 ? blk2c:NULL,
                          beta, C, pA, NULL, NULL, NULL, NULL, pC, pA);
#else

   #ifdef TCPLX
      a2blk(kb0, N, one, A, lda, rA, iA);
      #ifdef Conj_
         Mjoin(PATL,amsyrkK_b0)(nnu, nnu, KB0, iA, iA, rC, crA, ciA, iC);
         Mjoin(PATL,amsyrkK_b0)(nnu, nnu, KB0, crA, ciA, iC, rA, rA, rC);
         Mjoin(PATL,amsyrkK_b1)(nnu, nnu, KB0, rA, rA, rC, ciA, crA, iC);
         Mjoin(PATL,amsyrkK_bn)(nnu, nnu, KB0, ciA, crA, iC, iA, iA, rC);
      #else
         Mjoin(PATL,amsyrkK_b0)(nnu, nnu, KB0, iA, iA, rC, rA, iA, iC);
         Mjoin(PATL,amsyrkK_b0)(nnu, nnu, KB0, rA, iA, iC, rA, rA, rC);
         Mjoin(PATL,amsyrkK_bn)(nnu, nnu, KB0, rA, rA, rC, iA, rA, iC);
         Mjoin(PATL,amsyrkK_b1)(nnu, nnu, KB0, iA, rA, iC, iA, iA, rC);
      #endif
   #else
      a2blk(kb0, N, ATL_rone, A, lda, pA);
      Mjoin(PATL,amsyrkK_b0)(nnu, nnu, KB0, pA, pA, pC, pA, pA, pC);
   #endif
   A += (kb0 SHIFT) * (IS_COLMAJ(TA) ? lda : 1);
   for (k=0; k < nkb; k++, A += incAk)
   {
      #ifdef TCPLX
         a2blk(kb, N, one, A, lda, rA, iA);
         #ifdef Conj_
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, iA, iA, rC, crA, ciA, iC);
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, crA, ciA, iC, rA, rA, rC);
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, rA, rA, rC, ciA, crA, iC);
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, ciA, crA, iC, iA, iA, rC);
         #else
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, iA, iA, rC, rA, iA, iC);
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, rA, iA, iC, rA, rA, rC);
            Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, rA, rA, rC, iA, rA, iC);
            Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, iA, rA, iC, iA, iA, rC);
         #endif
      #else
         a2blk(kb, N, ATL_rone, A, lda, pA);
         Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, pA, pA, pC, pA, pA, pC);
      #endif
   }
   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb1);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1b1):Mjoin(PATL,SyrkIntoC_aXb1);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbN);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1bN):Mjoin(PATL,SyrkIntoC_aXbN);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb0);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1b0):Mjoin(PATL,SyrkIntoC_aXb0);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbX);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                    Mjoin(PATL,SyrkIntoC_a1bX):Mjoin(PATL,SyrkIntoC_aXbX);
      }
      #ifdef TCPLX
         blk2c(N, N, alpha, rC, iC, beta, C, ldc);
         #ifdef Conj_  /* must zero complex part of diagonal! */
            Mjoin(PATLU,zero)(N, C+1, (ldc+1)SHIFT);
         #endif
      #else
         blk2c(N, N, alpha, pC, beta, C, ldc);
      #endif
   }
   else /* Upper */
   {
      if (SCALAR_IS_ONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb1_L2UT);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b1_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb1_L2UT);
      }
      else if (SCALAR_IS_NONE(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbN_L2UT);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bN_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbN_L2UT);
      }
      else if (SCALAR_IS_ZERO(beta))
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNb0_L2UT);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1b0_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXb0_L2UT);
      }
      else
      {
         if (SCALAR_IS_NONE(alpha))
            blk2c = Mjoin(PATL,SyrkIntoC_aNbX_L2UT);
         else
            blk2c = SCALAR_IS_ONE(alpha) ?
                     Mjoin(PATL,SyrkIntoC_a1bX_L2UT)
                     :Mjoin(PATL,SyrkIntoC_aXbX_L2UT);
      }
      #ifdef TCPLX
         blk2c(N, N, alpha, rC, iC, beta, C, ldc);
         #ifdef Conj_  /* must zero complex part of diagonal! */
            Mjoin(PATLU,zero)(N, C+1, (ldc+1)SHIFT);
         #endif
      #else
         blk2c(N, N, alpha, pC, beta, C, ldc);
      #endif
   }
#endif
   free(vp);
}
#undef one
