#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#include Mstr(Mjoin(ATLAS_PRE,sysinfo.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_sum.h))
#include Mstr(Mjoin(ATLAS_PRE,ipmen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_syrk.h))
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"


static void DoWork(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tsyrk_tN_t *pd = pp->PD;
   const TYPE *A = pd->A;
   #ifdef TCPLX
      TYPE *iA = pd->w + vrank*pd->szW, *rA=iA+pd->szA, *iC=rA+pd->szA, *rC, *a;
      #ifdef Conj_
         TYPE *crA, *ciA;
      #endif
      ammkern_t syrkK_bn = Mjoin(PATL,amsyrkK_b0);
      ammkern_t syrkK_b1 = Mjoin(PATL,amsyrkK_b0);
      #define one pd->ONE
   #else
      TYPE *wA = pd->w + vrank*pd->szW, *wC = wA + pd->szA;
      ammkern_t syrkK = Mjoin(PATL,amsyrkK_b0);
      #define one ATL_rone
   #endif
   cm2am_t a2blk = pd->a2blk;
   const size_t njobs=pd->njobs, lda=pd->lda, nK1=pd->nK1, incAk=pd->incAk;
   ATL_CUINT kb=pd->kb, N=pd->N, nnu=pd->nnu;

   #ifdef TCPLX
      iC = ATL_AlignPtr(iC);
      rC = iC + pd->szC;
      #ifdef Conj_
         crA = (pd->NoTrans) ? rA : iA;
         ciA = (pd->NoTrans) ? iA : rA;
      #endif
   #else
      wC = ATL_AlignPtr(wC);
   #endif
/*
 * First do the big chunky jobs
 */
   if (njobs)
   {
      void *jobCtr = pd->jobCtr;
      size_t Jctr, njobs=pd->njobs;
      ATL_CUINT jobshift=pd->jobshift;
      while ( (Jctr = ATL_DecGlobalAtomicCount(jobCtr, vrank)) )
      {
         size_t k = (njobs-Jctr)<<jobshift;
         const size_t kend = k + (1<<jobshift);
         do
         {
            #ifdef TCPLX
               a2blk(kb, N, one, A+incAk*k, lda, rA, iA);
               #ifdef Conj_
                  syrkK_b1(nnu, nnu, kb, iA, iA, rC, crA, ciA, iC);
                  syrkK_bn(nnu, nnu, kb, crA, ciA, iC, rA, rA, rC);
                  Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, rA, rA, rC, ciA,crA,iC);
                  Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, ciA, crA, iC, iA,iA,rC);
               #else
                  syrkK_bn(nnu, nnu, kb, iA, iA, rC, rA, iA, iC);
                  syrkK_b1(nnu, nnu, kb, rA, iA, iC, rA, rA, rC);
                  Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, rA, rA, rC, iA, rA, iC);
                  Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, iA, rA, iC, iA, iA, rC);
               #endif
               syrkK_bn = Mjoin(PATL,amsyrkK_bn);
               syrkK_b1 = Mjoin(PATL,amsyrkK_b1);
            #else
               a2blk(kb, N, ATL_rone, A+incAk*k, lda, wA);
               syrkK(nnu, nnu, kb, wA, wA, wC, wA, wA, wC);
               syrkK = Mjoin(PATL,amsyrkK_b1);
            #endif
         }
         while (++k != kend);
      }
   }
/*
 * Now loop over individual full-kb computations
 */
   if (nK1)
   {
      const size_t k0 = (njobs << pd->jobshift);  /* # of k blks already done */
      size_t k, Kctr;
      void *K1ctr=pd->K1ctr;
      while ( (Kctr = ATL_DecGlobalAtomicCount(K1ctr, vrank)) )
      {
         k = k0 + nK1 - Kctr;
         #ifdef TCPLX
            a2blk(kb, N, one, A+incAk*k, lda, rA, iA);
            #ifdef Conj_
               syrkK_b1(nnu, nnu, kb, iA, iA, rC, crA, ciA, iC);
               syrkK_bn(nnu, nnu, kb, crA, ciA, iC, rA, rA, rC);
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, rA, rA, rC, ciA, crA, iC);
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, ciA, crA, iC, iA, iA, rC);
            #else
               syrkK_bn(nnu, nnu, kb, iA, iA, rC, rA, iA, iC);
               syrkK_b1(nnu, nnu, kb, rA, iA, iC, rA, rA, rC);
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, kb, rA, rA, rC, iA, rA, iC);
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, kb, iA, rA, iC, iA, iA, rC);
            #endif
            syrkK_bn = Mjoin(PATL,amsyrkK_bn);
            syrkK_b1 = Mjoin(PATL,amsyrkK_b1);
         #else
            a2blk(kb, N, ATL_rone, A+incAk*k, lda, wA);
            syrkK(nnu, nnu, kb, wA, wA, wC, wA, wA, wC);
            syrkK = Mjoin(PATL,amsyrkK_b1);
         #endif
      }
   }
/*
 * If last kb is partial, only one guy does that
 */
   if (pd->kb0)
   {
      if (ATL_DecAtomicCount(pd->KB0ctr))
      {
         size_t k = (njobs<<pd->jobshift) + nK1;
         #ifdef TCPLX
            ATL_CUINT KB=pd->KB0;
            a2blk(pd->kb0, N, one, A+incAk*k, lda, rA, iA);
            #ifdef Conj_
               syrkK_b1(nnu, nnu, KB, iA, iA, rC, crA, ciA, iC);
               syrkK_bn(nnu, nnu, KB, crA, ciA, iC, rA, rA, rC);
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, KB, rA, rA, rC, ciA, crA, iC);
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, KB, ciA, crA, iC, iA, iA, rC);
            #else
               syrkK_bn(nnu, nnu, KB, iA, iA, rC, rA, iA, iC);
               syrkK_b1(nnu, nnu, KB, rA, iA, iC, rA, rA, rC);
               Mjoin(PATL,amsyrkK_bn)(nnu, nnu, KB, rA, rA, rC, iA, rA, iC);
               Mjoin(PATL,amsyrkK_b1)(nnu, nnu, KB, iA, rA, iC, iA, iA, rC);
               syrkK_b1 = NULL;
            #endif
         #else
            a2blk(pd->kb0, N, ATL_rone, A+incAk*k, lda, wA);
            syrkK(nnu, nnu, pd->KB0, wA, wA, wC, wA, wA, wC);
            syrkK = NULL;
         #endif
      }
   }
/*
 * Write out if we did any work
 */
#ifdef TCPLX
   if (syrkK_b1 != Mjoin(PATL,amsyrkK_b0)) /* we did at least 1 blk */
#else
   if (syrkK != Mjoin(PATL,amsyrkK_b0))
#endif
   {
      if (pd->wT)  /* Upper triangle must transpose in workspace */
      {
         ablk2cmat_t blk2c=pd->blk2c_b1;
         TYPE *wc = pd->wT + (vrank SHIFT) * ((size_t)N)*N, *C=pd->C;
         ATL_CSZT ldc=pd->ldc;
         ATL_UINT k, APPBETA=0;

         if (blk2c == Mjoin(PATL,SyrkIntoC_a1_b1))
         blk2c = Mjoin(PATL,SyrkIntoC_a1_b0);
         else
            blk2c = (blk2c == Mjoin(PATL,SyrkIntoC_an_b1)) ?
               Mjoin(PATL,SyrkIntoC_an_b0) : Mjoin(PATL,SyrkIntoC_aX_b0);
         #ifdef TCPLX
         #else
            blk2c(N, N, pd->alpha, wC, ATL_rzero, wc, N);
         #endif
         ATL_mutex_lock(pd->Cmut);
         if (pd->appBeta)  /* only set if beta != 1 & we are first here */
         {
            const SCALAR beta = pd->beta;
            if (SCALAR_IS_ZERO(pd->beta))
               for (k=0; k < N; k++, wc += (1 SHIFT), C += (ldc SHIFT))
                   Mjoin(PATL,copy)(k+1, wc, N, C, 1);
            else
               for (k=0; k < N; k++, wc += (1 SHIFT), C += (ldc SHIFT))
                  Mjoin(PATL,axpby)(k+1, one, wc, N, beta, C, 1);
            pd->appBeta = 0;
         }
         else
            for (k=0; k < N; k++, wc += (1 SHIFT), C += (ldc SHIFT))
               Mjoin(PATL,axpby)(k+1, one, wc, N, one, C, 1);
         ATL_mutex_unlock(pd->Cmut);
      }
      else   /* Lower triangle can just write to C */
      {
         #ifdef Conj_
            ATL_UINT ZeroDiag=0;
         #endif
         ATL_mutex_lock(pd->Cmut);
         if (pd->appBeta)
         {
            pd->appBeta = 0;
            #ifdef TCPLX
               pd->blk2c(N, N, pd->alpha, rC, iC, pd->beta, pd->C, pd->ldc);
            #else
               pd->blk2c(N, N, pd->alpha, wC, pd->beta, pd->C, pd->ldc);
            #endif
         }
         else
         {
            #ifdef TCPLX
               pd->blk2c_b1(N, N, pd->alpha, rC, iC, pd->beta, pd->C, pd->ldc);
            #else
               pd->blk2c_b1(N, N, pd->alpha, wC, pd->beta, pd->C, pd->ldc);
            #endif
         }
         #ifdef Conj_
            if (!(--(pd->Pcnt)))
            #ifdef Conj_  /* must zero complex part of diagonal! */
               Mjoin(PATLU,zero)(N, pd->C+1, (pd->ldc+1)SHIFT);
            #endif
         #endif
         ATL_mutex_unlock(pd->Cmut);
      }
   }  /* end of if we did work */
}
#undef one

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
void Mjoin(PATL,therk_tN)
#else
void Mjoin(PATL,tsyrk_tN)
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
   ATL_tsyrk_tN_t pd;
   void *vp;
   #ifdef TCPLX
      const TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #endif
   #ifdef Conj_
      const TYPE alpha[2] = {ralpha, ATL_rzero};
      const TYPE beta[2] = {rbeta, ATL_rzero};
   #endif
   size_t szA, szC, szT, szE, szP, incAk, nkb, k;
   ATL_CUINT nnu = (N+ATL_SYRKK_NU-1)/ATL_SYRKK_NU, NN=nnu*ATL_SYRKK_NU;
   ATL_UINT kb, P;
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
/*
 * Quick return if problem isn't worth parallelizing
 */
   #if 0
   if (K < (Mmax(ATL_sqAMM_66KB,ATL_SYRKK_KU)*Mmin(ATL_NTHREADS,4)) ||
       N == 1)
   {
      #ifdef Conj_
         Mjoin(PATL,herk_IP)(Uplo, TA, N, K, ralpha, A, lda, rbeta, C, ldc);
      #else
         Mjoin(PATL,syrk_IP)(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
      #endif
       return;
   }
   #endif
   idx = N2Idx(NN);
   kb = ATL_GetVWsyrkKB(idx);
   szC = ((ATL_SYRKK_NU*ATL_SYRKK_NU+ATL_SYRKK_VLEN-1)/ATL_SYRKK_VLEN)
         * ATL_SYRKK_VLEN;
   szC *= ((nnu+1)*nnu)>>1;  /* only need lower tri blks, not full nnu*nnu */
   if (K > kb)
   {
      nkb = K / kb;
      pd.kb0 = K - nkb*kb;
      if (!kb0)
      {
         pd.kb0 = kb;
         nkb--;
      }
      KB0 = ((pd.kb0+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   }
   else
   {
      nkb = 0;
      pd.kb0 = K;
      kb = KB0 = ((pd.kb0+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
   }
   #if 0
/*
 * Now, see if we need to reduce kb to enhance parallelism, or if we have full
 */
   nkb = K / kb;
   if (nkb < ATL_NTHREADS && kb > ATL_sqAMM_66KB)  /* need to reduce */
   {
      do
      {
         kb -= ATL_SYRKK_KU;
         nkb = K / kb;
      }
      while (nkb < ATL_NTHREADS && kb > ATL_sqAMM_66KB);
   }
   #endif
   pd.kb = kb;
   pd.nnu = nnu;
   pd.szC = szC;
/*
 * Compute the number of jobs (consisting of multiple kblks) and number of
 * k blocks given out individually (for load balancing)
 */
   pd.N = N;
   #if 1
      if (nkb > ATL_NTHREADS*3)
      {
         size_t nj, nk1;
         for (k=1; k < 8; k++)
         {
            nk1 = (ATL_NTHREADS-1)<<k;
            if (nk1+nk1 > nkb)
               break;
            nj = (nkb - nk1) >> k;
            if (nj < 2*ATL_NTHREADS)
               break;
         }
         pd.jobshift = --k;
         pd.njobs = (nkb - ((ATL_NTHREADS-1)<<k)) >> k;
         pd.njobs = (pd.njobs/ATL_NTHREADS)*ATL_NTHREADS;
         pd.nK1 = nkb - (pd.njobs<<k);
      }
      else
      {
         pd.nK1 = nkb;
         pd.jobshift = pd.njobs = 0;
      }
      #if 0
      printf("nkb=%d, njob=%d, jobCnt=%d, nK1=%d\n",
             (int)nkb, (int)pd.njobs, 1<<pd.jobshift, (int)pd.nK1);
      #endif
   #elif 1
      pd.nK1 = nkb;
      pd.jobshift = pd.njobs = 0;
   #elif 0
      pd.njobs = nkb;
      pd.jobshift = 0;
      pd.nK1 = 0;
   #endif
   k = (pd.kb0) ? 1 : 0;
   k += pd.nK1 + pd.njobs << pd.jobshift;
   P = Mmin(ATL_NTHREADS, k);
   #if Conj_
      pd.Pcnt = P;
      pd.NoTrans = (TA == AtlasNoTrans);
   #endif
/*
 * Our SYRK C copy is always to Lower, so if output is Upper, will need
 * workspace to put the Lower part, before reflecting it to Upper.
 * Need max of this extra space, or preload distance.
 */
   if (Uplo == AtlasUpper)
   {
      szT = P * N * N;
      szE = (szT >= ATL_SYRKK_NU*ATL_SYRKK_NU*2)?ATL_SYRKK_NU*ATL_SYRKK_NU*2:0;
   }
   else
   {
      szT = 0;
      szE = ATL_SYRKK_NU*ATL_SYRKK_NU*2;
   }
   szA = nnu*ATL_SYRKK_NU * kb;
   szP = ATL_MulBySize(szA+szC) + ATL_Cachelen;
   szP = ATL_MulByCachelen(ATL_DivByCachelen(szP + ATL_Cachelen-1));
   vp = malloc(ATL_MulBySize(P*(szP+szT)+szE)+2*ATL_Cachelen);
   ATL_assert(vp);

   pd.szA = szA;
   pd.szW = szP = ATL_DivBySize(szP)SHIFT;
   pd.w = ATL_AlignPtr(vp);
   if (Uplo == AtlasUpper)
   {
      pd.wT = pd.w + szP*P;
      pd.wT = ATL_AlignPtr(pd.wT);
   }
   else
      pd.wT = NULL;
   pd.A = A;
   pd.C = C;
   pd.lda = lda;
   pd.ldc = ldc;
   pd.alpha = alpha;
   pd.beta = beta;
   #ifdef TCPLX
      pd.ONE = ONE;
   #endif
   pd.appBeta = !SCALAR_IS_ONE(beta);
   pd.Cmut = ATL_mutex_init();
   if (pd.kb0)
      pd.KB0ctr = ATL_SetAtomicCount(1);
   else
      pd.KB0ctr = NULL;
   if (pd.njobs)
      pd.jobCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.njobs,P),pd.njobs,0);
   else
      pd.jobCtr = NULL;
   if (pd.nK1)
      pd.K1ctr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nK1,P), pd.nK1, 0);
   else
      pd.K1ctr = NULL;


   if (IS_COLMAJ(TA))
   {
      pd.incAk = lda*(kb SHIFT);
      pd.a2blk = Mjoin(PATL,a2blk_syrkT);
   }
   else
   {
      pd.incAk = kb SHIFT;
      pd.a2blk = Mjoin(PATL,a2blk_syrkN);
   }
   if (Uplo == AtlasLower)
   {
      if (SCALAR_IS_NONE(alpha))
      {
         pd.blk2c_b1 = Mjoin(PATL,SyrkIntoC_aNb1);
         if (SCALAR_IS_ONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_aNb1);
         else if (SCALAR_IS_NONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_aNbN);
         else
            pd.blk2c = SCALAR_IS_ZERO(beta) ?  Mjoin(PATL,SyrkIntoC_an_b0) :
                       Mjoin(PATL,SyrkIntoC_aNbX);
      }
      else if (SCALAR_IS_ONE(alpha))
      {
         pd.blk2c_b1 = Mjoin(PATL,SyrkIntoC_a1b1);
         if (SCALAR_IS_ONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_a1b1);
         else if (SCALAR_IS_NONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_a1bN);
         else
            pd.blk2c = SCALAR_IS_ZERO(beta) ?  Mjoin(PATL,SyrkIntoC_a1b0) :
                       Mjoin(PATL,SyrkIntoC_a1bX);
      }
      else /* alpha = X */
      {
         pd.blk2c_b1 = Mjoin(PATL,SyrkIntoC_aXb1);
         if (SCALAR_IS_ONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_aXb1);
         else if (SCALAR_IS_NONE(beta))
            pd.blk2c = Mjoin(PATL,SyrkIntoC_aXbN);
         else
            pd.blk2c = SCALAR_IS_ZERO(beta) ?  Mjoin(PATL,SyrkIntoC_aX_b0) :
                       Mjoin(PATL,SyrkIntoC_aX_bX);
      }
   }
   else
   {
      if (SCALAR_IS_NONE(alpha))
         pd.blk2c = Mjoin(PATL,SyrkIntoC_an_b0);
      else
         pd.blk2c = SCALAR_IS_ONE(alpha) ?
                 Mjoin(PATL,SyrkIntoC_a1_b0):Mjoin(PATL,SyrkIntoC_aX_b0);
   }

   ATL_goParallel(P, DoWork, NULL, &pd, NULL);

   ATL_mutex_free(pd.Cmut);
   if (pd.KB0ctr)
      ATL_FreeAtomicCount(pd.KB0ctr);
   if (pd.K1ctr)
      ATL_FreeGlobalAtomicCount(pd.K1ctr);
   if (pd.jobCtr)
      ATL_FreeGlobalAtomicCount(pd.jobCtr);
   free(vp);
}
