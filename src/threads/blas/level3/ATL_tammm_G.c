#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
void Mjoin(PATL,ipCompBlk)
   (ipinfo_t *ip, size_t i, size_t j, size_t k, int be,
    TYPE *wA, TYPE *wB, TYPE *wC, TYPE *nA, TYPE *nB, TYPE *nC)
{
   #ifdef TCPLX
      TYPE *iC=wC, *rC=wC+ip->szC, *iA=wA, *iB=wB, *rA, *rB;
      ammkern_t amm_b1, amm_bn, amm0_b1, amm0_bn;
   #else
      ammkern_t amm;
   #endif
   ATL_CUINT kb = (k != ip->nfkblks) ? ip->kb : ip->KB0;
   const size_t nfmblks = ip->nfmblks, npmblks = ip->npmblks;
   const size_t nfnblks = ip->nfnblks, npnblks = ip->npnblks;
   unsigned int nmu, nnu;

   if (i+1 == nfmblks + npmblks)
      nmu = ip->nmuF;
   else
      nmu = (i < nfmblks) ? ip->nmu : ip->pnmu;
   if (j+1 == nfnblks + npnblks)
      nnu = ip->nnuF;
   else
      nnu = (j < nfnblks) ? ip->nnu : ip->pnnu;

   #ifdef TCPLX
      rA = iA + ((i < nfmblks) ? ip->szA : ip->pszA);
      rB = iB + ((j < nfnblks) ? ip->szB : ip->pszB);
      if (kb == ip->kb)  /* can use normal amm */
      {
         if (be)
         {
            amm0_b1 = amm_b1 = ip->amm_b1;
            amm0_bn = amm_bn = ip->amm_bn;
         }
         else
         {
            amm0_b1 = amm0_bn = ip->amm_b0;
            amm_b1 = ip->amm_b1;
            amm_bn = ip->amm_bn;
         }
      }
      else /* need K1 versions of everything */
      {
         if (be)
         {
            amm0_b1 = amm_b1 = ip->ammK1_b1;
            amm0_bn = amm_bn = ip->ammK1_bn;
         }
         else
         {
            amm0_b1 = amm0_bn = ip->ammK1_b0;
            amm_b1 = ip->ammK1_b1;
            amm_bn = ip->ammK1_bn;
         }
      }
      amm0_bn(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
      amm0_b1(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
      amm_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
      amm_b1(nmu, nnu, kb, iA, rB, iC, nA, nB, rC);
   #else
      if (kb == ip->kb)  /* can use normal amm */
         amm = (be) ? ip->amm_b1 : ip->amm_b0;
      else               /* must use ammK1 for K-cleanup */
         amm = (be) ? ip->ammK1_b1 : ip->ammK1_b0;
      amm(nmu, nnu, kb, wA, wB, wC, nA, nB, nC);
   #endif
}
void Mjoin(PATL,ipwriteC)(ipinfo_t *ip, size_t i, size_t j, const TYPE *beta,
                          TYPE *wC, TYPE *C)
{
   #ifdef TCPLX
      TYPE *rC = wC + ip->szC;
   #endif
   const size_t nfmblks=ip->nfmblks, nfnblks=ip->nfnblks;
   int mb, nb;

   C = IdxC_ip(ip, C, i, j);
   if (i+1 != nfmblks + ip->npmblks)
      mb = (i < nfmblks) ? ip->mb : ip->pmb;
   else
      mb = ip->mF;
   if (j+1 != nfnblks + ip->npnblks)
      nb = (j < nfnblks) ? ip->nb : ip->pnb;
   else
      nb = ip->nF;
   #ifdef TCPLX
      if (!beta)
         ip->blk2c_b1(mb, nb, ip->alpC, rC, wC, ip->ONE, C, ip->ldc);
      else
         ip->blk2c(mb, nb, ip->alpC, rC, wC, beta, C, ip->ldc);
   #else
      if (!beta)
         ip->blk2c_b1(mb, nb, ip->alpC, wC, ATL_rone, C, ip->ldc);
      else
         ip->blk2c(mb, nb, ip->alpC, wC, *beta, C, ip->ldc);
   #endif
}
/*
 * This routine computes one block C += A * B step
 */
static void DoBlkWtCpy
(
   ATL_tamm_gOOO_t *pd,              /* problem definition */
   ATL_CUINT rank,                   /* my (virtual) rank */
   size_t dblk,                      /* CpyBlk of C to compute */
   size_t k,                         /* which K-dim blk to do */
   int ibeta,                        /* 0: assume beta=0, else beta=1 */
   TYPE *wC                          /* private C workspace */
)
{
   ipinfo_t *ip=pd->ip;
   const size_t nmblks=pd->nmblks, nnblks=pd->nnblks, nkblks=ip->nfkblks+1;
   size_t i=dblk, j=dblk;
   TYPE *wA = pd->wA, *wB = pd->wB;

   if (dblk < nmblks)  /* need to copy A blk */
      wA = Mjoin(PATL,ipcopyA)(ip, pd->A, dblk, k, wA);
   else
   {
      i = nmblks-1;
      wA = IdxAw_ip(pd->ip, wA, i, k);
   }
   if (dblk < nnblks)  /* need to copy B blk */
      wB = Mjoin(PATL,ipcopyB)(ip, pd->B, k, dblk, wB);
   else
   {
      j = nnblks-1;
      wB = IdxBw_ip(pd->ip, wB, k, j);
   }
/*
 * If I just finished copying dblk K-panel for A & B, let folks know
 */
   if (ATL_DecGlobalAtomicCountDown((pd->KdonCtr)[dblk], rank) == 1)
   {
      ATL_mutex_lock(pd->cpmut);
      if (dblk < nmblks)
         ATL_SetBitBV(pd->cpyAdBV, dblk);
      if (!pd->cpyAdone)
         if (ATL_FindFirstUnsetBitBV(pd->cpyAdBV, 0) == -1)
            pd->cpyAdone = 1;

      if (dblk < nnblks)
         ATL_SetBitBV(pd->cpyBdBV, dblk);
      if (!pd->cpyBdone)
         if (ATL_FindFirstUnsetBitBV(pd->cpyBdBV, 0) == -1)
            pd->cpyBdone = 1;
      ATL_mutex_unlock(pd->cpmut);
   }
   if (dblk < Mmin(nmblks, nnblks))
      Mjoin(PATL,ipCompBlk)(ip, i, j, k, ibeta, wA, wB, wC, wA, wB, wC);
}
static void CompDiagC                /* old DoBlksWtCpy */
(
   ATL_tamm_gOOO_t *pd,              /* problem definition */
   unsigned const int rank,          /* my (virtual) rank */
   size_t dblk,                      /* CpyBlk of C to compute */
   size_t kctr,                      /* non-zero KbegCtr */
   TYPE *wC                          /* my private C workspace */
)
/*
 * This routine is called only at the beginning of the algorithm, and its
 * main purpose is to perform the data copy of both A & B.  In the process,
 * we compute the diagonal blocks of C, so that the thread doing the block
 * gets to use the data it copies at least once, and so that we cut bus
 * traffic at least a little bit by doing computation.  Large-scale systems
 * are likely to slow down to speed of memory during this step regardless.
 */
{
   const size_t nkblks = pd->ip->nfkblks + 1;
   const size_t minblks = Mmin(pd->nmblks, pd->nnblks);

   DoBlkWtCpy(pd, rank, dblk, nkblks-kctr, 0, wC);
   while (kctr = ATL_DecGlobalAtomicCount((pd->KbegCtr)[dblk], rank))
   {
      DoBlkWtCpy(pd, rank, dblk, nkblks-kctr, 1, wC);
   }
/*
 * If we computed a block of C (rather than copying A or B past C's diagonals),
 * seize mutex for block of original C, and add my part to it
 */
   if (dblk < minblks)
   {
      int COPIED=0;
      ATL_mutex_lock(pd->Cmuts[dblk]);
      if (!ATL_IsBitSetBV(pd->cbetaBV, dblk))
      {
         ATL_mutex_lock(pd->cbetamut);
         if (!ATL_IsBitSetBV(pd->cbetaBV, dblk))
         {
            ATL_SetBitBV(pd->cbetaBV, dblk);
            ATL_mutex_unlock(pd->cbetamut);
            Mjoin(PATL,ipwriteC)(pd->ip, dblk, dblk, SADD(pd->beta), wC, pd->C);
            COPIED=1;
         }
         else
            ATL_mutex_unlock(pd->cbetamut);
      }
      if (!COPIED)
         Mjoin(PATL,ipwriteC)(pd->ip, dblk, dblk, NULL, wC, pd->C);
      ATL_mutex_unlock(pd->Cmuts[dblk]);
   }
}

static void DoCompWithCopy(ATL_tamm_gOOO_t *pd, int rank, TYPE *wC)
/*
 * This routine loops over the diagonal blocks of C in order to copy A & B.
 * All threads first try to find a new diag blk to work on alone.  If no
 * diagonal blocks are left, then they find a diagonal block with work still
 * to be done and steal work from the current diagonal block owner.
 */
{
   int NEWBLK=1;
   const size_t ND = pd->nMNblks;  /* max K-panels in A or B for copying */
/*
 * Keep going as long as there is copy work to be done
 */
   while (!pd->NOCPWORK)
   {
       size_t d=0, k;
/*
 *     Find which diagonal block to work on, and then which k blk to use
 */
       if (NEWBLK)
       {
          d = ATL_DecGlobalAtomicCount(pd->ccCtr, rank);
          if (d)
          {
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[ND-d], rank);
             if (!k)     /* if no more K work to do */
                d = 0;   /* can't work on this diag after all */
          }
       }
/*
 *     If all diagonal blocks currently being worked on by threads, find
 *     one that I can help with.
 */
       if (!d)
       {
          size_t i;
          NEWBLK = 0;
          for (i=0; i < ND; i++)
          {
             size_t j = (i+rank)%ND;
             k = ATL_DecGlobalAtomicCount((pd->KbegCtr)[j], rank);
             d = ND-j;
             if (k)
                goto FOUNDDK;
          }
          pd->NOCPWORK = 1;   /* no work left to assign */
          return;             /* so done with copy of A&B and this func */
       }
/*
 *     If I reach here, I've got a valid d & k;  and I'll call a routine
 *     that continues to grab blocks from this diag along K until all K
 *     is done; it will then write the answer back to the original C, and
 *     return to this loop to see if it can help with another diag.
 */
       FOUNDDK:
          CompDiagC(pd, rank, ND-d, k, wC);
   }
}

static long FindCblk(const int rank, ATL_tamm_gOOO_t *pd)
/*
 * RETURNS: undone C block that has had its A and B k-panels copied or:
 *          -1: no work left to do;     -2: no work available
 * NOTE: this routine only called when A/B copying is still ongoing!
 */
{
   unsigned const int ncblks=pd->nCblks, nmblks=pd->nmblks, nnblks=pd->nnblks;
   int bb=0, ib;
   do
   {
      int i, j;
      ib = ATL_FindFirstUnsetBitBV(pd->cCblkBV, bb);
      if (ib == -1)
         return(-1);
      j = ib / nmblks;
      i = ib - j*nmblks;
      if (pd->cpyAdone || ATL_IsBitSetBV(pd->cpyAdBV, i))
      {
         if (pd->cpyBdone || ATL_IsBitSetBV(pd->cpyBdBV, j))
            return(ib);
      }
      bb = ib+1;
   }
   while(bb < ncblks);
   return(-2);
}
/*
 * Returns a number between [0,pd->ncK-1].  This is the Kctr to use for
 * work on this C, and it is the distance-1 from the last diagonal block in C.
 * Computes the (i,j) block coordinate of C block to use.
 * RETURNS: which KlastCtr to use.
 */
static INLINE int GetResCblk
   (ATL_CUINT rank, ATL_tamm_gOOO_t *pd, size_t *I, size_t *J)
{
   ipinfo_t *ip=pd->ip;
   const unsigned int nres = pd->ncK;  /* # of reserved C blks */
   const size_t nmblks=ip->nfmblks+ip->npmblks, nnblks=ip->nfnblks+ip->npnblks;
   unsigned int ctr, k;
/*
 * If copy still ongoing, assign C work based on dep info in cpy[A,B]dBV
 */
   while (!(pd->cpyAdone | pd->cpyBdone))
   {
      size_t i, j;
      for (k=0; k < nres; k++)
      {
         if (!ATL_GetGlobalAtomicCount(pd->KlastCtr[k], rank))
            continue;
         if (nmblks > nnblks) /* last row (non-diag) is reserved */
         {
            i = nmblks-1;
            j = k;
         }
         else               /* nnblks >= nmblks */
         {                  /* last col (non-diag) is reserved */
            i = k;
            j = nnblks-1;
         }
         if (!pd->cpyAdone && !ATL_IsBitSetBV(pd->cpyAdBV, i))
            continue;
         if (!pd->cpyBdone && !ATL_IsBitSetBV(pd->cpyBdBV, j))
            continue;
/*
 *       For weakly ordered caches, seize mutex to force sync
 */
         #if ATL_CBC_STRONG == 0
            ATL_mutex_lock(pd->cpmut);
            ATL_mutex_unlock(pd->cpmut);
         #endif
         *I = i;
         *J = j;
         return(k);
      }
   }
/*
 * If we reach here, done copying, so select first reserved C blk wt work
 */
   for (k=0; k < nres; k++)
   {
      ctr = (k + rank)%nres;
      if (ATL_GetGlobalAtomicCount(pd->KlastCtr[ctr], rank))
      {
         if (nmblks > nnblks) /* last row (non-diag) is reserved */
         {
            *I = nmblks-1;
            *J = ctr;
         }
         else               /* nnblks >= nmblks */
         {                  /* last col (non-diag) is reserved */
            *I = ctr;
            *J = nnblks-1;
         }
         return(ctr);
      }
   }
   return(-1);  /* no work left! */
}
/*
 * RETURNS: C block that has had its A and B k-panels copied
 */
static INLINE long GetCblk(ATL_CUINT rank, ATL_tamm_gOOO_t *pd)
{
   const int ncblks = pd->nCblks;
   int ib;
/*
 * If we are done copying, then use the counter to get the C block to work on
 * NOTE: thread that sets cpy[A,B]done must have cpmut to avoid race cond!
 */
   if (pd->cpyAdone & pd->cpyBdone)
   {
/*
 *    On weakly-ordered caches, checking cpyA/cpyB insufficient to ensure
 *    copied data is visible, so lock a BS mutex to force a memory barrier
 *    (the producer thread will have locked a mutex to update the copy bitvec)
 *    Use diagonal C mutex, which should have nobody using it.
 */
      #if !ATL_CBC_STRONG
         static char FIRST=1;
      if (FIRST)
      {
         FIRST=0;
         ib = (rank < pd->nMNblks)  ? rank : 0; /* spread across dead mutexes */
         ATL_mutex_lock(pd->Cmuts[ib]);
         ATL_mutex_unlock(pd->Cmuts[ib]);
      }
      #endif
      do
      {
         ib = ATL_DecGlobalAtomicCount(pd->cCtr, rank);
         if (!ib)
            return(-1);
         ib = ncblks - ib;
      }
      while (ATL_IsBitSetBV(pd->cCblkBV, ib));
      return(ib);
   }
/*
 * If we reach here, we must assign C work based on dep info in cpy[A,B]dBV
 */
   else
   {
      while(1)
      {
         ib = FindCblk(rank, pd);
/*
 *       If we have a candidate block, grab the mutex and make sure one is
 *       still there
 */
         if (ib >= 0)
         {
            ATL_mutex_lock(pd->cpmut);
/*
 *          If copy finished while locking, call ourselves to get fast answer
 */
            if (pd->cpyAdone & pd->cpyBdone)
            {
               ATL_mutex_unlock(pd->cpmut);
               return(GetCblk(rank, pd));
            }
/*
 *          Otherwise, see if we've still got a candidate block
 */
            ib = FindCblk(rank, pd);
            if (ib >= 0)  /* we've got a block! */
            {
               ATL_SetBitBV(pd->cCblkBV, ib);  /* claim the block */
               ATL_mutex_unlock(pd->cpmut);
               return(ib);
            }
            ATL_mutex_unlock(pd->cpmut);
         }
         if (ib == -1)
            return(-1);
/*
 *       If no work available, pause before trying again
 */
         if (ib == -2)
            ATL_thread_yield();
      }
   }
}

static void DoCompNoCopy(ATL_tamm_gOOO_t *pd, ATL_UINT rank, TYPE *wC)
{
   ipinfo_t *ip=pd->ip;
   TYPE *C=pd->C, *wA=pd->wA, *wB=pd->wB;
   #ifdef TCPLX
      const TYPE *beta=pd->beta;
      TYPE *rC = wC + ip->szC;
   #else
      const TYPE beta=pd->beta;
      #define rC wC
   #endif
   size_t nmblks = pd->nmblks;
   size_t ic;
   const ablk2cmat_t blk2c=ip->blk2c;

   while ((ic = GetCblk(rank, pd)) != -1)
   {
      int i, j;
      j = ic / nmblks;
      i = ic - j*nmblks;
      Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, IdxC_ip(ip, C, i, j), 3,
                           IdxAw_ip(ip, wA, i, 0), IdxBw_ip(ip, wB, 0, j),
                           rC, wC, beta, blk2c);
   }
}
#ifndef TCPLX
   #undef rC
#endif
static void DoLastComp(ATL_tamm_gOOO_t *pd, ATL_UINT rank, TYPE *wC)
{
   ipinfo_t *ip=pd->ip;
   TYPE *wA=pd->wA, *wB=pd->wB, *C=pd->C;
   void **Cmuts = pd->Cmuts + pd->nMNblks;
   const size_t nkblks = ip->nfkblks+1, nmblks=ip->nfmblks+ip->npmblks;
   size_t i, j;
   int cctr;

/*
 * As long as I can find C blks using the KlastCtr
 */
   cctr = GetResCblk(rank, pd, &i, &j);
   while (cctr >= 0)
   {
      size_t k;
      int WORKED=0;
/*
 *    Keep working on same C blk until we run out of A/B blocks in K-panel
 */
      while ( (k = ATL_DecGlobalAtomicCount(pd->KlastCtr[cctr], rank)) )
      {
         TYPE *a, *b;
         k = nkblks - k;
         a = IdxAw_ip(ip, wA, i, k);
         b = IdxBw_ip(ip, wB, k, j);
         Mjoin(PATL,ipCompBlk)(ip, i, j, k, WORKED, a, b, wC, a, b, wC);
         WORKED = 1;
      }
      if (WORKED)  /* if I did computation, will need to write C out */
      {
         int COPIED=0;
         ATL_mutex_lock(Cmuts[cctr]);
         if (!ATL_IsBitSetBV(pd->lbetaBV, cctr))
         {
            ATL_mutex_lock(pd->lbetamut);
            if (!ATL_IsBitSetBV(pd->lbetaBV, cctr))
            {
               ATL_SetBitBV(pd->lbetaBV, cctr);
               ATL_mutex_unlock(pd->lbetamut);
               Mjoin(PATL,ipwriteC)(ip, i, j, SADD(pd->beta), wC, C);
               COPIED=1;
            }
            else
               ATL_mutex_unlock(pd->lbetamut);
         }
         if (!COPIED)
            Mjoin(PATL,ipwriteC)(pd->ip, i, j, NULL, wC, C);
         ATL_mutex_unlock(Cmuts[cctr]);
/*printf("%d:wrote C(%d,%d), CPY=%d, ctr=%d\n", rank, i, j, COPIED, cctr);*/
      }
      cctr = GetResCblk(rank, pd, &i, &j);
   }
}

void Mjoin(PATL,DoWork_OOO)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_gOOO_t *pd = pp->PD;
   TYPE *myC = pd->wC + vrank*((pd->ip->szC)SHIFT);
   if (!pd->NOCPWORK)
      DoCompWithCopy(pd, vrank, myC);
   DoCompNoCopy(pd, vrank, myC);
   if (pd->ncK)
      DoLastComp(pd, vrank, myC);
}
/*
 * This algorithm can be used for any parallelizable case wt K>LASTKB.  Usually,
 * it deals out blocks of C, with K-dim divided only for initial copy of A & B
 * and some load balance at the end.  However, if C is small and K is large,
 * it can switch and do parallalize along K (with fairly high overhead).
 * It tries to copy all of A & B up front, so recursion may be needed to
 * use it for large matrices.
 * RETURNS: 0 if it did the operation, non-zero if it did not
 */

int Mjoin(PATL,tammm_G)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   size_t M,
   size_t N,
   size_t K,
   const SCALAR alpha,
   const TYPE *A,
   size_t lda,
   const TYPE *B,
   size_t ldb,
   const SCALAR beta,
   TYPE *C,
   size_t ldc
)
{
   ATL_tamm_gOOO_t pd;
   ipinfo_t ip;
   int P, i, j, k, idx, mb, pmb, nb, pnb;
   size_t nkcnt, nkblks, nKpar, nMNblks, sz, szA, szB;
   amminfo_t mminfo;
   void *vp;
   #ifdef ATL_PHI_SORTED
      const unsigned int nthr = ATL_TP_PTR ? ATL_TP_PTR->nthr : ATL_NTHREADS,
                         p4 = nthr>>2, p4_2 = nthr+nthr, p4_3 = p4_2+p4;
   #endif

   #if 1
   Mjoin(PATL,ipgenInfo)(&ip, 0, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   #else
   idx = Mjoin(PATL,tGetParCIndx)(&ip, ATL_NTHREADS, M, N, K);
   Mjoin(PATL,geFillInIPInfo)(&ip, idx, TA, TB, M, N, K, lda, ldb, ldc,
                              alpha, beta,
                              ip.nfmblks, ip.npmblks, ip.mb, ip.pmb,
                              ip.nfnblks, ip.npnblks, ip.nb, ip.pnb);
   #endif
/*   if (ip.kb > K)
        return(1); */

   pd.ip = &ip;
   pd.nmblks = ip.nfmblks + ip.npmblks;
   pd.nnblks = ip.nfnblks + ip.npnblks;
   pd.nCblks = sz = ((size_t)pd.nmblks)*pd.nnblks;
   if (pd.nCblks != sz)
      return(1);
   nkblks = ip.nfkblks + 1;
   pd.ncK = 0;
   if (pd.nmblks >= pd.nnblks)
   {
      pd.nMNblks = nMNblks = pd.nmblks;
      nKpar = pd.nnblks;
   }
   else
   {
      pd.nMNblks = nMNblks = pd.nnblks;
      nKpar = pd.nmblks;
   }
   pd.beta = beta;
/*
 * If K is large, dealing out C blks alone for parallelism can lead to idle
 * time at end of algorithm, as a long K-loop is given out as the last task.
 * Therefore, K is decent size, reserve roughly P-1 C blks that can parallize
 * along the K dim, providing for reduced idle time at the cost of extra
 * parallel idle time.
 */
   P = nkcnt = pd.nCblks - nKpar;  /* # of C blks given out as private work */
   P += nkblks * nKpar;            /* # of C blks parallelized along K */
   P = Mmin(ATL_NTHREADS,P);
   if (nkcnt > 1 && ip.nfkblks > 8)
   {        /* do we want to split K dim of some of these C blks? */
      size_t imin;
/*
 *    Our algorithm can only give out the non-diagonal elements of last row/col
 *    so we must limit ncK (# of non-diag c parallelized along K) to last
 *    row/col block count
 */
      if (pd.nmblks > pd.nnblks)
         imin = pd.nnblks;
      else
         imin = (pd.nmblks != pd.nnblks) ? pd.nmblks : pd.nmblks-1;
      if (pd.nCblks-nKpar < P)
         pd.ncK = imin;
      else
      {
         pd.ncK = (P > 2) ? P-2 : 1;
         pd.ncK = Mmin(pd.ncK,imin);
      }
      P += pd.ncK*nkblks;
      P = Mmin(ATL_NTHREADS,P);
   }
/*
 * Quick exit for problems too small to thread
 */
   if (P < 2)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
#if 0
printf("P=%d, ncK=%d, M=(%d,%d; %d,%d) N=(%d,%d; %d,%d)\n",
       P, (int)pd.ncK, (int)ip.nfmblks, (int)ip.npmblks, ip.mb, ip.pmb,
       (int)ip.nfnblks, (int)ip.npnblks, ip.nb, ip.pnb);
#endif
   szA = nkblks*(ip.nfmblks*ip.szA + ip.npmblks*ip.pszA);
   szB = nkblks*(ip.nfnblks*ip.szB + ip.npnblks*ip.pszB);
   sz = (nMNblks*3+2*pd.ncK)*sizeof(void*);/* K[beg,don,last]Ctr & Cmuts arrs */
   sz += ATL_MulBySize(ip.szC*P + szA + szB + (ip.mu<<1)*ip.nu);
   sz += 3*ATL_Cachelen;                  /* room for alignment */
   #ifdef ATL_PHI_SORTED
      pd.ncntxts = 1;
      pd.ncores = p4;
      if (P == p4_2)
         pd.ncntxts = 2;
      else if (P == p4_3)
         pd.ncntxts = 3;
      else if (P == nthr)
         pd.ncntxts = 4;
      else
         pd.ncores = P;
      if (pd.ncntxts > 1)
         sz += ATL_Cachelen*(pd.ncores+1);
   #endif
   if (sz > ATL_PTMAXMALLOC)
      return(2);
   vp = malloc(sz);
   if (!vp)
      return(2);
   #ifdef ATL_PHI_SORTED
   if (pd.ncntxts > 1)
   {
      pd.chkin = (volatile int*)ATL_AlignPtr(vp);
      pd.KbegCtr = (void*)(((char*)pd.chkin)+pd.ncores*ATL_Cachelen);
      for (i=0; i < pd.ncores; i++)
      {
         int *ip;
         char *cp = (char*)pd.chkin;

         cp += i*ATL_Cachelen;
         ip = (int*) cp;
         *ip = ip[1] = ip[2] = ip[3] = -2;
      }
   }
   else
   #endif
   pd.KbegCtr = vp;
   pd.cpyAdone = pd.cpyBdone = pd.NOCPWORK = 0;
   pd.ccCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nMNblks, P),nMNblks, 0);
   pd.cCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nCblks, P),pd.nCblks, 0);
   nkcnt = ATL_EstNctr(nkblks, P);
   pd.KdonCtr = pd.KbegCtr + nMNblks;
   pd.Cmuts = pd.KdonCtr + nMNblks;
   if (pd.ncK)
   {
      pd.KlastCtr = pd.Cmuts+nMNblks+pd.ncK;
      pd.wA = (TYPE*)(pd.KlastCtr+pd.ncK);
   }
   else
   {
      pd.KlastCtr = NULL;
      pd.wA = (TYPE*)(pd.Cmuts+nMNblks);
   }
   pd.wA = ATL_AlignPtr(pd.wA);
   pd.wB = pd.wA + (szA SHIFT);
   pd.wB = ATL_AlignPtr(pd.wB);
   pd.wC = pd.wB + (szB SHIFT);
   pd.wC = ATL_AlignPtr(pd.wC);
   pd.A = A;
   pd.B = B;
   pd.C = C;
   for (i=0; i < nMNblks; i++)
   {
      pd.KbegCtr[i] = ATL_SetGlobalAtomicCount(nkcnt, nkblks, 0);
      pd.KdonCtr[i] = ATL_SetGlobalAtomicCountDown(nkcnt, nkblks);
      pd.Cmuts[i] = ATL_mutex_init();
   }
   for (i=0; i < pd.ncK; i++)
   {
      pd.KlastCtr[i] = ATL_SetGlobalAtomicCount(nkcnt, nkblks, 0);
      pd.Cmuts[i+nMNblks] = ATL_mutex_init();
   }
   pd.cpmut = ATL_mutex_init();
   pd.cbetamut = ATL_mutex_init();
   if (pd.ncK)
      pd.lbetamut = ATL_mutex_init();
   else
      pd.lbetamut = NULL;
   pd.cpyAdBV = ATL_NewBV(pd.nmblks);
   pd.cpyBdBV = ATL_NewBV(pd.nnblks);
   pd.cCblkBV = ATL_NewBV(pd.nCblks);
   pd.cbetaBV = ATL_NewBV(nMNblks);
   if (pd.ncK)
      pd.lbetaBV = ATL_NewBV(pd.ncK);
   else
      pd.lbetaBV = NULL;
/*
 * Initialize cCblkBV so that all diagonal blocks are shown as already complete
 * This BV is used for assigning work to non-copy blocks.
 */
   k = Mmin(pd.nmblks, pd.nnblks);
   for (i=0; i < k; i++)
      ATL_SetBitBV(pd.cCblkBV, i*(pd.nmblks+1));
/*
 * Reserve blks for load balancing at end; take largest of last row/col-panel
 */
   if (pd.nmblks > pd.nnblks) /* reserve last row panel of C */
   {
      size_t k = pd.nmblks-1;
      for (i=0; i < pd.ncK; i++, k += pd.nmblks)
         ATL_SetBitBV(pd.cCblkBV, k);
   }
   else                      /* nnblks >= nmblks */
   {                         /* reserve last col panel of C */
      size_t k = pd.nCblks-pd.nmblks;
      for (i=0; i < pd.ncK; i++)
         ATL_SetBitBV(pd.cCblkBV, k+i);
   }
   ATL_goParallel(P, Mjoin(PATL,DoWork_OOO), NULL, &pd, NULL);
/*
 * Free allocated structures and return;
 */
   ATL_FreeBV(pd.cpyAdBV);
   ATL_FreeBV(pd.cpyBdBV);
   ATL_FreeBV(pd.cCblkBV);
   ATL_FreeBV(pd.cbetaBV);
   ATL_mutex_free(pd.cpmut);
   ATL_mutex_free(pd.cbetamut);
   ATL_FreeGlobalAtomicCount(pd.cCtr);
   ATL_FreeGlobalAtomicCount(pd.ccCtr);
   for (i=0; i < nMNblks; i++)
   {
      ATL_FreeGlobalAtomicCount(pd.KbegCtr[i]);
      ATL_FreeGlobalAtomicCountDown(pd.KdonCtr[i]);
      ATL_mutex_free(pd.Cmuts[i]);
   }
   for (i=0; i < pd.ncK; i++)
   {
      ATL_FreeGlobalAtomicCount(pd.KlastCtr[i]);
      ATL_FreeGlobalAtomicCount(pd.Cmuts[i+nMNblks]);
   }
   free(vp);
   return(0);
}
