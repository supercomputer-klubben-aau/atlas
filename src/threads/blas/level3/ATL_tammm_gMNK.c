#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_cbc.h"
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"

TYPE *Mjoin(PATL,ipcopyA)(ipinfo_t *ip, const TYPE *A, size_t i, size_t k,
                          TYPE *w)
{
   ATL_CUINT K = (k != ip->nfkblks) ? ip->kb : ip->kb0;
   ATL_UINT M = (i < ip->nfmblks) ? ip->mb : ip->pmb;
   if (i == ip->nfmblks + ip->npmblks - 1)
      M = ip->mF;
   #ifdef TCPLX
      TYPE *iw, *rw;
      A = IdxA_ip(ip, A, i, k);
      iw = IdxAw_ip(ip, w, i, k);
      rw = iw + ((i < ip->nfmblks) ? ip->szA : ip->pszA);
      ip->a2blk(K, M, ip->alpA, A, ip->lda, rw, iw);
      return(iw);
   #else
      w = IdxAw_ip(ip, w, i, k);
      A = IdxA_ip(ip, A, i, k);
      ip->a2blk(K, M, ip->alpA, A, ip->lda, w);
      return(w);
   #endif
}

TYPE *Mjoin(PATL,ipcopyB)(ipinfo_t *ip, const TYPE *B, size_t k, size_t j,
                          TYPE *w)
{
   ATL_CUINT K = (k != ip->nfkblks) ? ip->kb : ip->kb0;
   ATL_UINT N = (j < ip->nfnblks) ? ip->nb : ip->pnb;
   if (j == ip->nfnblks + ip->npnblks - 1)
      N = ip->nF;
   #ifdef TCPLX
      TYPE *iw, *rw;
      B = IdxB_ip(ip, B, k, j);
      iw = IdxBw_ip(ip, w, k, j);
      rw = iw + ((j < ip->nfnblks) ? ip->szB : ip->pszB);
      ip->b2blk(K, N, ip->alpB, B, ip->ldb, rw, iw);
      return(iw);
   #else
      w = IdxBw_ip(ip, w, k, j);
      B = IdxB_ip(ip, B, k, j);
      ip->b2blk(K, N, ip->alpB, B, ip->ldb, w);
      return(w);
   #endif
}

void Mjoin(PATL,DoWork_amm_gMNK)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_gMNK_t *pd = pp->PD;
   ipinfo_t *ip = pd->ip;
   const TYPE *A=pd->A, *B=pd->B;
   TYPE *pB = pd->wB, *pA = pd->wA;
   TYPE *pC = pd->wC + vrank*((ip->szC)SHIFT), *C=pd->C;
   ablk2cmat_t blk2c=ip->blk2c;
   #ifdef TCPLX
      TYPE *rC=pC+ip->szC;
      const TYPE *beta=pd->beta;
   #else
      const TYPE beta=pd->beta;
      #define rC pC
   #endif
   ATL_CUINT nByBlks=pd->nByBlks, nByRows=pd->nByRows;
   const size_t nmblks=ip->nfmblks+ip->npmblks, nnblks=ip->nfnblks+ip->npnblks;
   const size_t nkblks = ip->nfkblks + 1;
   const size_t nAblks = pd->nAblks, nBblks = pd->nBblks;
   size_t ctr;
/*
 * First, copy all of B & A
 */
   while ( (ctr = ATL_DecGlobalAtomicCount(pd->asgBctr, vrank)) )
   {
      size_t j, k=nBblks-ctr;
      j = k / nkblks;
      k -= j*nkblks;
      Mjoin(PATL,ipcopyB)(ip, B, k, j, pB);
      #if ATL_CBC_STRONG
         ATL_DecGlobalAtomicCount(pd->donBctr, vrank);
      #endif
   }
   while ( (ctr = ATL_DecGlobalAtomicCount(pd->asgActr, vrank)) )
   {
      size_t i, k=nAblks-ctr;
      i = k / nkblks;
      k -= i*nkblks;
      Mjoin(PATL,ipcopyA)(ip, A, i, k, pA);
      #if ATL_CBC_STRONG
         ATL_DecGlobalAtomicCount(pd->donActr, vrank);
      #endif
   }
/*
 * We now have global A&B copied for everyone's use.  For weakly-ordered caches,
 * we sync all thr to to make sure we can all see each others' copies;
 * Strongly-ordered caches need to hang-fire until A&B copy is complete.
 */
   #if ATL_CBC_STRONG
      while (ATL_GetGlobalAtomicCount(pd->donBctr, vrank))
         ATL_thread_yield();      /* await B cpy finish */
      while (ATL_GetGlobalAtomicCount(pd->donActr, vrank))
         ATL_thread_yield();      /* await A cpy finish */
   #else
      ATL_cbc_barrier(pp->nworkers, vrank, NULL);  /* barrier & memory fence */
   #endif
/*
 * Now loop over rowpans from [0,nByRows-1], wt thread doing entire col,
 * and global A & B copy known to be already in access-major in workspace
 */
   if (nByRows)
   {
      while ( (ctr = ATL_DecGlobalAtomicCount(pd->RowCtr, vrank)) )
      {
         int i = nByRows - ctr, j;
         const TYPE *a;
         TYPE *wa, *c;

         wa = IdxAw_ip(ip, pA, i, 0);
         for (j=0; j < nnblks; j++)
         {
            TYPE *wb, *c;
            wb = IdxBw_ip(ip, pB, 0, j);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, wa, wb, rC, pC,
                                 beta, blk2c);
         }
      }
   }
/*
 * Finally, load balance end of computation by using block-level scheduling
 * for last nByBlks rows
 */
   if (nByBlks)
   {
      while (1)
      {
         ATL_UINT max, i, imax=0;
         TYPE *wa;
/*
 *       Find row with maximum remaining blocks, and work on that one
 */
         max = ATL_GetGlobalAtomicCount(pd->BlkCtrs[0], vrank);
         for (i=1; i < nByBlks; i++)
         {
            int k;
            k = ATL_GetGlobalAtomicCount(pd->BlkCtrs[i], vrank);
            if (k >= max)
            {
               max = k;
               imax = i;
            }
         }
         if (!max)   /* if no blocks are left in any of the nByBlks rowpans */
            break;   /* we are done */
/*
 *       For chosen rowpan, work on individual blks of C with other threads.
 *       This is the last nByBlks rows, so add 1st and nByRows to get glob i
 *       Both the shared A and the global B have been copied at start of alg.
 */
         i = imax + nByRows;
         wa = IdxAw_ip(ip, pA, i, 0);
         while ( (ctr = ATL_DecGlobalAtomicCount(pd->BlkCtrs[imax], vrank)) )
         {
            TYPE *wb, *c;
            int j = nnblks - ctr;
            wb = IdxBw_ip(ip, pB, 0, j);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, wa, wb, rC, pC,
                                 beta, blk2c);
         }
      }
   }
}
#ifndef TCPLX
   #undef rC
#endif

int Mjoin(PATL,tammm_gMNK)
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
   int idx;
   ipinfo_t ip;
   ATL_tamm_gMNK_t pd;
   size_t nmblks, nnblks, nkblks, szB, szA, sz;
   void *vp=NULL;
   ATL_UINT P;

   pd.beta = beta;
   idx = Mjoin(PATL,geGetAmmmIndx)(M, N, K);
   Mjoin(PATL,geComputeIPInfo)(&ip, idx, TA, TB, M, N, K, lda, ldb, ldc,
                               alpha, beta);
/*
 * Compute how many columns to handle wt block-level and colpan-lvl sheduling
 */
   nmblks = ip.nfmblks + ip.npmblks;
   nnblks = ip.nfnblks + ip.npnblks;
   nkblks = ip.nfkblks + 1;
   pd.nAblks = nmblks * nkblks;
   pd.nBblks = nkblks * nnblks;
   #if 1
      pd.nByBlks = Mmin(32,ATL_NTHREADS-1); /* huge causes too many ctrs */
      pd.nByBlks = Mmin(pd.nByBlks, nmblks);
   #elif 0
      pd.nByBlks = 0;         /* just for testing, bad load balance */
   #else
      pd.nByBlks = nmblks;    /* just for testing, bad parallel overhead */
   #endif
   pd.nByRows = nmblks - pd.nByBlks;
/*
 * Compute max parallelism, and call serial if inadequate
 */
   P = pd.nByRows + pd.nByBlks*nnblks;
   P = Mmin(P, ATL_NTHREADS);
   if (0 && P < 2)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
   #if 0
   printf("P=%d, nByRows=%d, nByBlks=%d, nfnblks=%d, nb=%d npnblks=%d, pnb=%d\n",
          (int)P, (int)pd.nByRows, (int)pd.nByBlks, (int)ip.nfnblks, ip.nb,
          (int)ip.npnblks, ip.pnb);
   #endif
   pd.ip = &ip;
   pd.A = A;
   pd.B = B;
   pd.C = C;
/*
 * Get wrkspc
 */
   szA = (ip.szA*ip.nfmblks + ip.pszA*ip.npmblks)*nkblks;
   szB = (ip.szB*ip.nfnblks + ip.pszB*ip.npnblks)*nkblks;
   sz = (ip.szC SHIFT)*P + (ip.mu<<1)*ip.nu;
   sz = ATL_MulBySize(szA + szB + sz) + 3*ATL_Cachelen;
   if (sz <= ATL_PTMAXMALLOC)
      vp = malloc(sz);
   if (!vp)
      return(1);
   pd.wA = ATL_AlignPtr(vp);
   pd.wB = pd.wA + (szA SHIFT);
   pd.wB = ATL_AlignPtr(pd.wB);
   pd.wC = pd.wB + (szB SHIFT);
   pd.wC = ATL_AlignPtr(pd.wC);

   if (pd.nByRows)
      pd.RowCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nByRows, P),
                                            pd.nByRows, 0);
   else
      pd.RowCtr = NULL;
   pd.asgActr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nAblks,P), pd.nAblks,0);
   pd.asgBctr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nBblks,P), pd.nBblks,0);
   #if ATL_CBC_STRONG
      pd.donActr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nAblks,P),
                                            pd.nAblks, 0);
      pd.donBctr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nBblks,P),
                                            pd.nBblks, 0);
   #endif
   if (pd.nByBlks)
   {
      int nat, i;
      pd.BlkCtrs = malloc(pd.nByBlks * sizeof(void*));
      ATL_assert(pd.BlkCtrs);
      nat = ATL_EstNctr(nnblks,P);
      for (i=0; i < pd.nByBlks; i++)
         pd.BlkCtrs[i] = ATL_SetGlobalAtomicCount(nat, nnblks, 0);
   }
   else
      pd.BlkCtrs = NULL;
   ATL_goParallel(P, Mjoin(PATL,DoWork_amm_gMNK), NULL, &pd, NULL);

   if (pd.nByBlks)
   {
      int i;
      for (i=0; i < pd.nByBlks; i++)
         ATL_FreeGlobalAtomicCount(pd.BlkCtrs[i]);
      free(pd.BlkCtrs);
   }
   ATL_FreeGlobalAtomicCount(pd.asgActr);
   ATL_FreeGlobalAtomicCount(pd.asgBctr);
   #if ATL_CBC_STRONG
      ATL_FreeGlobalAtomicCount(pd.donActr);
      ATL_FreeGlobalAtomicCount(pd.donBctr);
   #endif
   if (pd.nByRows)
      ATL_FreeGlobalAtomicCount(pd.RowCtr);
   free(vp);
   return(0);
}
