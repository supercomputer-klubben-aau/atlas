#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"

/*
#undef ATL_CBC_STRONG
#define ATL_CBC_STRONG 0
*/

void Mjoin(PATL,DoWork_amm_sNK)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_sNK_t *pd = pp->PD;
   ipinfo_t *ip = pd->ip;
   const TYPE *B = pd->B, *A = pd->A;
   TYPE *pB = pd->wB, *pA = pd->w + vrank*pd->wsz, *pC, *C=pd->C;
   ablk2cmat_t blk2c=ip->blk2c;
   #ifdef TCPLX
      TYPE *rC;
      const TYPE *beta=pd->beta;
   #else
      const TYPE beta=pd->beta;
      #define rC pC
   #endif
   ATL_CUINT nByBlks=pd->nByBlks, nByRows=pd->nByRows;
   const size_t nmblks=ip->nfmblks+ip->npmblks, nnblks=ip->nfnblks+ip->npnblks;
   const size_t nkblks = ip->nfkblks + 1, szAp=pd->szAp;
   int ictr;
/*
 * First, copy last nByRows A operands
 */
   pC = pA + szAp;
   pC = ATL_AlignPtr(pC);
   #ifdef TCPLX
      rC = pC + ip->szC;
   #endif
   if (nByBlks)
   {
      while ( (ictr = ATL_DecAtomicCount(pd->begABlksCtr)) )
      {
          size_t i=nmblks-ictr, k;
          ATL_UINT nb, nnu;
          TYPE *wa = pd->wAb+(nByBlks-ictr)*szAp;
          Mjoin(PATL,iploopsK)(ip, i, 0, IdxA_ip(ip, A, i, 0), NULL, NULL,
                               3, wa, NULL, NULL, NULL, ip->alpA, NULL);
          #if ATL_CBC_STRONG
             ATL_DecAtomicCount(pd->donABlksCtr);
          #endif
      }
   }
/*
 * Now, compute 1st rowpan of C while copying global B;
 * strongly-ordered cachces uses single shared A, while weak redundant A
 */
   ictr = ATL_DecGlobalAtomicCount(pd->begBCtr, vrank);
   if (ictr)
   {
      int ACOPIED = 0;
      int j = nnblks - ictr;
      TYPE *wb, *wa;
      #if ATL_CBC_STRONG
         wa = pd->wAb + nByBlks*szAp;
         if (ATL_DecAtomicCount(pd->begACtr))
         {
             Mjoin(PATL,iploopsK)(ip, 0, 0, A, NULL, NULL, 3, wa, NULL,
                                  NULL, NULL, ip->alpA, NULL);
            ACOPIED = ATL_DecAtomicCount(pd->donACtr);
            ATL_assert(ACOPIED);
         }
/*
 *       Copy B, then await A done ACK
 */
         wb = IdxBw_ip(ip, pB, 0, j);
         Mjoin(PATL,iploopsK)(ip, 0, j, NULL, IdxB_ip(ip, B, 0, j), NULL, 3,
                              NULL, wb, NULL, NULL, ip->alpB, NULL);
         ATL_DecGlobalAtomicCount(pd->donBCtr, vrank);
         if (!ACOPIED)
         {
            while (ATL_GetAtomicCount(pd->donACtr))  /* await A cpy finish */
               ATL_thread_yield();
         }
/*
 *       Now multiply local(weak)/glob B * local A, and write to my piece of C
 */
         Mjoin(PATL,iploopsK)(ip, 0, j, NULL, NULL, IdxC_ip(ip, C, 0, j),
                              3, wa, wb, rC, pC, beta, blk2c);
      #else
         wa = pA;
/*
 *       Copy chosen block to correct place in global B workspace & do multiply
 */
         wb = IdxBw_ip(ip, pB, 0, j);
         Mjoin(PATL,iploopsK)(ip, 0, j, A, IdxB_ip(ip, B, 0, j),
                              IdxC_ip(ip, C, 0, j), 3, wa, wb, rC, pC,
                              beta, blk2c);
      #endif
/*
 *    Now loop over any remaining B blks with known-good wa
 */
      while ( (ictr = ATL_DecGlobalAtomicCount(pd->begBCtr, vrank)) )
      {
         j = nnblks - ictr;
         wb = IdxBw_ip(ip, pB, 0, j);
         #if ATL_CBC_STRONG
            Mjoin(PATL,iploopsK)(ip, 0, j, NULL, IdxB_ip(ip, B, 0, j), NULL, 3,
                                 NULL, wb, NULL, NULL, ip->alpB, NULL);
            ATL_DecGlobalAtomicCount(pd->donBCtr, vrank);
            Mjoin(PATL,iploopsK)(ip, 0, j, NULL, NULL, IdxC_ip(ip, C, 0, j), 3,
                                 wa, wb, rC, pC, beta, blk2c);
         #else
            Mjoin(PATL,iploopsK)(ip, 0, j, NULL, IdxB_ip(ip, B, 0, j),
                                 IdxC_ip(ip, C, 0, j), 3, wa, wb, rC, pC,
                                 beta, blk2c);
         #endif
      }
   }
/*
 * We now have global B copied for everyone's use.  For weakly-ordered caches,
 * we sync all thr to to make sure we can all see each others' copies;
 * Strongly-ordered caches need to hang-fire until B copy is complete.
 */
   #if ATL_CBC_STRONG
      while (ATL_GetGlobalAtomicCount(pd->donBCtr, vrank))
         ATL_thread_yield();      /* await B cpy finish */
   #else
      ATL_cbc_barrier(pp->nworkers, vrank, NULL);  /* barrier & memory fence */
   #endif
/*
 * Now loop over rowpans from [1,nByRows+1], wt thread doing entire col,
 * and global B copy known to be ready to use.
 */
   if (nByRows)
   {
      while ( (ictr = ATL_DecGlobalAtomicCount(pd->RowCtr, vrank)) )
      {
         int i = nByRows - ictr + 1, j;
         const TYPE *a;
         TYPE *c, *wb;
/*
 *       For top block, use copied B, and copy this k-panel of A
 */
         a = IdxA_ip(ip, A, i, 0);
         c = IdxC_ip(ip, C, i, 0);
         Mjoin(PATL,iploopsK)(ip, i, 0, a, NULL, c, 3, pA, pB, rC, pC,
                              beta, blk2c);
/*
 *       Remaining blocks don't need to copy A or B
 */
         for (j=1; j < nnblks; j++)
         {
            TYPE *bw, *c;
            bw = IdxBw_ip(ip, pB, 0, j);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, pA, bw, rC, pC,
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
      #if ATL_CBC_STRONG
         while (ATL_GetAtomicCount(pd->donABlksCtr))  /* await A cpy finish */
            ATL_thread_yield();
      #endif
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
         i = imax + 1 + nByRows;
         wa = pd->wAb + imax*szAp;
         while ( (ictr = ATL_DecGlobalAtomicCount(pd->BlkCtrs[imax], vrank)) )
         {
            TYPE *bw, *c;
            int j = nnblks - ictr;
            bw = IdxBw_ip(ip, pB, 0, j);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, wa, bw, rC, pC,
                                 beta, blk2c);
         }
      }
   }
}
#ifndef TCPLX
   #undef rC
#endif
/*
 * This routine handles the case where K = i*maxKB, and M is reasonably large
 * so that the it can spend most of its time giving out entire rowpanels of C,
 * and only go to block-level syncs for a few rowpans at end.  It requires
 * enough workspace to allocate common B (whole matrix), and P*(szApan+szC).
 * The need to allocate entire B means that N can't be too large either.
 * It is particularly important for recursive panel LU/QR factorization.
 */
int Mjoin(PATL,tammm_sNK)
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
   ATL_tamm_sNK_t snk;
   size_t nmblks, nnblks, nkblks, szBt, szAp, szAb, sz;
   void *vp=NULL;
   ATL_UINT P;

   snk.beta = beta;
   #if 1
   Mjoin(PATL,ipgenInfo)(&ip, 0, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   #elif 1
   idx = Mjoin(PATL,tGetParCIndx)(&ip, ATL_NTHREADS, M, N, K);
   Mjoin(PATL,geFillInIPInfo)(&ip, idx, TA, TB, M, N, K, lda, ldb, ldc,
                              alpha, beta,
                              ip.nfmblks, ip.npmblks, ip.mb, ip.pmb,
                              ip.nfnblks, ip.npnblks, ip.nb, ip.pnb);
   #else
   idx = Mjoin(PATL,geGetAmmmIndx)(M, N, K);
   Mjoin(PATL,geComputeIPInfo)(&ip, idx, TA, TB, M, N, K, lda, ldb, ldc,
                               alpha, beta);
   #endif
/*
 * Compute how many columns to handle wt block-level and colpan-lvl sheduling
 */
   nmblks = ip.nfmblks + ip.npmblks;
   nnblks = ip.nfnblks + ip.npnblks;
   nkblks = ip.nfkblks + 1;
   #if 1
      snk.nByBlks = ATL_NTHREADS - 1;
      snk.nByBlks = Mmin(snk.nByBlks, 32); /* huge causes too many Ctrs */
      snk.nByBlks = Mmin(snk.nByBlks, nmblks-1);
   #elif 0
      snk.nByBlks = 0;         /* just for testing, bad load balance */
   #else
      snk.nByBlks = nmblks-1;  /* just for testing, bad parallel overhead */
   #endif
   snk.nByRows = nmblks - snk.nByBlks - 1;
/*
 * Compute max parallelism, and call serial if inadequate
 */
   P = snk.nByRows + (snk.nByBlks+1)*nnblks;
   P = Mmin(P, ATL_NTHREADS);
   if (P < 2)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
   #if 1
printf("snk:P=%d, M=(%d,%d; %d,%d) N=(%d,%d; %d,%d)\n",
       P, (int)ip.nfmblks, (int)ip.npmblks, ip.mb, ip.pmb,
       (int)ip.nfnblks, (int)ip.npnblks, ip.nb, ip.pnb);
   printf("P=%d, nByRows=%d, nByBlks=%d, nfnblks=%d, nb=%d npnblks=%d, pnb=%d\n",
          (int)P, (int)snk.nByRows, (int)snk.nByBlks, (int)ip.nfnblks, ip.nb,
          (int)ip.npnblks, ip.pnb);
   #endif
   snk.ip = &ip;
   snk.A = A;
   snk.B = B;
   snk.C = C;
/*
 * Get wrkspc
 */
   szBt = (ip.szB*ip.nfnblks + ip.pszB*ip.npnblks)*(ip.nfkblks+1);
   szAp = (ip.nfmblks) ? ip.szA:ip.pszA;
   szAp *= nkblks;
   snk.szAp = szAp SHIFT;
   #if ATL_CBC_STRONG
      szAb = szAp*(1+snk.nByBlks);
   #else
      szAb = szAp*snk.nByBlks;
   #endif
   snk.wsz = ATL_MulBySize(szAp + ip.szC) + ATL_Cachelen;
   snk.wsz = ATL_MulByCachelen(ATL_DivByCachelen(snk.wsz + ATL_Cachelen-1));
   sz = ATL_MulBySize(szBt+szAb) + P*snk.wsz + 3*ATL_Cachelen;
   if (sz <= ATL_PTMAXMALLOC)
      vp = malloc(sz);
   if (!vp)
      return(1);
   snk.wsz = ATL_DivBySize(snk.wsz)SHIFT;
   snk.wB = ATL_AlignPtr(vp);
   snk.wAb = snk.wB + (szBt SHIFT);
   snk.wAb = ATL_AlignPtr(snk.wAb);
   snk.w = snk.wAb + (szAb SHIFT);
   snk.w = ATL_AlignPtr(snk.w);
   if (snk.nByRows)
      snk.RowCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(snk.nByRows, P),
                                            snk.nByRows, 0);
   else
      snk.RowCtr = NULL;
   snk.begBCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nnblks,P), nnblks, 0);
   if (snk.nByBlks)
      snk.begABlksCtr = ATL_SetAtomicCount(snk.nByBlks);
   else
      snk.begABlksCtr = NULL;
   #if ATL_CBC_STRONG
      snk.donBCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nnblks,P), nnblks, 0);
      if (snk.nByBlks)
         snk.donABlksCtr = ATL_SetAtomicCount(snk.nByBlks);
      else
         snk.donABlksCtr = NULL;
      snk.begACtr = ATL_SetAtomicCount(1);
      snk.donACtr = ATL_SetAtomicCount(1);
   #endif
   if (snk.nByBlks)
   {
      int nat, i;
      snk.BlkCtrs = malloc(snk.nByBlks * sizeof(void*));
      ATL_assert(snk.BlkCtrs);
      nat = ATL_EstNctr(nnblks,P);
      for (i=0; i < snk.nByBlks; i++)
         snk.BlkCtrs[i] = ATL_SetGlobalAtomicCount(nat, nnblks, 0);
   }
   else
      snk.BlkCtrs = NULL;
   ATL_goParallel(P, Mjoin(PATL,DoWork_amm_sNK), NULL, &snk, NULL);

   if (snk.nByBlks)
   {
      int i;
      for (i=0; i < snk.nByBlks; i++)
         ATL_FreeGlobalAtomicCount(snk.BlkCtrs[i]);
      free(snk.BlkCtrs);
   }
   ATL_FreeGlobalAtomicCount(snk.begBCtr);
   if (snk.nByBlks)
      ATL_FreeAtomicCount(snk.begABlksCtr);
   #if ATL_CBC_STRONG
      ATL_FreeAtomicCount(snk.begACtr);
      ATL_FreeAtomicCount(snk.donACtr);
      ATL_FreeGlobalAtomicCount(snk.donBCtr);
      if (snk.nByBlks)
         ATL_FreeAtomicCount(snk.donABlksCtr);
   #endif
   if (snk.nByRows)
      ATL_FreeGlobalAtomicCount(snk.RowCtr);
   free(vp);
   return(0);
}
