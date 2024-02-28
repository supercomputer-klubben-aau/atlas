#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"

void Mjoin(PATL,DoWork_amm_sMK)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_sMK_t *pd = pp->PD;
   ipinfo_t *ip = pd->ip;
   const TYPE *B = pd->B, *A = pd->A;
   TYPE *pA = pd->wA, *pB = pd->w + vrank*pd->wsz, *pC, *C=pd->C;
   ablk2cmat_t blk2c=ip->blk2c;
   #ifdef TCPLX
      TYPE *rC;
      const TYPE *beta=pd->beta;
   #else
      const TYPE beta=pd->beta;
      #define rC pC
   #endif
   ATL_CUINT nByBlks=pd->nByBlks, nByCols=pd->nByCols;
   const size_t nmblks=ip->nfmblks+ip->npmblks, nnblks=ip->nfnblks+ip->npnblks;
   const size_t nkblks = ip->nfkblks + 1, szBp=pd->szBp;
   int ictr;
/*
 * First, copy last nByCols B operands
 */
   pC = pB + szBp;
   pC = ATL_AlignPtr(pC);
   #ifdef TCPLX
      rC = pC + ip->szC;
   #endif
   if (nByBlks)
   {
      while ( (ictr = ATL_DecAtomicCount(pd->begBBlksCtr)) )
      {
          size_t j=nnblks-ictr, k;
          ATL_UINT nb, nnu;
          TYPE *b = pd->wBb+(nByBlks-ictr)*szBp;
          Mjoin(PATL,iploopsK)(ip, 0, j, NULL, IdxB_ip(ip, B, 0, j), NULL,
                               2, NULL, b, NULL, NULL, ip->alpB, NULL);
          #if ATL_CBC_STRONG
             ATL_DecAtomicCount(pd->donBBlksCtr);
          #endif
      }
   }
/*
 * Now, compute 1st colpan of C while copying global A;
 * strongly-ordered cachces uses shared B, while weak duplicate B
 */
   ictr = ATL_DecGlobalAtomicCount(pd->begACtr, vrank);
   if (ictr)
   {
      int BCOPIED = 0;
      int i = nmblks - ictr;
      TYPE *b, *a;
      #if ATL_CBC_STRONG
         b = pd->wBb + nByBlks*szBp;
         if (ATL_DecAtomicCount(pd->begBCtr))
         {
             Mjoin(PATL,iploopsK)(ip, 0, 0, NULL, B, NULL, 2, NULL, b,
                                  NULL, NULL, ip->alpB, NULL);
            BCOPIED = ATL_DecAtomicCount(pd->donBCtr);
            ATL_assert(BCOPIED);
         }
/*
 *       Copy A, then await B done ACK
 */
         a = IdxAw_ip(ip, pA, i, 0);
         Mjoin(PATL,iploopsK)(ip, i, 0, IdxA_ip(ip, A, i, 0), NULL, NULL, 3,
                              a, NULL, NULL, NULL, ip->alpA, NULL);
         ATL_DecGlobalAtomicCount(pd->donACtr, vrank);
         if (!BCOPIED)
         {
            while (ATL_GetAtomicCount(pd->donBCtr))  /* await B cpy finish */
               ATL_thread_yield();
         }
/*
 *       Now multiply local(weak)/glob B * local A, and write to my piece of C
 */
         Mjoin(PATL,iploopsK)(ip, i, 0, NULL, NULL, IdxC_ip(ip, C, i, 0),
                              3, a, b, rC, pC, beta, blk2c);
      #else
         b = pB;
/*
 *       Copy chosen block to correct place in global A workspace & do multiply
 */
         a = IdxAw_ip(ip, pA, i, 0);
         Mjoin(PATL,iploopsK)(ip, i, 0, IdxA_ip(ip, A, i, 0), B,
                              IdxC_ip(ip, C, i, 0), 3, a, b, rC, pC,
                              beta, blk2c);
      #endif
/*
 *    Now loop over any remaining A blks with known-good b
 */
      while ( (ictr = ATL_DecGlobalAtomicCount(pd->begACtr, vrank)) )
      {
         i = nmblks - ictr;
         a = IdxAw_ip(ip, pA, i, 0);
         #if ATL_CBC_STRONG
            Mjoin(PATL,iploopsK)(ip, i, 0, IdxA_ip(ip, A, i, 0), NULL, NULL, 3,
                                 a, NULL, NULL, NULL, ip->alpA, NULL);
            ATL_DecGlobalAtomicCount(pd->donACtr, vrank);
            Mjoin(PATL,iploopsK)(ip, i, 0, NULL, NULL, IdxC_ip(ip, C, i, 0), 3,
                                 a, b, rC, pC, beta, blk2c);
         #else
            Mjoin(PATL,iploopsK)(ip, i, 0, IdxA_ip(ip, A, i, 0), NULL,
                                 IdxC_ip(ip, C, i, 0), 3, a, b, rC, pC,
                                 beta, blk2c);
         #endif
      }
   }
/*
 * We now have global A copied for everyone's use.  For weakly-ordered caches,
 * we sync all thr to to make sure we can all see each others' copies;
 * Strongly-ordered caches need to hang-fire until A copy is complete.
 */
   #if ATL_CBC_STRONG
      while (ATL_GetGlobalAtomicCount(pd->donACtr, vrank))
         ATL_thread_yield();      /* await A cpy finish */
   #else
      ATL_cbc_barrier(pp->nworkers, vrank, NULL);  /* barrier & memory fence */
   #endif
/*
 * Now loop over colpans from [1,nByCols+1], wt thread doing entire col,
 * and global A copy known to be ready to use.
 */
   if (nByCols)
   {
      ATL_CUINT nfnblks = ip->nfnblks;
      while ( (ictr = ATL_DecGlobalAtomicCount(pd->ColCtr, vrank)) )
      {
         int i, j = nByCols - ictr + 1;
         const TYPE *b;
         TYPE *c, *wa;
/*
 *       For top block, use copied A, and copy this k-panel of B
 */
         b = IdxB_ip(ip, B, 0, j);
         c = IdxC_ip(ip, C, 0, j);
         Mjoin(PATL,iploopsK)(ip, 0, j, NULL, b, c, 3, pA, pB, rC, pC,
                              beta, blk2c);
/*
 *       Remaining blocks don't need to copy A or B
 */
         for (i=1; i < nmblks; i++)
         {
            TYPE *aw, *c;
            aw = IdxAw_ip(ip, pA, i, 0);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, aw, pB, rC, pC,
                                 beta, blk2c);
         }
      }
   }
/*
 * Finally, load balance end of computation by using block-level scheduling
 * for last nByBlks cols
 *
 * Currently, I have this designed to take any block in the last nByBlks cols.
 * However, this would result in a given thread moving pB unnecessarily.
 * A better idea is probably to nByBlks-len array of nmblks ctrs, and guys
 * can start at col+rank to spread around usage.  This will allow us to reuse
 * a given pB maximally.
 */
   if (nByBlks)
   {
      #if ATL_CBC_STRONG
         while (ATL_GetAtomicCount(pd->donBBlksCtr))  /* await B cpy finish */
            ATL_thread_yield();
      #endif
      while (1)
      {
         ATL_UINT max, j, jmax=0;
         TYPE *wb;
/*
 *       Find column with maximum remaining blocks, and work on that one
 */
         max = ATL_GetGlobalAtomicCount(pd->BlkCtrs[0], vrank);
         for (j=1; j < nByBlks; j++)
         {
            int k;
            k = ATL_GetGlobalAtomicCount(pd->BlkCtrs[j], vrank);
            if (k >= max)
            {
               max = k;
               jmax = j;
            }
         }
         if (!max)   /* if no blocks are left in any of the nByBlks colpans */
            break;   /* we are done */
/*
 *       For chosen colpan, work on individual blks of C with other threads.
 *       This is the last nByBlks columns, so add 1st and nByCols to get glob j
 *       Both the shared B and the global A have been copied at start of alg.
 */
         j = jmax + 1 + nByCols;
         wb = pd->wBb + jmax*szBp;
         while ( (ictr = ATL_DecGlobalAtomicCount(pd->BlkCtrs[jmax], vrank)) )
         {
            TYPE *aw, *c;
            int i = nmblks - ictr;
            aw = IdxAw_ip(ip, pA, i, 0);
            c = IdxC_ip(ip, C, i, j);
            Mjoin(PATL,iploopsK)(ip, i, j, NULL, NULL, c, 3, aw, wb, rC, pC,
                                 beta, blk2c);
         }
      }
   }
}
#ifndef TCPLX
   #undef rC
#endif
/*
 * This routine handles the case where K = i*maxKB, and N is reasonably large
 * so that the it can spend most of its time giving out entire colpans of C,
 * and only go to block-level syncs for a few colpans at end.  It requires
 * enough workspace to allocate common A (whole matrix), and P*(szB+szC).
 * The need to allocate entire A means that M can't be too large either.
 * It is particularly important for the statically blocked right-looking LU
 * or QR factorizations.
 */
int Mjoin(PATL,tammm_sMK)
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
   ATL_tamm_sMK_t smk;
   size_t nmblks, nnblks, nkblks, szAt, szBp, szBb, sz;
   void *vp=NULL;
   ATL_UINT P;

   smk.beta = beta;
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
      smk.nByBlks = ATL_NTHREADS - 1;
      smk.nByBlks = Mmin(smk.nByBlks, 32); /* huge causes too many Ctrs */
      smk.nByBlks = Mmin(smk.nByBlks, nnblks-1);
   #elif 1
      smk.nByBlks = 0;         /* just for testing, bad load balance */
   #else
      smk.nByBlks = nnblks-1;  /* just for testing, bad parallel overhead */
   #endif
   smk.nByCols = nnblks - smk.nByBlks - 1;
/*
 * Compute max parallelism, and call serial if inadequate
 */
   P = smk.nByCols + (smk.nByBlks+1)*nmblks;
   P = Mmin(P, ATL_NTHREADS);
   if (P < 2)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
#if 0
printf("smk:P=%d, M=(%d,%d; %d,%d) N=(%d,%d; %d,%d)\n",
       P, (int)ip.nfmblks, (int)ip.npmblks, ip.mb, ip.pmb,
       (int)ip.nfnblks, (int)ip.npnblks, ip.nb, ip.pnb);
#endif
   smk.ip = &ip;
   smk.A = A;
   smk.B = B;
   smk.C = C;
/*
 * Get wrkspc
 */
   szAt = (ip.szA*ip.nfmblks + ip.pszA*ip.npmblks)*(ip.nfkblks+1);
   szBp = (ip.nfnblks) ? ip.szB:ip.pszB;
   szBp *= nkblks;
   smk.szBp = szBp SHIFT;
   #if ATL_CBC_STRONG
      szBb = szBp*(1+smk.nByBlks);
   #else
      szBb = szBp*smk.nByBlks;
   #endif
   smk.wsz = ATL_MulBySize(szBp + ip.szC) + ATL_Cachelen;
   smk.wsz = ATL_MulByCachelen(ATL_DivByCachelen(smk.wsz + ATL_Cachelen-1));
   sz = ATL_MulBySize(szAt+szBb) + P*smk.wsz + 3*ATL_Cachelen;
   if (sz <= ATL_PTMAXMALLOC)
      vp = malloc(sz);
   if (!vp)
      return(1);
   smk.wsz = ATL_DivBySize(smk.wsz)SHIFT;
   smk.wA = ATL_AlignPtr(vp);
   smk.wBb = smk.wA + (szAt SHIFT);
   smk.wBb = ATL_AlignPtr(smk.wBb);
   smk.w = smk.wBb + (szBb SHIFT);
   smk.w = ATL_AlignPtr(smk.w);
   if (smk.nByCols)
      smk.ColCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(smk.nByCols, P),
                                            smk.nByCols, 0);
   else
      smk.ColCtr = NULL;
   smk.begACtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks,P), nmblks, 0);
   if (smk.nByBlks)
      smk.begBBlksCtr = ATL_SetAtomicCount(smk.nByBlks);
   else
      smk.begBBlksCtr = NULL;
   #if ATL_CBC_STRONG
      smk.donACtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks,P), nmblks, 0);
      if (smk.nByBlks)
         smk.donBBlksCtr = ATL_SetAtomicCount(smk.nByBlks);
      else
         smk.donBBlksCtr = NULL;
      smk.begBCtr = ATL_SetAtomicCount(1);
      smk.donBCtr = ATL_SetAtomicCount(1);
   #endif
   if (smk.nByBlks)
   {
      int nat, i;
      smk.BlkCtrs = malloc(smk.nByBlks * sizeof(void*));
      ATL_assert(smk.BlkCtrs);
      nat = ATL_EstNctr(nmblks,P);
      for (i=0; i < smk.nByBlks; i++)
         smk.BlkCtrs[i] = ATL_SetGlobalAtomicCount(nat, nmblks, 0);
   }
   else
      smk.BlkCtrs = NULL;
   ATL_goParallel(P, Mjoin(PATL,DoWork_amm_sMK), NULL, &smk, NULL);

   if (smk.nByBlks)
   {
      int i;
      for (i=0; i < smk.nByBlks; i++)
         ATL_FreeGlobalAtomicCount(smk.BlkCtrs[i]);
      free(smk.BlkCtrs);
   }
   ATL_FreeGlobalAtomicCount(smk.begACtr);
   if (smk.nByBlks)
      ATL_FreeAtomicCount(smk.begBBlksCtr);
   #if ATL_CBC_STRONG
      ATL_FreeAtomicCount(smk.begBCtr);
      ATL_FreeAtomicCount(smk.donBCtr);
      ATL_FreeGlobalAtomicCount(smk.donACtr);
      if (smk.nByBlks)
         ATL_FreeAtomicCount(smk.donBBlksCtr);
   #endif
   if (smk.nByCols)
      ATL_FreeGlobalAtomicCount(smk.ColCtr);
   free(vp);
   return(0);
}
