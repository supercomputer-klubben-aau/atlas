#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
void Mjoin(PATL,DoWork_amm_tK)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_tK_t *pd = pp->PD;
   opinfo_t *op = pd->rp;
   const TYPE *B = pd->B, *A = pd->A;
   TYPE *pA = pd->wA, *pB = pd->w + vrank*pd->wsz, *pC, *C=pd->C;
   #ifdef TCPLX
      TYPE *rC;
   #else
      #define rC pC
   #endif
   ATL_CUINT nByBlks=pd->nByBlks, nByCols=pd->nByCols;
   const size_t nmblks=op->nfmblks+op->npmblks, nnblks=op->nfnblks+op->npnblks;
   int ictr, ict2;
/*
 * First, copy last nByCols B operands
 */
   pC = pB + ((op->szB)SHIFT);
   pC = ATL_AlignPtr(pC);
   #ifdef TCPLX
      rC = pC + op->szC;
   #endif
   if (nByBlks)
   {
      ATL_UINT szB=((op->szB)SHIFT);
      while ( (ictr = ATL_DecAtomicCount(pd->begBBlksCtr)) )
      {
          Mjoin(PATL,opblk)(op, 0, nnblks-ictr, NULL, B, NULL, NULL, NULL,
                            pd->wBb+(nByBlks-ictr)*szB, NULL, NULL, NULL);
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
      TYPE *b, *a, *an;
      #if ATL_CBC_STRONG
         b = pd->wBb + (nByBlks)*((op->szB)SHIFT);
         if (ATL_DecAtomicCount(pd->begBCtr))
         {
            Mjoin(PATL,opblk)(op, 0, 0, NULL, B, NULL, NULL, NULL, b, NULL,
                              NULL, NULL);
            BCOPIED = ATL_DecAtomicCount(pd->donBCtr);
            ATL_assert(BCOPIED);
         }
/*
 *       Copy A, then await B done ACK
 */
         a = IdxAw_rkK(op, pA, i);
         Mjoin(PATL,opblk)(op, i, 0, A, NULL, NULL, a, NULL, NULL, NULL,
                           NULL, NULL);
         ATL_DecGlobalAtomicCount(pd->donACtr, vrank);
         if (!BCOPIED)
         {
            while (ATL_GetAtomicCount(pd->donBCtr))  /* await B cpy finish */
               ATL_thread_yield();
         }
         ict2 = ATL_DecGlobalAtomicCount(pd->begACtr, vrank);
         if (ict2)
            an = IdxAw_rkK(op, pA, nmblks-ict2);
         else
            an = pB;
/*
 *       Now multiply local(weak)/glob B * local A, and write to my piece of C
 */
         Mjoin(PATL,opblk)(op, i, 0, NULL, NULL, C, a, an, b, b, rC, pC);
      #else
         b = pd->w + vrank*pd->wsz;
/*
 *       All threads make local copy of B to avoid costly sync
 */
         Mjoin(PATL,opblk)(op, 0, 0, NULL, B, NULL, NULL, NULL, b, NULL,
                           NULL, NULL);
/*
 *       Copy chosen block to correct place in global A workspace & do multiply
 */
         a = IdxAw_rkK(op, pA, i);
         ict2 = ATL_DecGlobalAtomicCount(pd->begACtr, vrank);
         if (ict2)
            an = IdxAw_rkK(op, pA, nmblks-ict2);
         else
            an = pB;
         Mjoin(PATL,opblk)(op, i, 0, A, NULL, C, a, an, b, b, rC, pC);
      #endif
/*
 *    Now loop over any remaining A blks with known-good b
 */
      while (ict2)
      {
         ictr = ict2;
         i = nmblks - ictr;
         a = an;
         #if ATL_CBC_STRONG
            Mjoin(PATL,opblk)(op, i, 0, A, NULL, NULL, a, NULL, NULL, NULL,
                              NULL, NULL);
            ATL_DecGlobalAtomicCount(pd->donACtr, vrank);
            ict2 = ATL_DecGlobalAtomicCount(pd->begACtr, vrank);
            if (an)
               an = IdxAw_rkK(op, pA, nmblks-ict2);
            else
               an = pB;
            Mjoin(PATL,opblk)(op, i, 0, NULL, NULL, C, a, an, b, b, rC, pC);
         #else
            ict2 = ATL_DecGlobalAtomicCount(pd->begACtr, vrank);
            if (an)
               an = IdxAw_rkK(op, pA, nmblks-ict2);
            else
               an = pB;
            Mjoin(PATL,opblk)(op, i, 0, A, NULL, C, a, an, b, b, rC, pC);
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
      ATL_CUINT nfnblks = op->nfnblks;
      while ( (ictr = ATL_DecGlobalAtomicCount(pd->ColCtr, vrank)) )
      {
         int i, j = nByCols - ictr + 1;
         TYPE *an;
/*
 *       For top block, use copied A, and copy this colblk of B
 */
         an = IdxAw_rkK(op, pA, 1);
         Mjoin(PATL,opblk)(op, 0, j, NULL, B, C, pA, an, pB, pB, rC, pC);
/*
 *       Remaining blocks don't need to copy A or B
 */
         for (i=1; i < nmblks; i++)
         {
            TYPE *a=an;
            an = IdxAw_rkK(op, pA, i+1);
            Mjoin(PATL,opblk)(op, i, j, NULL, NULL, C, a, an, pB, pB, rC, pC);
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
         wb = pd->wBb + jmax*((op->szB)SHIFT);
         ictr = ATL_DecGlobalAtomicCount(pd->BlkCtrs[jmax], vrank);
         if (ictr)
         {
            TYPE *a, *an;
            int i = nmblks - ictr;
            ict2 = ATL_DecGlobalAtomicCount(pd->BlkCtrs[jmax], vrank);
            a  = IdxAw_rkK(op, pA, i);
            an = IdxAw_rkK(op, pA, nmblks-ict2);
            Mjoin(PATL,opblk)(op, i, j, NULL, NULL, C, a, an, wb, wb, rC, pC);
            while (ict2)
            {
               i = nmblks-ict2;
               a = an;
               ict2 = ATL_DecGlobalAtomicCount(pd->BlkCtrs[jmax], vrank);
               an = IdxAw_rkK(op, pA, nmblks-ict2);
               Mjoin(PATL,opblk)(op, i, j, NULL, NULL, C, a, an, wb, wb, rC,pC);
            }
         }
      }
   }
}
#ifndef TCPLX
   #undef rC
#endif

/*
 * This routine handles the case where K <= maxKB, and N is reasonably large
 * so that the it can spend most of its time giving out entire colpans of C,
 * and only go to block-level syncs for a few colpans at end.  It requires
 * enouch workspace to allocate common A (whole matrix), and P*(szB+szC).
 * It is particularly important for the statically blocked right-looking LU
 * or QR factorizations.
 */
int Mjoin(PATL,tammm_tK)
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
   opinfo_t rki;
   ATL_tamm_tK_t tnk;
   size_t nmblks, nnblks, szAp, szBb, sz;
   void *vp=NULL;

   Mjoin(PATL,opinfo)(&rki, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   {
      sz = rki.nfmblks * rki.nfnblks;
      tnk.P = Mmin(ATL_NTHREADS, sz);
   }
   if (tnk.P < 2)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
/*   printf("P=%d\n", tnk.P); */
   tnk.rp = &rki;
   tnk.A = A;
   tnk.B = B;
   tnk.C = C;
/*
 * Compute how many columns to handle wt block-level and colpan-lvl sheduling
 */
   nmblks = rki.nfmblks + rki.npmblks;
   nnblks = rki.nfnblks + rki.npnblks;
   #if 1
   tnk.nByBlks = tnk.P - 1;
   tnk.nByBlks = Mmin(tnk.nByBlks, 32); /* allowing huge causes too many Ctrs */
   tnk.nByBlks = Mmin(tnk.nByBlks, nnblks-1);
   #elif 0
      tnk.nByBlks = 0;         /* just for testing, bad load balance */
   #else
      tnk.nByBlks = nnblks-1;  /* just for testing, bad parallel overhead */
   #endif
   tnk.nByCols = nnblks - tnk.nByBlks - 1;
/*
 * Get wrkspc
 */
   szAp = rki.szA*rki.nfmblks + rki.pszA*rki.npmblks;
   #if ATL_CBC_STRONG
      szBb = rki.szB*(1+tnk.nByBlks);
   #else
      szBb = rki.szB*tnk.nByBlks;
   #endif
   tnk.wsz = ATL_MulBySize(rki.szB + rki.szC) + ATL_Cachelen;
   tnk.wsz = ATL_MulByCachelen(ATL_DivByCachelen(tnk.wsz + ATL_Cachelen-1));
   sz = ATL_MulBySize(szAp+szBb+rki.exsz) + tnk.P*tnk.wsz + 3*ATL_Cachelen;
   if (sz <= ATL_PTMAXMALLOC)
      vp = malloc(sz);
   if (!vp)
      return(1);
   tnk.wsz = ATL_DivBySize(tnk.wsz)SHIFT;
   tnk.wA = ATL_AlignPtr(vp);
   tnk.wBb = tnk.wA + (szAp SHIFT);
   tnk.wBb = ATL_AlignPtr(tnk.wBb);
   tnk.w = tnk.wBb + (szBb SHIFT);
   tnk.w = ATL_AlignPtr(tnk.w);
   if (tnk.nByCols)
      tnk.ColCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(tnk.nByCols, tnk.P),
                                            tnk.nByCols, 0);
   else
      tnk.ColCtr = NULL;
   tnk.begACtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks,tnk.P), nmblks, 0);
   if (tnk.nByBlks)
      tnk.begBBlksCtr = ATL_SetAtomicCount(tnk.nByBlks);
   else
      tnk.begBBlksCtr = NULL;
   #if ATL_CBC_STRONG
      tnk.donACtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks,tnk.P),
                                             nmblks, 0);
      if (tnk.nByBlks)
         tnk.donBBlksCtr = ATL_SetAtomicCount(tnk.nByBlks);
      else
         tnk.donBBlksCtr = NULL;
      tnk.begBCtr = ATL_SetAtomicCount(1);
      tnk.donBCtr = ATL_SetAtomicCount(1);
   #endif
   if (tnk.nByBlks)
   {
      int nat, i;
      tnk.BlkCtrs = malloc(tnk.nByBlks * sizeof(void*));
      ATL_assert(tnk.BlkCtrs);
      nat = ATL_EstNctr(nmblks,tnk.P);
      for (i=0; i < tnk.nByBlks; i++)
         tnk.BlkCtrs[i] = ATL_SetGlobalAtomicCount(nat, nmblks, 0);
   }
   else
      tnk.BlkCtrs = NULL;
   ATL_goParallel(tnk.P, Mjoin(PATL,DoWork_amm_tK), NULL, &tnk, NULL);

   if (tnk.nByBlks)
   {
      int i;
      for (i=0; i < tnk.nByBlks; i++)
         ATL_FreeGlobalAtomicCount(tnk.BlkCtrs[i]);
      free(tnk.BlkCtrs);
   }
   ATL_FreeGlobalAtomicCount(tnk.begACtr);
   if (tnk.nByBlks)
      ATL_FreeAtomicCount(tnk.begBBlksCtr);
   #if ATL_CBC_STRONG
      ATL_FreeAtomicCount(tnk.begBCtr);
      ATL_FreeAtomicCount(tnk.donBCtr);
      ATL_FreeGlobalAtomicCount(tnk.donACtr);
      if (tnk.nByBlks)
         ATL_FreeAtomicCount(tnk.donBBlksCtr);
   #endif
   if (tnk.nByCols)
      ATL_FreeGlobalAtomicCount(tnk.ColCtr);
   free(vp);
   return(0);
}
