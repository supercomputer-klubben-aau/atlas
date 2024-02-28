#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
void Mjoin(PATL,DoWork_tamm_tNK)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp = vp;
   ATL_tamm_tNK_t *pd = lp->opstruct;  /* parallel prob def struct */
   opinfo_t *rp = pd->rp;
   size_t ldc=rp->ldc, lda=rp->lda;
   ATL_CUINT rank = tp->rank, nfmblks=rp->nfmblks, npmblks=rp->npmblks;
   ATL_CUINT nmblks=nfmblks+npmblks, kb = rp->KB, K=rp->kb;
   TYPE *pB = pd->wB, *pA = pd->w + rank*pd->wsz, *pC, *C=pd->C;
   const TYPE *A=pd->A;
   #ifdef TCPLX
      TYPE *rC;
   #else
      #define rC pC
   #endif

   ATL_UINT BCOPIED=0, mb, nb, nmu, nnu;
   ATL_CUINT N = (rp->nfnblks) ? rp->nb : rp->nF;
   ablk2cmat_t blk2c = rp->blk2C;
   cm2am_t a2blk = rp->a2blk;
   int ictr;
   #ifdef TCPLX
      ammkern_t amm_b0=rp->amm_b0, amm_b1=rp->amm_b1, amm_bn=rp->amm_bn;
   #else
      ammkern_t amm=rp->amm_b0;
   #endif

   pC = pA + ((rp->szA)SHIFT);
   pC = ATL_AlignPtr(pC);
   #ifdef TCPLX
      rC = pC + rp->szC;
   #endif
/*
 * First guy here starts to copy B
 */
   if (ATL_DecAtomicCount(pd->BassgCtr))
   {
/*
 *    Copy B, which is known to be only one block in this algorithm
 */
      Mjoin(PATL,opblk)(rp, 0, 0, NULL, pd->B, NULL, NULL, NULL, pB, NULL,
                        NULL, NULL);
/*
 *    Let waiting threads know B is ready for use, memsync for weakly-ordered
 */
      #if ATL_CBC_STRONG
         BCOPIED = ATL_DecAtomicCount(pd->BdoneCtr);
         ATL_assert(BCOPIED);
      #else
         ATL_cbc_barrier(pd->P, rank, NULL);
         BCOPIED = 1;
      #endif
   }
/*
 * For first A block I work on, I must await completion of B copy
 */
   if (!BCOPIED)
   {

      ictr = ATL_DecGlobalAtomicCount(pd->MbCtr, rank);
      if (ictr)
      {
         int iblk = nmblks - ictr;
/*
 *       Copy A, then await B done ACK
 */
         Mjoin(PATL,opblk)(rp, iblk, 0, A, NULL, NULL, pA, NULL, NULL, NULL,
                           NULL, NULL);
         #if ATL_CBC_STRONG
            while (ATL_GetAtomicCount(pd->BdoneCtr))  /* await B cpy finish */
               ATL_thread_yield();
         #else
            ATL_cbc_barrier(pd->P, rank, NULL);
         #endif
/*
 *       Now multiply global B * local A, and write to my piece of C
 */
         Mjoin(PATL,opblk)(rp, iblk, 0, NULL, NULL, C, pA, pA, pB, pB, rC, pC);
      }
   }
/*
 * Now, B is ready, so just go to town on remaining blocks
 */
   while ((ictr = ATL_DecGlobalAtomicCount(pd->MbCtr, rank)))
   {
      Mjoin(PATL,opblk)(rp, nmblks-ictr, 0, A, NULL, C, pA, pA, pB, pB, rC, pC);
   }
}
#ifndef TCPLX
   #undef rC
#endif
/*
 * This routine handles the case where N <= maxNB && K <= maxKB, so B is
 * only one block.  It is particularly important for the panel factorizations
 * of both LU and QR.
 */
int Mjoin(PATL,tammm_tNK)
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
   ATL_tamm_tNK_t tnk;
   size_t nmblks;
   void *vp;

   if (N > ATL_rkAMM_LASTNB || K >= ATL_VWopgen_MAX_KB || M < ATL_rkAMM_LASTMB ||
       M < Mmin(8,ATL_NTHREADS)*ATL_rkAMM_LASTMB)
      return(1);
   tnk.P = Mjoin(PATL,tGetopinfo_tNK)(&rki, ATL_NTHREADS, TA, TB, M, N, K,
                                  lda, ldb, ldc, alpha, beta);
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
   tnk.wsz = ATL_MulBySize(rki.szA + rki.szC) + ATL_Cachelen;
   tnk.wsz = ATL_MulByCachelen(ATL_DivByCachelen(tnk.wsz + ATL_Cachelen-1));
   vp = malloc(ATL_MulBySize(rki.szB) + tnk.P*tnk.wsz + 2*ATL_Cachelen);
   if (!vp)
      return(1);
   tnk.wsz = ATL_DivBySize(tnk.wsz)SHIFT;
   tnk.wB = ATL_AlignPtr(vp);
   tnk.w = tnk.wB + (rki.szB SHIFT);
   tnk.w = ATL_AlignPtr(tnk.w);
   nmblks = rki.nfmblks + rki.npmblks;
   tnk.MbCtr = ATL_SetGlobalAtomicCount(ATL_EstNctr(nmblks, tnk.P), nmblks, 0);
   tnk.BassgCtr = ATL_SetAtomicCount(1);
   #if ATL_CBC_STRONG
      tnk.BdoneCtr = ATL_SetAtomicCount(1);
   #endif
   ATL_goparallel(tnk.P, Mjoin(PATL,DoWork_tamm_tNK), &tnk, NULL);

   #if ATL_CBC_STRONG
   ATL_FreeAtomicCount(tnk.BdoneCtr);
   #endif
   ATL_FreeAtomicCount(tnk.BassgCtr);
   ATL_FreeGlobalAtomicCount(tnk.MbCtr);
   free(vp);
   return(0);
}
