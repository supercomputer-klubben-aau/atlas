#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(Mjoin(ATLAS_PRE,amm),_sum.h))
/*
 * This routine called in degenerate case where all dims less than max block,
 * so we can do entire operation with one kernel call
 */
int Mjoin(PATL,ammm_1b)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   size_t szA, szB, szC;
   int i;
   int nmu, nnu, mu, nu, ku;
   void *vp;
   #ifdef TCPLX
      int KK;
      TYPE *iA, *iB, *iC, *rA, *rB, *rC, *p;
      ammkern_t amm;
   #else
      TYPE *pA, *pB, *pC, *p;
   #endif
   opinfo_t mminfo;

   if (K > ATL_rkAMM_LASTKB || M > ATL_rkAMM_LASTMB || N > ATL_rkAMM_LASTNB)
      return(1);
   if (Mjoin(PATL,GetInfo_1b_oprk)(&mminfo, TA, TB, M, N, K, lda, ldb, ldc,
                                   alpha, beta) == -1)
      return(-1);  /* can't do this problem with 1b! */
/*
 * These kernels all take runtime M/N, and do well with near-square, so
 * blindly use this kernel with nM = CEIL(M/mu)*mu, nN = CEIL(N/nu)*nu
 */
   mu = mminfo.mu;
   nu = mminfo.nu;
   ku = mminfo.ku;
   nmu = (mminfo.nfmblks) ? mminfo.nmu : mminfo.pnmu;
   nnu = (mminfo.nfnblks) ? mminfo.nnu : mminfo.pnnu;

   szA = mminfo.szA;
   szB = mminfo.szB;
   szC = mminfo.szC;
   vp = malloc(ATL_MulBySize(szC+szA+szB + (mu+mu)*nu*ku) + 3*ATL_Cachelen);
   if (!vp)
      return(1);
   #ifdef TCPLX
      iB = ATL_AlignPtr(vp);
      rB = iB + szB;
      iA = rB + szB;
      iA = ATL_AlignPtr(iA);
      rA = iA + szA;
      iC = rA + szA;
      iC = ATL_AlignPtr(iC);
      rC = iC + szC;
   #else
      pB = ATL_AlignPtr(vp);
      pA = pB + szB;
      pA = ATL_AlignPtr(pA);
      pC = pA + szA;
      pC = ATL_AlignPtr(pC);
   #endif
/*
 * Copy A & B into workspace, and pad K its if necessary
 */
   #ifdef TCPLX
      amm = mminfo.amm_b0;
      KK = mminfo.KB;

      mminfo.a2blk(K, M, mminfo.alpA, A, lda, rA, iA);
      mminfo.b2blk(K, N, mminfo.alpB, B, ldb, rB, iB);
      amm(nmu, nnu, KK, iA, iB, rC, rA, iB, iC);
      amm(nmu, nnu, KK, rA, iB, iC, rA, rB, rC);
      mminfo.amm_bn(nmu, nnu, KK, rA, rB, rC, iA, rB, iC);
      mminfo.amm_b1(nmu, nnu, KK, iA, rB, iC, iA, rB, iC);
      mminfo.blk2C(M, N, mminfo.ONE, rC, iC, beta, C, ldc);
   #else
      mminfo.a2blk(K, M, mminfo.alpA, A, lda, pA);
      mminfo.b2blk(K, N, mminfo.alpB, B, ldb, pB);
      mminfo.amm_b0(nmu, nnu, mminfo.KB, pA, pB, pC, pA, pB, pC);
      mminfo.blk2C(M, N, ATL_rone, pC, beta, C, ldc);
   #endif

   free(vp);
   return(0);
}
