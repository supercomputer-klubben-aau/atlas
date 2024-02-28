#include "atlas_amm.h"
/*
 * This routine handles N <= MAXN, K & M large (left-looking shape)
 */
int Mjoin(PATL,ammm_tN)
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
   size_t szB, szC;
   TYPE *a, *b, *c;
   void *vp;
   int idx;
   ipinfo_t ip;


   Mjoin(PATL,ipgenInfo)(&ip, 2, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   szC = ip.szC;
   szB = ip.pszB*ip.npnblks + ip.szB*ip.nfnblks;
   szB = ip.szB*(ip.nfkblks+1);
   vp = malloc(ATL_MulBySize(szC+ip.szA+szB+ip.exsz) + 3*ATL_Cachelen);
   if (!vp)
      return(1);
   a = ATL_AlignPtr(vp);
   b = a + (ip.szA SHIFT);
   b = ATL_AlignPtr(b);
   c = b + (szB SHIFT);
   c = ATL_AlignPtr(c);
   Mjoin(PATL,iploopsMK)(&ip, 0, 0, A, B, C, 2, a, b,
   #ifdef TCPLX
                         c+szC, c, beta, ip.blk2c);
   #else
                         c, c, beta, ip.blk2c);
   #endif

   free(vp);
   return(0);
}
