#include "atlas_misc.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
#include "atlas_cache.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
/*
 * This routine handles M <= MAXM && N <= MAXN && very long K, or
 * the inner-product GEMM form.  It appears in the GEMM-based SYRK, which
 * is important for Cholesky.  It is typically the worst-case for ATLAS,
 * since the copy of A and B are of the same order as the computation.
 * It is a minimal workspace routine.
 */
int Mjoin(PATL,ammm_IP)
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
   ipinfo_t ip;
   TYPE *a, *b, *c;
   void *vp, *ipp;
   unsigned int i, imm, mu, nu, mb, nb;

   Mjoin(PATL,ipgenInfo)(&ip, 3, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   vp = malloc(ATL_MulBySize(ip.szC+ip.szA+ip.szB+ip.exsz)
               + 3*ATL_Cachelen);
   if (!vp)
      return(1);

   a = ATL_AlignPtr(vp);
   b = a + (ip.szA SHIFT);
   b = ATL_AlignPtr(b);
   c = b + (ip.szB SHIFT);
   c = ATL_AlignPtr(c);
   Mjoin(PATL,iploopsK)(&ip, 0, 0, A, B, C, 0, a, b,
   #ifdef TCPLX
                        c+ip.szC, c, beta, ip.blk2c);
   #else
                        c, c, beta, ip.blk2c);
   #endif
   free(vp);
   return(0);
}
