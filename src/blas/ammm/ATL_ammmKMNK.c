#include "atlas_misc.h"
#include "atlas_amm.h"

int Mjoin(PATL,ammmKMNK)
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
   TYPE *a, *b, *c;
   void *vp=NULL;
   size_t szA, szB, nfkblks, nkblksP;
   #ifdef TCPLX
      TYPE *rC;
      const TYPE *bet=beta;
   #else
      TYPE bet=beta;
      #define rC c
   #endif
   int i, kb, extra;
   ipinfo_t ip;

   Mjoin(PATL,ipgenInfo)(&ip, 0, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
/*
 * Figure out a partition of K that should allow malloc to succeed.
 * In this routine, we need to allocate space for a block of C,
 * a row-panel of A, and the full matrix B.  If we can't get that
 * much space, we partition K into Kp, and reduce Kp until the malloc
 * succeeds.  If Kp drops below 3*LASTKB, we give up and return non-zero,
 * to indicate that recursion should be used to cut all dims until a malloc
 * can succeed.
 */
   kb = ip.kb;
   nkblksP = nfkblks = ip.nfkblks;
/*
 * Set up extra to have all sizes that are independent of K.
 * sz[A,B] will be multiplied by number of K blks to get final answer.
 */
   extra = ATL_MulBySize(ip.szC+ip.exsz) + 3 * ATL_Cachelen;
   szB = ip.szB*ip.nfnblks + ip.pszB*ip.npnblks;
   szA = ip.szA;
   while(1)
   {
      size_t sz;
      sz = ATL_MulBySize((szA+szB)*(nkblksP+1)) + extra;
      if (sz <= ATL_MaxMalloc)
         vp = malloc(sz);
      if (vp)
         break;
      nkblksP = (nkblksP>>1);
      if (nkblksP < 3)
         return(1);
   }
   szA *= (nkblksP+1);
   szB *= (nkblksP+1);
   a = ATL_AlignPtr(vp);
   b = a + (szA SHIFT);
   b = ATL_AlignPtr(b);
   c = b + (szB SHIFT);
   c = ATL_AlignPtr(c);
   #ifdef TCPLX
      rC=c+ip.szC;
   #endif

   if (nkblksP == nfkblks)
      Mjoin(PATL,iploopsMNK)(&ip, 0, 0, A, B, C, 1|2|8, a, b,
                             rC, c, beta, ip.blk2c);
   else   /* have partitioned K into areas of at most (nkblksP+1)*kb */
   {
      const size_t incAk=ip.incAk*(nkblksP+1);
      const size_t incBk=ip.incBk*(nkblksP+1), Kp=(nkblksP+1)*kb;
      size_t k;
      ablk2cmat_t blk2c = ip.blk2c;

      for (k=Kp; k < K; k += Kp)
      {
         Mjoin(PATL,ipgenInfo)(&ip,0,TA,TB, M, N, K, lda, ldb, ldc, alpha, bet);
         Mjoin(PATL,iploopsMNK)(&ip, 0, 0, A, B, C, 1|2|8, a, b,
                                rC, c, bet, ip.blk2c);
         #ifdef TCPLX
            bet = ip.ONE;
         #else
            bet = ATL_rone;
         #endif
         A += incAk;
         B += incBk;
      }
      k = K - k + Kp;
      Mjoin(PATL,ipgenInfo)(&ip, 0, TA, TB, M, N, K, lda, ldb, ldc, alpha, bet);
      Mjoin(PATL,iploopsMNK)(&ip, 0, 0, A, B, C, 1|2|8, a, b,
                             rC, c, bet, ip.blk2c);
   }

   free(vp);
   return(0);
}
#ifndef TCPLX
   #undef rC
#endif
