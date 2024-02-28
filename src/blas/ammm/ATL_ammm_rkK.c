#include "atlas_amm.h"
/*
 * This routine called when 2 < K <= MAXK
 */
int Mjoin(PATL,ammm_rkK)
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
   void *vp;
   TYPE *a, *b, *iC, *rC;
   opinfo_t op;
   size_t nnblks;
   #ifdef TCPLX
      const size_t ldc2=ldc+ldc;
   #else
      #define ldc2 ldc
   #endif
void Mjoin(PATL,opinfo)(opinfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);


   if (K < 3)
   {
      Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(0);
   }
   Mjoin(PATL,opinfo)(&op, TA, TB, M, N, K, lda, ldb, ldc, alpha, beta);
   nnblks = op.nfnblks + op.npnblks;
   #if 0
   printf("MB=(%d,%d), mb=(%d,%d); NB=(%d,%d), nb=(%d,%d)\n",
          op.mb, (int)op.nfmblks, op.pmb, (int)op.npmblks,
          op.nb, (int)op.nfnblks, op.pnb, (int)op.npnblks);
   #endif
/*
 *  Usual case, compute column-panel of C at a time, need workspace of
 *  1 block of B and C, and 1 col-panel of A.
 */
   if (op.nfmblks+op.npmblks != 1)
   {
      size_t sz, szA, j;
      if (nnblks > 1)
         szA = op.szA*op.nfmblks + op.pszA*op.npmblks;
      else
         szA = op.szA;
      sz = szA + op.szC + op.szB + (op.mu<<1)*op.nu;
      vp = malloc(ATL_MulBySize(sz) + 3*ATL_Cachelen);
      if (!vp)
         return(1);
      a = ATL_AlignPtr(vp);
      b = a + (szA SHIFT);
      b = ATL_AlignPtr(b);
      iC = b + (op.szB SHIFT);
      iC = ATL_AlignPtr(iC);
      #ifdef TCPLX
         rC = iC + op.szC;
      #else
         rC = iC;
      #endif

      if (nnblks == 1)
         Mjoin(PATL,oploopsM)(&op, 0, 0, A, B, C, 0, a, b, rC, iC);
      else
      {
         for (j=0; j < nnblks; j++)
         {
            Mjoin(PATL,oploopsM)(&op, 0, j, A, B, C, 1, a, b, rC, iC);
            A = NULL;
         }
      }
   }
   else  /* only 1 block of A, traverse C row-wise */
   {
      size_t sz;
      sz = op.szA + op.szB + op.szC + op.exsz;
      vp = malloc(ATL_MulBySize(sz) + 3*ATL_Cachelen);
      if (!vp)
         return(1);
      a = ATL_AlignPtr(vp);
      b = a + (op.szA SHIFT);
      b = ATL_AlignPtr(b);
      iC = b + (op.szB SHIFT);
      iC = ATL_AlignPtr(iC);
      #ifdef TCPLX
         rC = iC + op.szC;
      #else
         rC = iC;
      #endif
      if (nnblks > 1)
         Mjoin(PATL,oploopsN)(&op, 0, 0, A, B, C, 0, a, b, rC, iC);
      else
         Mjoin(PATL,opblk)(&op, 0, 0, A, B, C, a, a, b, b, rC, iC);
   }
   free(vp);
   return(0);
}
#ifdef ldc2
   #undef ldc2
#endif
