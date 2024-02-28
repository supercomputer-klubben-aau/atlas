#include "atlas_trsmK.h"
/*
 * Take column major L, put it in block-major/amm form for use in trsm_amm.
 * In particular, diagonal blocks are copied and padded to match rank-4-based
 * trsm kernel, while blocks below the diagonal are copied to amm format.
 */
void Mjoin(PATL,trsmK_L2blk)
(
   enum ATLAS_DIAG Diag,
   const int mb0,   /* size of first triangular block */
   const int b,     /* size of all other blocks */
   int nmblks,      /* number of blocks of size b (not incl mb0) */
   const TYPE *L,
   size_t ldl,
   TYPE *W,
   size_t incW,
   cm2am_t a2blk    /* copy rect blks to amm storage */
)
{
   const int ldW0=(mb0+3)&(~3), ldW=(b+3)&(~3);
   const size_t incL = ldl*b, incL0 = ldl*mb0;
   ATL_UINT i;

   #ifdef DEBUG
      ATL_assert(incW >= ldW0*ldW0 && incW >= ldW*ldW);
   #endif

   Mjoin(PATL,trcpypad4L)(Diag, mb0, L, ldl, W, ldW0);
   L += mb0;
   W += incW;
   for (i=1; i < nmblks; i++)
   {
      const TYPE *l=L;
      ATL_UINT j, kb=mb0;

      a2blk(mb0, b, ATL_rnone, l, ldl, W);
      l += incL0; W += incW;

      for (j=1; j < i; j++, l += incL, W += incW)
         a2blk(b, b, ATL_rnone, l, ldl, W);

      Mjoin(PATL,trcpypad4L)(Diag, b, l, ldl, W, ldW);
      W += incW;
      L += b;
   }
}
