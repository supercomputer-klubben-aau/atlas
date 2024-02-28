/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_lvl2.h"
#include "atlas_lvl3.h"

void Mjoin(PATL,mvnk_smallN)
   (ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *X, ATL_CINT incX, const SCALAR beta, TYPE *Y, ATL_CINT incY)
/*
 * y = alpha*A*x + beta*y, A is MxN, len(X) = N, len(Y) = M
 */
#ifdef TREAL
{
   ATL_INT j;

   if (alpha == ATL_rzero)
   {
      if (beta == ATL_rzero)
         Mjoin(PATL,zero)(M, Y, incY);
      else if (beta != ATL_rone)
         Mjoin(PATL,scal)(M, beta, Y, incY);
      return;
   }
   if (beta == ATL_rzero)
   {
      Mjoin(PATL,cpsc)(M, *X * alpha, A, 1, Y, incY);
      A += lda;
      j = 1;
      X += incX;
   }
   else if (beta != ATL_rone)
   {
      Mjoin(PATL,axpby)(M, *X * alpha, A, 1, beta, Y, incY);
      A += lda;
      j = 1;
      X += incX;
   }
   else
      j=0;
   for (; j < N; j++, A += lda, X += incX)
      Mjoin(PATL,axpy)(M, *X * alpha, A, 1, Y, incY);
}
#else
{
   TYPE xx[2];
   const register TYPE ral=(*alpha), ial=alpha[1], rbe=(*beta), ibe=beta[1];
   register TYPE rx, ix;
   ATL_INT j;
   ATL_CINT lda2 = lda+lda, incX2 = incX+incX;

/*
 * Handle alpha=0 case where A & X don't matter
 */
   if (ral == ATL_rzero && ial == ATL_rzero)
   {
      if (ibe == ATL_rzero)
      {
         if (rbe == ATL_rzero)
            Mjoin(PATL,zero)(M, Y, incY);
         else if (rbe != ATL_rone)
            Mjoin(PATL,scal)(M, beta, Y, incY);
      }
      else
         Mjoin(PATL,scal)(M, beta, Y, incY);
      return;
   }
/*
 * Apply beta only on first write of Y
 */
   if (rbe == ATL_rzero && ibe == ATL_rzero)
   {
      rx = *X; ix = X[1];
      xx[0] = rx*ral - ix*ial;
      xx[1] = rx*ial + ix*ral;
      Mjoin(PATL,cpsc)(M, xx, A, 1, Y, incY);
      A += lda2;
      j = 1;
      X += incX2;
   }
   else if (rbe != ATL_rone || ibe != ATL_rzero)
   {
      rx = *X; ix = X[1];
      xx[0] = rx*ral - ix*ial;
      xx[1] = rx*ial + ix*ral;
      Mjoin(PATL,axpby)(M, xx, A, 1, beta, Y, incY);
      A += lda2;
      j = 1;
      X += incX2;
   }
   else
      j=0;
   for (; j < N; j++, A += lda2, X += incX2)
   {
      rx = *X;
      ix = X[1];
      xx[0] = rx*ral - ix*ial;
      xx[1] = rx*ial + ix*ral;
      Mjoin(PATL,axpy)(M, xx, A, 1, Y, incY);
   }
}
#endif
