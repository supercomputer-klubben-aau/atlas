/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_lvl2.h"
#include "atlas_lvl3.h"

void Mjoin(PATL,mvtk_smallN)
   (ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *X, ATL_CINT incX, const SCALAR beta, TYPE *Y, ATL_CINT incY)
/*
 * y = alpha*A*x + beta*y
 * A is MxN, so X is of length M and Y is of length M
 * NOTE: this routines is usually called for cleanup or when N is too small
 *       to allow for vector copies
 * This routine calls dot product, in the hope that the Level 1 has been
 * at least somewhat tuned to the architecture.
 */
{
#ifdef TREAL
   TYPE y0, dot;
   ATL_INT j;
   ATL_CINT BetaIsNonZero = (beta != ATL_rzero);

   for (j=0; j < N; j++, A += lda, Y += incY)
   {
      y0 = (BetaIsNonZero) ? *Y * beta : ATL_rzero;
      dot = alpha * Mjoin(PATL,dot)(M, A, 1, X, incX);
      *Y = y0 + dot;
   }
#else
   TYPE ry, iy, rd, id, tmp;
   const TYPE rbe = *beta, ibe = beta[1], ral = *alpha, ial = alpha[1];
   ATL_CINT lda2 = lda+lda, incY2 = incY+incY;
   ATL_INT j;

   if (ibe == ATL_rzero) /* real scalar */
   {
      if (rbe == ATL_rzero)     /* beta is zero */
      {
         for (j=0; j < N; j++, A += lda2, Y += incY2)
         {
            Mjoin(PATL,dotu_sub)(M, A, 1, X, incX, Y);
            rd = Y[0];
            id = Y[1];
            tmp = rd*ral - id*ial;
            id  = rd*ial + id*ral;
            *Y = tmp;
            Y[1] = id;
         }
      }
      else if (rbe == ATL_rone) /* beta is one */
      {
         for (j=0; j < N; j++, A += lda2, Y += incY2)
         {
            ry = *Y;
            iy = Y[1];
            Mjoin(PATL,dotu_sub)(M, A, 1, X, incX, Y);
            rd = Y[0];
            id = Y[1];
            tmp = rd*ral - id*ial;
            id  = rd*ial + id*ral;
            *Y = ry+tmp;
            Y[1] = iy+id;
         }
      }
      else                      /* beta is arbitrary real scalar */
      {
         for (j=0; j < N; j++, A += lda2, Y += incY2)
         {
            ry = *Y * rbe;
            iy = Y[1] * rbe;
            Mjoin(PATL,dotu_sub)(M, A, 1, X, incX, Y);
            rd = Y[0];
            id = Y[1];
            tmp = rd*ral - id*ial;
            id  = rd*ial + id*ral;
            *Y = ry+tmp;
            Y[1] = iy+id;
         }
      }
   }
   else /* beta is a complex scalar */
   {
      for (j=0; j < N; j++, A += lda2, Y += incY2)
      {
         ry = *Y;
         iy = Y[1];
         tmp = ry*rbe - iy*ibe;
         iy  = ry*ibe + iy*rbe;
         ry  = tmp;
         Mjoin(PATL,dotu_sub)(M, A, 1, X, incX, Y);
         rd = Y[0];
         id = Y[1];
         tmp = rd*ral - id*ial;
         id  = rd*ial + id*ral;
         *Y = ry+tmp;
         Y[1] = iy+id;
      }
   }
#endif
}
