/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009, 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl2.h"

void ATL_UGERK
   (ATL_CINT M, ATL_CINT N, const TYPE *X, const TYPE *Y, TYPE *A, ATL_CINT lda)
{
   int i, j;
   const TYPE *x;
   TYPE *A0 = A, *A1 = A + lda, *A2 = A1 + lda, *A3 = A2 + lda;
   ATL_CINT M8 = ((M-1)>>3)<<3, mr = M-M8-1;
   ATL_CINT N4 = (N>>2)<<2, incAn = (lda<<2) - M + 1;
   register TYPE m0, m1, m2, m3, x0, y0, y1, y2, y3;

   if (M8)
   {
      for (j=N4; j; j -= 4)
      {
         y0 = *Y; y1 = Y[1]; y2 = Y[2]; y3 = Y[3];
         x0 = *X; x = X + 1;
         m0 = y0 * x0; Y += 4;
         m1 = y1 * x0;
         m2 = y2 * x0;
         m3 = y3 * x0;
         for (i=M8; i; i -= 8)
         {
            x0 = *x;
            *A0 += m0; m0 = y0 * x0;
            *A1 += m1; m1 = y1 * x0;
            *A2 += m2; m2 = y2 * x0;
            *A3 += m3; m3 = y3 * x0;
            x0 = x[1];
            A0[1] += m0; m0 = y0 * x0;
            A1[1] += m1; m1 = y1 * x0;
            A2[1] += m2; m2 = y2 * x0;
            A3[1] += m3; m3 = y3 * x0;
            x0 = x[2];
            A0[2] += m0; m0 = y0 * x0;
            A1[2] += m1; m1 = y1 * x0;
            A2[2] += m2; m2 = y2 * x0;
            A3[2] += m3; m3 = y3 * x0;
            x0 = x[3];
            A0[3] += m0; m0 = y0 * x0;
            A1[3] += m1; m1 = y1 * x0;
            A2[3] += m2; m2 = y2 * x0;
            A3[3] += m3; m3 = y3 * x0;
            x0 = x[4];
            A0[4] += m0; m0 = y0 * x0;
            A1[4] += m1; m1 = y1 * x0;
            A2[4] += m2; m2 = y2 * x0;
            A3[4] += m3; m3 = y3 * x0;
            x0 = x[5];
            A0[5] += m0; m0 = y0 * x0;
            A1[5] += m1; m1 = y1 * x0;
            A2[5] += m2; m2 = y2 * x0;
            A3[5] += m3; m3 = y3 * x0;
            x0 = x[6];
            A0[6] += m0; m0 = y0 * x0;
            A1[6] += m1; m1 = y1 * x0;
            A2[6] += m2; m2 = y2 * x0;
            A3[6] += m3; m3 = y3 * x0;
            x0 = x[7]; x += 8;
            A0[7] += m0; m0 = y0 * x0; A0 += 8;
            A1[7] += m1; m1 = y1 * x0; A1 += 8;
            A2[7] += m2; m2 = y2 * x0; A2 += 8;
            A3[7] += m3; m3 = y3 * x0; A3 += 8;
         }
         if (mr)
         {
            for (i=mr; i; i--)
            {
               x0 = *x++;
               *A0++ += m0; m0 = y0 * x0;
               *A1++ += m1; m1 = y1 * x0;
               *A2++ += m2; m2 = y2 * x0;
               *A3++ += m3; m3 = y3 * x0;
            }
         }
         *A0 += m0; A0 += incAn;
         *A1 += m1; A1 += incAn;
         *A2 += m2; A2 += incAn;
         *A3 += m3; A3 += incAn;
      }
      if (N-N4)
         Mjoin(PATL,gerk_axpy)(M, N-N4, ATL_rone, X, 1, Y, 1, A0, lda);
   }
   else
      Mjoin(PATL,gerk_Mlt16)(M, N, ATL_rone, X, 1, Y, 1, A, lda);
}
