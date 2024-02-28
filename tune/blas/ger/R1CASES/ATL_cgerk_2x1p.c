/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010, 2009, 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lvl2.h"
#include "atlas_prefetch.h"

void ATL_UGERK
   (ATL_CINT M, ATL_CINT N, const TYPE *X, const TYPE *Y, TYPE *A, ATL_CINT lda)
{
   ATL_CINT mr = M - ((M>>1)<<1);
   ATL_CINT incA = ((lda-M+mr)<<1), lda2=lda+lda;
   const TYPE *stY = Y + N+N, *stX = X + ((M>>1)<<2), *x;
   register TYPE ry, iy, rx, ix, ra, ia;

   if ( N > 0)
   {
      if (M > 1)
      {
         do
         {
            #ifdef Conj_
               ry = *Y;
               iy = -Y[1];
            #else
               ry = *Y;
               iy = Y[1];
            #endif
            x = X;
            do
            {
               rx = *x; ix = x[1];
               ra = *A; ia = A[1];
               ra += rx * ry; ATL_pfl1W(A+16);
               ia += rx * iy;
               ra -= ix * iy;
               ia += ix * ry;
               *A = ra;
               A[1] = ia;
               rx = x[2]; ix = x[3]; x += 4;
               ra = A[2]; ia = A[3];
               ra += rx * ry;
               ia += rx * iy;
               ra -= ix * iy;
               ia += ix * ry;
               A[2] = ra;
               A[3] = ia; A += 4;
            }
            while (x != stX);
            if (!mr) goto L1;
            else
            {
               rx = *x; ix = x[1];
               ra = *A; ia = A[1];
               ra += rx * ry;
               ra -= ix * iy;
               ia += rx * iy;
               ia += ix * ry;
               *A = ra;
               A[1] = ia;
            }
L1:         Y += 2;
            A += incA;
         }
         while (Y != stY);
      }
      else if (M == 1)
      {
         #ifdef Conj_
            rx = *X; ix = X[1];
            do
            {
               ry = *Y; iy = Y[1];
               ra = *A; ia = A[1];
               ra += rx * ry;
               ia -= rx * iy;
               ra += ix * iy;
               ia += ix * ry;
               *A = ra;
               A[1] = ia;
               Y += 2;
               A += lda2;
            }
            while (Y != stY);
         #else
            Mjoin(PATL,axpy)(N, X, Y, 1, A, lda);
         #endif
      }
   }
}
