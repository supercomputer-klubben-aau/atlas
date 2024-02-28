/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"

static void my_swap
   (const int N, TYPE *X, const int incx, TYPE *Y, const int incy)
{
   const int incX = incx SHIFT, incY = incy SHIFT;
   int i;
   TYPE t0;
   for (i=N; i; i--, X += incX, Y += incY)
   {
      t0 = *X;
      *X = *Y;
      *Y = t0;
      #ifdef TCPLX
         t0 = X[1];
         X[1] = Y[1];
         Y[1] = t0;
      #endif
   }
}
void Mjoin(PATL,tstsqtran)(const int N, TYPE *A, const int lda0)
/*
 * transposes the square matrix A, easiest (and slowest) algorithm possible
 */
{
   const int lda=(lda0 SHIFT), ldap1=lda0+1;
   int i;

   for (i=1; i < N; i++)
      my_swap(N-i, A+(i SHIFT), ldap1, A+i*lda, ldap1);
}
