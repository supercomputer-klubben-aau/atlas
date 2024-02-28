/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

#define f77ger Mjoin(PATL,f77ger)
#define F77GER F77ger
void f77ger(const int M, const int N, const SCALAR alpha0,
            const TYPE *X, const int incX, const TYPE *Y, const int incY,
            TYPE *A, const int lda)
{
   #ifdef TCPLX
      TYPE alpha[2];
   #else
      TYPE alpha = alpha0;
   #endif
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77M=M, F77N=N, F77lda=lda, F77incx=incX, F77incy=incY;
   #else
      #define F77M    M
      #define F77N    N
      #define F77lda  lda
      #define F77incx incX
      #define F77incy incY
   #endif
   #ifdef TCPLX
      alpha[0] = *alpha0; alpha[1] = alpha0[1];
   #endif
   if (incX < 0) X -= (1-M)*incX SHIFT;
   if (incY < 0) Y -= (1-N)*incY SHIFT;

   F77GER(&F77M, &F77N, SADD alpha, X, &F77incx, Y, &F77incy, A, &F77lda);
}
