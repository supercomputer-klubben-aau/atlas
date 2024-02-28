/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"

int ATL_getrf(const enum CBLAS_ORDER Order, const int M, const int N,
              TYPE *A, const int lda, int *ipiv)
{
   if (Order == CblasColMajor) return(ATL_getrfC(M, N, A, lda, ipiv));
   else return(ATL_getrfR(M, N, A, lda, ipiv));
}
