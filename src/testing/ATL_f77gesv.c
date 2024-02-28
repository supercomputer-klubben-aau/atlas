/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2007 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GESV Mjoin(PRE,gesv)
#elif defined (UpCase)
   #define F77GESV Mjoin(PREU,GESV)
#elif defined (Add_) || defined(Add__)
   #define F77GESV Mjoin(PRE,gesv_)
#endif
#define f77gesv Mjoin(PATL,f77gesv)

int f77gesv(const int N, const int NRHS, TYPE *A, const int lda,
            int *ipiv, TYPE *B, const int ldb)
{
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77N=N, F77lda=lda, F77ldb=ldb, F77NRHS=NRHS;
      F77_INTEGER info;
      F77_INTEGER *F77ipiv;
   #else
      int info;
      #define F77N N
      #define F77NRHS NRHS
      #define F77lda lda
      #define F77ldb ldb
      #define F77ipiv ipiv
   #endif
   int i;
   #ifdef ATL_FunkyInts
      F77ipiv = malloc(N*sizeof(F77_INTEGER));
      ATL_assert(F77ipiv);
   #endif
   F77GESV(&F77N, &F77NRHS, A, &F77lda, F77ipiv, B, &F77ldb, &info);
   #ifdef ATL_FunkyInts
      for (i=0; i < N; i++) ipiv[i] = F77ipiv[i] - 1;
      free(F77ipiv);
   #else
      for (i=0; i < N; i++) ipiv[i]--;
   #endif
   return(info);
}
