/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GETRF Mjoin(PRE,getrf)
#elif defined (UpCase)
   #define F77GETRF Mjoin(PREU,GETRF)
#elif defined (Add_) || defined(Add__)
   #define F77GETRF Mjoin(PRE,getrf_)
#endif
#define f77getrf Mjoin(PATL,f77getrf)

int f77getrf(const enum ATLAS_ORDER Order, const int M, const int N,
             TYPE *A, const int lda, int *ipiv)
{
   int i;
   const int MN=Mmin(M,N);
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77M=M, F77N=N, F77lda=lda;
      F77_INTEGER info, *F77ipiv;
   #else
      int info;
      #define F77M M
      #define F77N N
      #define F77lda lda
      #define F77ipiv ipiv
   #endif
   ATL_assert(Order == AtlasColMajor);
   #ifdef ATL_FunkyInts
      F77ipiv = malloc(MN * sizeof(F77_INTEGER));
      ATL_assert(F77ipiv);
   #else
      #define F77ipiv ipiv
   #endif

   F77GETRF(&F77M, &F77N, A, &F77lda, F77ipiv, &info);

   #ifdef ATL_FunkyInts
      for (i=0; i < MN; i++) ipiv[i] = F77ipiv[i] - 1;
      free(F77ipiv);
   #else
      for (i=0; i < MN; i++) ipiv[i]--;
   #endif
   return(info);
}
