/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2001 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GETRI Mjoin(PRE,getri)
#elif defined (UpCase)
   #define F77GETRI Mjoin(PREU,GETRI)
#elif defined (Add_) || defined(Add__)
   #define F77GETRI Mjoin(PRE,getri_)
#endif
#define f77getri Mjoin(PATL,f77getri)

int f77getri(const enum ATLAS_ORDER Order, const int N, TYPE *A, const int lda,
             int *ipiv, TYPE *work, int lwork)
{
   int i;
   const int MN=N;
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77N=N, F77lda=lda, F77lwork=lwork;
      F77_INTEGER info, *F77ipiv;
   #else
      int info;
      #define F77M M
      #define F77N N
      #define F77lda lda
      #define F77ipiv ipiv
      #define F77lwork lwork
   #endif
   ATL_assert(Order == AtlasColMajor);
   #ifdef ATL_FunkyInts
      F77ipiv = malloc(MN * sizeof(F77_INTEGER));
      ATL_assert(F77ipiv);
   #else
      #define F77ipiv ipiv
   #endif

   #ifdef ATL_FunkyInts
      F77lwork = lwork;
      for (i=0; i < MN; i++) F77ipiv[i] = ipiv[i] + 1;
   #else
      for (i=0; i < MN; i++) ipiv[i]++;
   #endif
   F77GETRI(&F77N, A, &F77lda, F77ipiv, work, &F77lwork, &info);
   #ifdef ATL_FunkyInts
      for (i=0; i < MN; i++) ipiv[i] = F77ipiv[i] - 1;
      free(F77ipiv);
   #else
      for (i=0; i < MN; i++) ipiv[i]--;
   #endif
   return(info);
}
