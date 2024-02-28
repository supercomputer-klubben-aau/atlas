/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GEQLF Mjoin(PRE,geqlf)
#elif defined (UpCase)
   #define F77GEQLF Mjoin(PREU,GEQLF)
#elif defined (Add_) || defined(Add__)
   #define F77GEQLF Mjoin(PRE,geqlf_)
#endif
#define f77geqlf Mjoin(PATL,f77geqlf)

int f77geqlf(const enum ATLAS_ORDER Order, const int M, const int N,
             TYPE *A, const int lda, TYPE *tau, TYPE *work, int lwork)
{
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77M=M, F77N=N, F77lda=lda, F77lwork=lwork, F77info;
   #else
      int F77info;
      #define F77M M
      #define F77N N
      #define F77lda lda
      #define F77lwork lwork
   #endif
   ATL_assert(Order == AtlasColMajor);

   F77GEQLF(&F77M, &F77N, A, &F77lda, tau, work, &F77lwork, &F77info);
   #ifndef ATL_FunkyInts
   #endif

   return(F77info);
}
