/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GERQF Mjoin(PRE,gerqf)
#elif defined (UpCase)
   #define F77GERQF Mjoin(PREU,GERQF)
#elif defined (Add_) || defined(Add__)
   #define F77GERQF Mjoin(PRE,gerqf_)
#endif
#define f77gerqf Mjoin(PATL,f77gerqf)

int f77gerqf(const enum ATLAS_ORDER Order, const int M, const int N,
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

   F77GERQF(&F77M, &F77N, A, &F77lda, tau, work, &F77lwork, &F77info);
   #ifndef ATL_FunkyInts
   #endif

   return(F77info);
}
