/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2006 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77GELS Mjoin(PRE,gels)
#elif defined (UpCase)
   #define F77GELS Mjoin(PREU,GELS)
#elif defined (Add_) || defined(Add__)
   #define F77GELS Mjoin(PRE,gels_)
#endif
#define f77gels Mjoin(PATL,f77gels)

int f77gels(const enum ATLAS_TRANS TA, const int M, const int N, const int NRHS,
            TYPE *A, const int lda, TYPE *B, const int ldb)
{
   #if defined(StringSunStyle)
      #if defined(ATL_FunkyInts)
         F77_INTEGER ONE=1;
      #else
         int ONE=1;
      #endif
   #elif defined(StringStructVal) || defined(StringStructPtr) || defined(StringCrayStyle)
      F77_CHAR ftrans;
   #endif
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77N=N, F77lda=lda, F77ldb=ldb, F77M=M, F77NRHS=NRHS;
      F77_INTEGER lwork, info;
   #else
      int info, lwork;
      #define F77M M
      #define F77N N
      #define F77NRHS NRHS
      #define F77lda lda
      #define F77ldb ldb
   #endif
   char ctrans;
   TYPE *work, wrk[2];

   if (TA == AtlasNoTrans) ctrans = 'N';
   else if (TA == AtlasTrans) ctrans = 'T';
   else ctrans = 'C';
   #if defined(StringSunStyle)
      #define args &ctrans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, \
                   work, &lwork, &info, ONE
   #elif defined(StringCrayStyle)
      ftrans = ATL_C2F_TransChar(cuplo);
      #define args ftrans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, \
                   work, &lwork, &info
   #elif defined(StringStructVal)
      ftrans.len = 1;
      ftrans.cp = &ctrans;
      #define args ftrans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, \
                   work, &lwork, &info
   #elif defined(StringStructPtr)
      ftrans.len = 1;
      ftrans.cp = &ctrans;
      #define args &ftrans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, \
                   work, &lwork, &info
   #else
      #define args NULL
      fprintf(stderr, "\n\nF77/C interface not defined!!\n\n");
      exit(-1);
   #endif
/*
 * Query routine for optimal workspace, allocate it, and call routine with it
 */
   work = wrk;
   lwork = -1;
   F77GELS(args);
   lwork = wrk[0];
   work = malloc(ATL_MulBySize(lwork));
   ATL_assert(work);
   info = 0;
   F77GELS(args);
   free(work);
   return(info);
}
