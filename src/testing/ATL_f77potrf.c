/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77POTRF Mjoin(PRE,potrf)
#elif defined (UpCase)
   #define F77POTRF Mjoin(PREU,POTRF)
#elif defined (Add_) || defined(Add__)
   #define F77POTRF Mjoin(PRE,potrf_)
#endif
#define f77potrf Mjoin(PATL,f77potrf)

int f77potrf(const enum ATLAS_UPLO Uplo, const int N, TYPE *A, const int lda)
{
   #if defined(StringSunStyle)
      #if defined(ATL_FunkyInts)
         F77_INTEGER ONE=1;
      #else
         int ONE=1;
      #endif
   #elif defined(StringStructVal) || defined(StringStructPtr) || defined(StringCrayStyle)
      F77_CHAR fuplo;
   #endif
   #ifdef ATL_FunkyInts
      const F77_INTEGER F77N=N, F77lda=lda;
      F77_INTEGER info;
   #else
      int info;
      #define F77N N
      #define F77lda lda
   #endif
   char cuplo;

   if (Uplo == AtlasUpper) cuplo = 'U';
   else cuplo = 'L';
   #if defined(StringSunStyle)
      F77POTRF(&cuplo, &F77N, A, &F77lda, &info, ONE);
   #elif defined(StringCrayStyle)
      fuplo = ATL_C2F_TransChar(cuplo);
      F77POTRF(fuplo, &F77N, A, &F77lda, &info);
   #elif defined(StringStructVal)
      fuplo.len = 1;
      fuplo.cp = &cuplo;
      F77POTRF(fuplo, &F77N, A, &F77lda, &info);
   #elif defined(StringStructPtr)
      fuplo.len = 1;
      fuplo.cp = &cuplo;
      F77POTRF(&fuplo, &F77N, A, &F77lda, &info);
   #else
      fprintf(stderr, "\n\nF77/C interface not defined!!\n\n");
      exit(-1);
   #endif
   return(info);
}
