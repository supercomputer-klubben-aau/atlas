/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2007 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77.h"

#if defined(NoChange)
   #define F77POSV Mjoin(PRE,posv)
#elif defined (UpCase)
   #define F77POSV Mjoin(PREU,POSV)
#elif defined (Add_) || defined(Add__)
   #define F77POSV Mjoin(PRE,posv_)
#endif
#define f77posv Mjoin(PATL,f77posv)

int f77posv(const enum ATLAS_UPLO UPLO, const int N, const int NRHS,
            TYPE *A, const int lda, TYPE *B, const int ldb)
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
      const F77_INTEGER F77N=N, F77lda=lda, F77ldb=ldb, F77NRHS=NRHS;
      F77_INTEGER info;
   #else
      int info;
      #define F77N N
      #define F77NRHS NRHS
      #define F77lda lda
      #define F77ldb ldb
   #endif
   char cuplo;
   cuplo = (UPLO == AtlasUpper) ? 'U' : 'L';
   #if defined(StringSunStyle)
      #define args &cuplo, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, &info, ONE
   #elif defined(StringCrayStyle)
      ftrans = ATL_C2F_TransChar(cuplo);
      #define args cuplo, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, &info
   #elif defined(StringStructVal)
      fuplo.len = 1;
      fuplo.cp = &cuplo;
      #define args fuplo, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, &info
   #elif defined(StringStructPtr)
      fuplo.len = 1;
      fuplo.cp = &cuplo;
      #define args &fuplo, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb, &info
   #else
      #define args NULL
      fprintf(stderr, "\n\nF77/C interface not defined!!\n\n");
      exit(-1);
   #endif
   F77POSV(args);
   return(info);
}
