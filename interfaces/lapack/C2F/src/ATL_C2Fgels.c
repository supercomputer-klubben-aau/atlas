/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_f77.h"
#include "atlas_C2Flapack.h"

#if defined(NoChange)
   #define F77GELS Mjoin(PRE,gels)
#elif defined (UpCase)
   #define F77GELS Mjoin(PREU,GELS)
#elif defined (Add_) || defined(Add__)
   #define F77GELS Mjoin(PRE,gels_)
#endif
#define PC2F Mjoin(ATL_C2F,PRE)

int Mjoin(PC2F,gels_wrk)(const enum CBLAS_TRANSPOSE TA, ATL_CINT M, ATL_CINT N,
                         ATL_CINT NRHS, TYPE *A, ATL_CINT lda,
                         TYPE *B, ATL_CINT ldb, TYPE *wrk, ATL_INT lwrk)
{
#if defined(StringSunStyle)
   F77_INTEGER ONE=1;
#elif defined(StringStructVal) || defined(StringStructPtr) || \
      defined(StringCrayStyle)
   F77_CHAR F77trans;
#endif
   F77_INTEGER F77M=M, F77N=N, F77NRHS=NRHS, F77lda=lda, F77ldb=ldb,
               F77lwrk = lwrk, F77info;
#if defined(StringSunStyle)
   void F77GELS(char*, F77_INTEGER*, F77_INTEGER*, F77_INTEGER*, TYPE*,
               F77_INTEGER*, TYPE*, F77_INTEGER*, TYPE*, F77_INTEGER*,
               F77_INTEGER*, F77_INTEGER);
#elif defined(StringStructPtr)
   void F77GELS(F77_CHAR*, F77_INTEGER*, F77_INTEGER*, F77_INTEGER*, TYPE*,
               F77_INTEGER*, TYPE*, F77_INTEGER*, TYPE*, F77_INTEGER*,
               F77_INTEGER*);
#else
   void F77GELS(F77_CHAR, F77_INTEGER*, F77_INTEGER*, F77_INTEGER*, TYPE*,
               F77_INTEGER*, TYPE*, F77_INTEGER*, TYPE*, F77_INTEGER*,
               F77_INTEGER*);
#endif
   char ctrans;

   if (TA == CblasNoTrans) ctrans = 'N';
   else if (TA == CblasTrans) ctrans = 'T';
   else ctrans = 'C';
#if defined(StringCrayStyle)
   f77trans = ATL_C2F_TransChar(cuplo);
#elif defined(StringStructVal) || defined(StringStructPtr)
   f77trans.len = 1;
   f77trans.cp = &ctrans;
#elif !defined(StringSunStyle)
   fprintf(stderr, "\n\nF77/C interface not defined!!\n\n");
   exit(-1);
#endif

#if defined(StringSunStyle)
   F77GELS(&ctrans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb,
           wrk, &F77lwrk, &F77info, ONE);
#elif defined(StringStructPtr)
   F77GELS(&f77trans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb,
           wrk, &F77lwrk, &F77info);
#elif defined(StringStructVal) || defined(StringCrayStyle)
   F77GELS(f77trans, &F77M, &F77N, &F77NRHS, A, &F77lda, B, &F77ldb,
           wrk, &F77lwrk, &F77info);
#endif
   return(F77info);
}

int Mjoin(PC2F,gels)(const enum CBLAS_TRANSPOSE TA, ATL_CINT M, ATL_CINT N,
                     ATL_CINT NRHS, TYPE *A, ATL_CINT lda,
                     TYPE *B, ATL_CINT ldb)
{
   TYPE work[2];
   TYPE *wrk;
   ATL_INT lwrk;
   int iret;
/*
 * Query routine for optimal workspace, allocate it, and call routine with it
 */
   ATL_assert(!Mjoin(PC2F,gels_wrk)(TA, M, N, NRHS, A, lda, B, ldb, work, -1));
   lwrk = work[0];
   wrk = malloc(ATL_MulBySize(lwrk));
   ATL_assert(wrk);
   iret = Mjoin(PC2F,gels_wrk)(TA, M, N, NRHS, A, lda, B, ldb, wrk, lwrk);
   free(wrk);
   return(iret);
}
