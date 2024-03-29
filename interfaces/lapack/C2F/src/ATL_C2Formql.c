/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2008 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_f77.h"
#include "atlas_C2Flapack.h"

#if defined(NoChange)
   #define F77ORMQL Mjoin(PRE,ormql)
#elif defined (UpCase)
   #define F77ORMQL Mjoin(PREU,ORMQL)
#elif defined (Add_) || defined(Add__)
   #define F77ORMQL Mjoin(PRE,ormql_)
#endif
#define PC2F Mjoin(ATL_C2F,PRE)

int Mjoin(PC2F,ormql_wrk)
   (const enum CBLAS_SIDE Side, const enum CBLAS_TRANSPOSE TA,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, TYPE *A, ATL_CINT lda, TYPE *TAU,
    TYPE *C, ATL_CINT ldc, TYPE *wrk, ATL_INT lwrk)
{
#if defined(StringSunStyle)
   F77_INTEGER ONE=1;
#elif defined(StringStructVal) || defined(StringStructPtr) || \
      defined(StringCrayStyle)
   F77_CHAR F77trans, F77Side;
#endif
   F77_INTEGER F77M=M, F77N=N, F77K=K, F77lda=lda, F77ldc=ldc,
               F77lwrk = lwrk, F77info;
#if defined(StringSunStyle)
   void F77ORMQL(char*, char*, F77_INTEGER*, F77_INTEGER*, F77_INTEGER*,
                 TYPE*, F77_INTEGER*, TYPE*, TYPE*, F77_INTEGER*, TYPE*,
                 F77_INTEGER*, F77_INTEGER*, F77_INTEGER, F77_INTEGER);
#elif defined(StringStructPtr)
   void F77ORMQL(F77_CHAR*, F77CHAR*,  F77_INTEGER*, F77_INTEGER*,
                 F77_INTEGER*, TYPE*, F77_INTEGER*, TYPE*, TYPE*, F77_INTEGER*,
                 TYPE*, F77_INTEGER*, F77_INTEGER*);
#else
   void F77ORMQL(F77_CHAR, F77_CHAR, F77_INTEGER*, F77_INTEGER*,
                 F77_INTEGER*, TYPE*, F77_INTEGER*, TYPE*, TYPE*, F77_INTEGER*,
                 TYPE*, F77_INTEGER*, F77_INTEGER*);
#endif
   char cside, ctrans;

   if (TA == CblasNoTrans) ctrans = 'N';
   else if (TA == CblasTrans) ctrans = 'T';
   else ctrans = 'C';
   if (Side == CblasLeft) cside = 'L';
   else cside = 'R';
#if defined(StringCrayStyle)
   f77side  = ATL_C2F_TransChar(cside);
   f77trans = ATL_C2F_TransChar(cuplo);
#elif defined(StringStructVal) || defined(StringStructPtr)
   f77trans.len = 1;
   f77trans.cp = &ctrans;
   f77side.len = 1;
   f77side.cp = &cside;
#elif !defined(StringSunStyle)
   fprintf(stderr, "\n\nF77/C interface not defined!!\n\n");
   exit(-1);
#endif

#if defined(StringSunStyle)
   F77ORMQL(&cside, &ctrans, &F77M, &F77N, &F77K, A, &F77lda, TAU,
            C, &F77ldc, wrk, &F77lwrk, &F77info, ONE, ONE);
#elif defined(StringStructPtr)
   F77ORMQL(&f77side, &f77trans, &F77M, &F77N, &F77K, A, &F77lda, TAU,
            C, &F77ldc, wrk, &F77lwrk, &F77info);
#elif defined(StringStructVal) || defined(StringCrayStyle)
   F77ORMQL(&f77side, &f77trans, &F77M, &F77N, &F77K, A, &F77lda, TAU,
            C, &F77ldc, wrk, &F77lwrk, &F77info);
#endif
   return(F77info);
}

int Mjoin(PC2F,ormql)
   (const enum CBLAS_SIDE Side, const enum CBLAS_TRANSPOSE TA,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, TYPE *A, ATL_CINT lda, TYPE *TAU,
    TYPE *C, ATL_CINT ldc)
{
   TYPE work[2];
   void *vp;
   TYPE *wrk;
   ATL_INT lwrk;
   int iret;
/*
 * Query routine for optimal workspace, allocate it, and call routine with it
 */
   ATL_assert(!Mjoin(PC2F,ormql_wrk)(Side, TA, M, N, K, A, lda, TAU, C, ldc,
                                     work, -1));
   lwrk = work[0];
   vp = malloc(ATL_MulBySize(lwrk) + ATL_Cachelen);
   ATL_assert(vp);
   wrk = ATL_AlignPtr(vp);
   iret = Mjoin(PC2F,ormql_wrk)(Side, TA, M, N, K, A, lda, TAU, C, ldc,
                                wrk, lwrk);
   free(vp);
   return(iret);
}
