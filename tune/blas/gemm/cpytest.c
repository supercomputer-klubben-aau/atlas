/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "atlas_bitvec.h"
#ifndef TYPE
   #define ATL_CSZT const size_t
   #if defined(SREAL) || defined(SCPLX)
      #define TYPE float
   #elif defined(DREAL) || defined(DCPLX)
      #define TYPE double
   #endif
   #ifdef SREAL
      #define ATL_MulBySize(i_) ((i_)<<2)
      #define ATL_DivBySize(i_) ((i_)>>2)
   #elif defined(DREAL) || defined(SCPLX)
      #define ATL_MulBySize(i_) ((i_)<<3)
      #define ATL_DivBySize(i_) ((i_)>>3)
   #else
      #define ATL_MulBySize(i_) ((i_)<<4)
      #define ATL_DivBySize(i_) ((i_)>>4)
   #endif
   #if defined(SCPLX) || defined(DCPLX)
      #define SHIFT << 1
      #define TCPLX 1
   #else
      #define TREAL 1
      #define SHIFT
   #endif
   #ifdef TCPLX
      #define SCALAR TYPE*
      #define SCALAR_IS_ZERO(s_) ((s_)[0] == 0.0 && (s_)[1] == 0.0)
   #else
      #define SCALAR TYPE
      #define SCALAR_IS_ZERO(s_) ((s_) == 0.0)
   #endif
   #define ATL_Cachelen 128
   #define ATL_AlignPtr(vp) \
      (void*) ( ATL_Cachelen + (((((size_t) (vp))+127)>>7)<<7) )
#endif
#if defined(SREAL) || defined(SCPLX)
   #define FEPS 5.0e-7
#else
   #define FEPS 1.0e-15
#endif
#ifndef ATL_MU
   #define ATL_MU 1
#endif
#ifndef ATL_NU
   #define ATL_NU 1
#endif
#if TO_BLK
   #ifdef FromBlk_
      #undef FromBlk_
   #endif
#else
   #define FromBlk_
#endif
#ifdef Trans_
   #if Trans_ == 0
      #undef Trans_
   #endif
#endif
#ifdef Conj_
   #if Conj_ == 0
      #undef Conj_
   #endif
#endif
static int MA=0, NA=0, LDA=0, IA=0, JA=0;
#ifdef TCPLX
   #ifdef ALPHA1
      TYPE ALPHA[2] = {1.0, 0.0};
   #elif defined(ALPHAN) || defined(ALPHAN1)
      TYPE ALPHA[2] = {-1.0, 0.0};
   #elif defined(ALPHA0)
      TYPE ALPHA[2] = {0.0, 0.0};
   #else
      TYPE ALPHA[2] = { 1.55, -0.35};
   #endif
   #ifdef BETA1
      TYPE BETA[2] = {1.0, 0.0};
   #elif defined(BETAN) || defined(BETAN1)
      TYPE BETA[2] = {-1.0, 0.0};
   #elif defined(BETA0)
      TYPE BETA[2] = {0.0, 0.0};
   #else
      TYPE BETA[2] = {-0.25, 0.75};
   #endif
#else
   #ifdef ALPHA1
      TYPE ALPHA = 1.0;
   #elif defined(ALPHAN) || defined(ALPHAN1)
      TYPE ALPHA = -1.0;
   #elif defined(ALPHA0)
      TYPE ALPHA = 0.0;
   #else
      TYPE ALPHA = 1.08;
   #endif
   #ifdef BETA1
      TYPE BETA = 1.0;
   #elif defined(BETAN) || defined(BETAN1)
      TYPE BETA = -1.0;
   #elif defined(BETA0)
      TYPE BETA = 0.0;
   #else
      TYPE BETA = 0.25;
   #endif
#endif

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr,"  -A <M> <N> <LDA>: set dimensions of col-major matrix\n");
   fprintf(stderr,"  -b <M> <N> : set dims of access-major block\n");
   fprintf(stderr,"  -B # <M1> <N1> <M2> <N2> ... <M#> <N#>\n");
   fprintf(stderr,"  -s # <sz> allocate sz for each block\n");
   fprintf(stderr,"  -S # <sz1> ... <szN>: allocate sz for each block\n");
   fprintf(stderr,"  -i I J : offset in array to access\n");
   exit(ierr ? ierr : 1);
}

unsigned int *GetFlags (int nargs, char **args)
/*
 * RETURNS: int array, 1st elt is # of blocks, then M,N pairs of block dims
 */
{
   unsigned int nblk=0, i, SZ=0;
   unsigned int *bD=NULL, *szs=NULL;

   IA = JA = -1;
   for (i=1; i < nargs; i++)
   {
      char ch;
      int k;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 's':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         SZ = atoi(args[i]);
         break;
      case 'S':
         if (szs)
            free(szs);
         szs = NULL;
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         nblk = atoi(args[i]);
         if (bD)
            assert(bD[0] == nblk);
         else
         {
            bD = calloc(((nblk<<1)+nblk)+1,sizeof(int));
            assert(bD);
            bD[0] = nblk;
         }
         for (k=0; k < nblk; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            bD[(k<<1)+k+1] = atoi(args[i]);
         }
         break;
      case 'B': /* -B # m1 n1 ... m# n# */
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         nblk = atoi(args[i]);
         if (bD)
         {
            if (bD[3])
               assert(bD[0] == nblk);
            else
            {
               free(bD);
               bD = NULL;
            }
         }
         goto GET_BLKS;
      case 'b': /* -b m n */
         nblk = 1;
      GET_BLKS:
         if (!bD)
         {
            bD = calloc(((nblk<<1)+nblk)+1,sizeof(int));
            assert(bD);
         }
         *bD = nblk;
         nblk += nblk<<1;
         for (k=0; k < nblk; k += 3)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            bD[k+1] = atoi(args[i]);
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
            bD[k+2] = atoi(args[i]);
         }
         break;
      case 'i': /* -i <I> <J> */
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         IA = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         JA = atoi(args[i]);
         break;
      case 'A': /* -A M N lda */
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         MA = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         NA = atoi(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         LDA = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (IA == -1)
   {
      #ifdef FromBlk_  /* want to check pre/post pad! */
         IA = 8;
         JA = 6;
      #else
         IA = JA = 0;
      #endif
   }
   if (bD == NULL)
   {
      bD = malloc(4*sizeof(int));
      assert(bD);
      bD[0] = 1;
      bD[1] = 3*ATL_MU;
      bD[2] = 5*ATL_NU;
      bD[3] = bD[1] * bD[2];
   }
   else if (!bD[3])
   {
      int nblk = bD[0];
      nblk += nblk<<1;
      if (SZ)
         for (i=0; i < nblk; i += 3)
            bD[3+i] = SZ;
      else
         for (i=0; i < nblk; i += 3)
         {
            int M=bD[1+i], N=bD[2+i];
            M = ((M+ATL_MU-1)/ATL_MU)*ATL_MU;
            N = ((N+ATL_NU-1)/ATL_NU)*ATL_NU;
            bD[3+i] = N*M;
         }
   }
   if (!NA)
   {
      int k, nblk=bD[0];
      nblk += nblk<<1;
      MA = bD[1];
      NA = bD[2];
      for (k=3; k < nblk; k += 3)
      {
         MA = Mmax(MA, bD[k+1]);
         NA = Mmax(NA, bD[k+2]);
      }
      #ifdef Trans_
         k = MA;
         MA = NA;
         NA = k;
      #endif
      #ifdef FromBlk_
         NA += NA;
      #endif
      MA += IA;
      NA += JA;
      LDA = MA + 7;
   }
   else
   {
      int k, nblk=bD[0];
      nblk += nblk<<1;
      for (k=0; k < nblk; k += 3)
      {
         #ifdef Trans_
            assert(IA + bD[k+2] <= MA);
            assert(JA + bD[k+1] <= NA);
         #else
            assert(IA + bD[k+1] <= MA);
            assert(JA + bD[k+2] <= NA);
         #endif
      }
      assert(LDA >= MA);
   }
   return(bD);
}

TYPE *getMat(size_t M, size_t N, size_t lda, int seed)
{
   int i, j;
   const int lda2 = lda SHIFT;
   const double mul = 1.0 / (lda2*N);
   TYPE *A, *A0;

   A0 = A = malloc(N*lda2*sizeof(TYPE));
   assert(A);
   for (j=0; j < N; j++, A += lda2)
      for (i=0; i != lda2; i++)
         A[i] = 1.0 + (mul*j)*lda2+(i*1e-1);
   return(A0);
}

#define AlignPtr(vp_) (void*) (((((size_t)(vp_))+127)>>7)<<7)
TYPE *getBlock(int sz, int pad)
{
   void *vp;
   TYPE *b;
   size_t ip;
   int i, SZ = (sz + pad + pad)SHIFT;
   double mul = 1.0 / SZ;

   vp = malloc(SZ*sizeof(TYPE) + 128);
   b = AlignPtr(vp);
   for (i=0; i < SZ; i++)
      b[i] = -(1.0+mul*i);
   return(vp);
}

#ifdef matC_
   #ifdef Trans_
      #error "No test of Transpose C matrix right now!"
   #endif
   #ifdef FromBlk_
      #ifdef TCPLX
         void ATL_USERCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *r,
             const TYPE *c, const SCALAR beta, TYPE *C, const size_t ldc);
         void ATL_GOODCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *r,
             const TYPE *c, const SCALAR beta, TYPE *C, const size_t ldc);
      #else
         void ATL_USERCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *b,
             const SCALAR beta, TYPE *C, const size_t ldc);
         void ATL_GOODCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *b,
             const SCALAR beta, TYPE *C, const size_t ldc);
      #endif
   #else /* ToBlock */
      #ifdef TCPLX
         void ATL_USERCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *C,
             const size_t ldc, const SCALAR beta, TYPE *rb, TYPE *ib);
         void ATL_GOODCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *C,
             const size_t ldc, const SCALAR beta, TYPE *rb, TYPE *ib);
      #else
         void ATL_USERCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *C,
             const size_t ldc, const SCALAR beta, TYPE *b);
         void ATL_GOODCPMM
            (const size_t M, const size_t N, const SCALAR alpha, const TYPE *C,
             const size_t ldc, const SCALAR beta, TYPE *b);
      #endif
   #endif
#else  /* Trans & NoTrans have same prototype */
   #ifdef FromBlk_
      #ifdef TCPLX
         void ATL_USERCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           TYPE *rb, TYPE *ib, const TYPE *A, ATL_CSZT lda);
         void ATL_GOODCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           TYPE *rb, TYPE *ib, const TYPE *A, ATL_CSZT lda);
      #else
         void ATL_USERCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           TYPE *b, const TYPE *A, ATL_CSZT lda);
         void ATL_GOODCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           TYPE *b, const TYPE *A, ATL_CSZT lda);
      #endif
   #else
      #ifdef TCPLX
         void ATL_USERCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           const TYPE *A, ATL_CSZT lda, TYPE *rb, TYPE *ib);
         void ATL_GOODCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           const TYPE *A, ATL_CSZT lda, TYPE *rb, TYPE *ib);
      #else
         void ATL_USERCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           const TYPE *A, ATL_CSZT lda, TYPE *b);
         void ATL_GOODCPMM(ATL_CSZT K, ATL_CSZT N, const SCALAR alpha,
                           const TYPE *A, ATL_CSZT lda, TYPE *b);
      #endif
   #endif
#endif
void initMat(SCALAR beta, int M, int N, TYPE *A, size_t lda)
{
   size_t lda2 = lda SHIFT;
   unsigned int M2 = M SHIFT, i, j;
   if (SCALAR_IS_ZERO(beta))
   {
      TYPE nan = 0.0 / 0.0;
      for (j=0; j < N; j++, A += lda2)
         for (i=0; i < M2; i++)
            A[i] = nan;
   }
   else
   {
      const double mul = 1.0 / (((double)M2)*N);
      for (j=0; j < N; j++, A += lda2)
         for (i=0; i < M2; i++)
            A[i] = (i + j*lda2)*mul + 1.0;
   }
}

#ifdef matC_
   #define C_MAT 'C'
#else
   #define C_MAT 'A'
#endif
#ifdef FromBlk_
   #define C_DIR 'F'
#else
   #define C_DIR 'T'
#endif
#ifdef Trans_
   #ifdef Conj_
      #define C_TRANS 'C'
   #else
      #define C_TRANS 'T'
   #endif
#else
   #ifdef Conj_
      #define C_TRANS 'H'
   #else
      #define C_TRANS 'N'
   #endif
#endif

#ifdef FromBlk_
int doTests(unsigned int ntest, unsigned int *bD)
{
   int t;
   int nerr=0;
   const unsigned int I2 = IA SHIFT, M2 = MA SHIFT, LDA2 = LDA SHIFT;
   const double mul = 1.0 / (LDA2*NA);

   printf("RUN %u TESTS: MAT=%c, TA=%c MD=(%u,%u), LDA=%u, (%u,%u)\n",
          ntest, C_MAT, C_TRANS, MA, NA, LDA, I2, JA);
   for (t=0; t < ntest; t++)
   {
      unsigned int T = (t<<1)+t, mb = bD[T], nb = bD[T+1], sz=bD[T+2];
      #ifdef Trans_
         const unsigned int nrow = nb, ncol=mb;
      #else
         const unsigned int nrow = mb, ncol=nb;
      #endif
      const unsigned int M0=(IA+nrow)SHIFT, N0=JA+ncol, N1=N0+ncol;
      const size_t offT = LDA2*ncol;
      unsigned int j;
      void *v0;
      TYPE *b0, *b1;
      TYPE *A, *a0, *a1, *A0;
      #ifdef TCPLX
         TYPE *rb, *ib;
      #endif

      assert(N0 <= NA);
      assert(N1 <= NA);
      assert(M0 <= (MA SHIFT));
      printf("   Test %u, D=(%u,%u), SZ=%u\n", t, mb, nb, sz);
      fflush(stdout);
      v0 = getBlock(sz, 0);
      b0 = AlignPtr(v0);
      A0 = A = getMat(MA, NA, LDA, MA+NA*LDA);
      a0 = A+((IA+JA*LDA)SHIFT);
      a1 = a0 + offT;

      initMat(BETA, nrow, ncol, a0, LDA);
      initMat(BETA, nrow, ncol, a1, LDA);
      #ifdef TCPLX
         #ifdef matC_
            ATL_GOODCPMM(mb, nb, ALPHA, b0+sz, b0, BETA, a0, LDA);
            ATL_USERCPMM(mb, nb, ALPHA, b0+sz, b0, BETA, a1, LDA);
         #else
            ATL_GOODCPMM(nb, mb, ALPHA, b0+sz, b0, a0, LDA);
            ATL_USERCPMM(nb, mb, ALPHA, b0+sz, b0, a1, LDA);
         #endif
      #else
         #ifdef matC_
            ATL_GOODCPMM(mb, nb, ALPHA, b0, BETA, a0, LDA);
            ATL_USERCPMM(mb, nb, ALPHA, b0, BETA, a1, LDA);
         #else
            ATL_GOODCPMM(nb, mb, ALPHA, b0, a0, LDA);
            ATL_USERCPMM(nb, mb, ALPHA, b0, a1, LDA);
         #endif
      #endif
      a0 = A;
      for (j=0; j < NA; j++, A += LDA2)
      {
         unsigned int i;
         for (i=0; i < M2; i++)
         {
            double got, exp;
            if (i >= I2 && i < M0)
            {
               if (j >= JA && j < N0)
               {
                  int nanE, nanG;
                  exp = A[i];
                  got = A[i+offT];
                  nanE = exp != exp;
                  nanG = got != got;
                  if (nanE || nanG)
                  {
                     if (!nanE || !nanG)
                     {
                        nerr++;
                        fprintf(stderr,
                                "      A(%u,%u),b(%u,%u)=%e, expected=%e\n",
                                i, j, i-I2, j-JA, got, exp);
                     }
                  }
                  else
                  {
                     double diff;
                     diff = exp - got;
                     if (diff < 0.0)
                        diff = -diff;
                     if (diff > FEPS)
                     {
                        fprintf(stderr,
                                "      A(%u,%u),b(%u,%u)=%e, expected=%e\n",
                                i, j, i-I2, j-JA, got, exp);
                        nerr++;
                     }
                     continue;
                  }
               }
               else if (j <= N1)  /* don't scope 2nd overwrittes submat */
                  continue;       /* since was scoped right above */
            }
            exp = 1.0 + (mul*j)*LDA2+(i*1e-1);
            got = A[i];
            if (exp != got)
            {
               fprintf(stderr, "      PAD ERR, A(%u,%u)=%e, expected=%e\n",
                       i, j, got, exp);
               nerr++;
            }
         }
      }
      free(A0);
      free(v0);
      if (nerr > 1)
         return(nerr);
      printf("   DONE Test %u\n", t);
   }
   return(nerr);
}
#else
int doTests(unsigned int ntest, unsigned int *bD)
{
   int t;
   TYPE *A, *a;
   int nerr=0, pad=32, pad2=pad SHIFT;

   A = getMat(MA, NA, LDA, MA+NA*LDA);
   a = A+((IA+JA*LDA)SHIFT);
   printf("RUN %u TESTS: MAT=%c, TA=%c MD=(%u,%u), LDA=%u\n", ntest, C_MAT,
          C_TRANS, MA, NA, LDA);
   for (t=0; t < ntest; t++)
   {
      void *v0, *v1;
      TYPE *b0, *b1, *b2, *b3, *B0, *B1;
      unsigned int T = (t<<1)+t, mb = bD[T], nb = bD[T+1], sz=bD[T+2];
      unsigned int i, N;
      #ifdef TCPLX
         TYPE *rB0, *iB0, *rB1, *iB1;
      #endif
      printf("   Test %u, D=(%u,%u), SZ=%u\n", t, mb, nb, sz);
      v0 = getBlock(sz, pad);
      v1 = getBlock(sz, pad);
      b0 = AlignPtr(v0);
      b1 = AlignPtr(v1);
      B0 = b0 + pad2;
      B1 = b1 + pad2;
      b2 = B0 + (sz SHIFT);
      b3 = B1 + (sz SHIFT);
      #ifdef TCPLX
         #ifdef matC_
            ATL_GOODCPMM(mb, nb, ALPHA, a, LDA, BETA, B0+sz, B0);
            ATL_USERCPMM(mb, nb, ALPHA, a, LDA, BETA, B1+sz, B1);
         #else
            ATL_GOODCPMM(nb, mb, ALPHA, a, LDA, B0+sz, B0);
            ATL_USERCPMM(nb, mb, ALPHA, a, LDA, B1+sz, B1);
         #endif
      #else
         #ifdef matC_
            ATL_GOODCPMM(mb, nb, ALPHA, a, LDA, BETA, B0);
            ATL_USERCPMM(mb, nb, ALPHA, a, LDA, BETA, B1);
         #else
            ATL_GOODCPMM(nb, mb, ALPHA, a, LDA, B0);
            ATL_USERCPMM(nb, mb, ALPHA, a, LDA, B1);
         #endif
      #endif
      for (i=0; i < pad2; i++)
      {
         if (b0[i] != b1[i])
         {
            printf("      PREPAD  %d, expected=%e, got=%e\n", i, b0[i], b1[i]);
            nerr++;
         }
         if (b3[i] != b2[i])
         {
            printf("      POSTPAD %d, expected=%e, got=%e\n", i, b0[i], b1[i]);
            nerr++;
         }
      }
      for (N=sz SHIFT, i=0; i < N; i++)
      {
         TYPE diff = B0[i] - B1[i];

         diff = (diff >= 0.0) ? diff : -diff;
         if (B1[i] != B1[i])
         {
            fprintf(stderr, "      B[%u]: got NaN, expected=%e\n", i, B0[i]);
            if (B0[i] == B0[i])
               nerr++;
         }
         else if (diff > FEPS)
         {
            fprintf(stderr, "      B[%u]: got %e, expected=%e\n", i,
                    B1[i], B0[i]);
            nerr++;
         }
      }
      free(v0);
      free(v1);
      if (nerr > 1)
         return(nerr);
      printf("   DONE Test %u\n", t);
   }
   free(A);
   return(nerr);
}
#endif

int main(int nargs, char **args)
{
   int *bD;
   int nerr = 0;
   bD = GetFlags(nargs, args);
   nerr = doTests(bD[0], bD+1);
   if (!nerr)
      printf("\nSUCCESS: PASS %u ALL TESTS\n", bD[0]);
   else
      printf("\nFAILED WITH %u ERRORS\n", nerr);
   free(bD);
   return(nerr);
}
