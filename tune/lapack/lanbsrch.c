/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2009 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_lapack.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias_lapack.h"
#endif
#ifdef ATL_ilaenv
   #undef ATL_ilaenv
#endif
#include "cblas.h"
#include "atlas_cblastypealias.h"
#include "atlas_tst.h"
#include "atlas_level1.h"
#include "atlas_lvl3.h"
#include <string.h>
#include <ctype.h>
#ifdef GCCWIN
   ___main(){} __main(){} MAIN__(){} _MAIN_(){}
#endif

#define CBP Mjoin(cblas_,PRE)


double time00();

#define CAN_NB 32
#include "C_lapack.h"

static int Uplo_LA2ATL(int i)
{
   return( (i & LAUpper) ? AtlasUpper : AtlasLower);
}
static int Side_LA2ATL(int i)
{
   return( (i & LARight) ? AtlasRight : AtlasLeft);
}

#ifdef TimeF77
   #define test_potrf(Uplo_, N_, A_, lda_) \
      ATL_assert(Mjoin(PATL,f77potrf)(Uplo_LA2ATL(Uplo_), N_, A_, lda_) == 0)
#elif defined(TimeC)
   #include "clapack.h"
   #define Cpotrf Mjoin(Mjoin(clapack_,PRE),potrf)
   #define test_potrf(Uplo_, N_, A_, lda_) \
      ATL_assert(Cpotrf(AtlasColMajor, Uplo_LA2ATL(Uplo_), N_, A_, lda_) == 0)
#else
   #define test_potrf(Uplo_, N_, A_, lda_) \
      ATL_assert(ATL_potrf(AtlasColMajor, Uplo_LA2ATL(Uplo_), N_, A_, lda_)==0)
#endif
#ifdef TimeF77
   #define test_getrf(Major_, M_, N_, A_, lda_, ipiv_) \
      ATL_assert(Mjoin(PATL,f77getrf)(Major_, M_, N_, A_, lda_, ipiv_) == 0)
#elif defined(TimeC)
   #include "clapack.h"
   #define Cgetrf Mjoin(Mjoin(clapack_,PRE),getrf)
   #define test_getrf(Major_, M_, N_, A_, lda_, ipiv_) \
      ATL_assert(Cgetrf(Major_, M_, N_, A_, lda_, ipiv_) == 0)
#elif defined(TimeBCAMM)
   #include "atlas_bcamm.h"
   #define test_getrf(Major_, M_, N_, A_, lda_, ipiv_) \
      ATL_assert(Mjoin(PATL,tgetrf_bcAmm)(Major_, M_, N_, A_, lda_, ipiv_) == 0)
#else
   #define test_getrf(Major_, M_, N_, A_, lda_, ipiv_) \
      ATL_assert(ATL_getrf(Major_, M_, N_, A_, lda_, ipiv_) == 0)
#endif
#ifndef TimeF77
   #include "atlas_lapack.h"
#endif
#ifdef TimeF77
   #define test_gelqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
   ATL_assert(!Mjoin(PATL,f77gelqf)(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_))
#elif defined(TimeC)
   #define test_gelqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_gelqf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#else
   #include "clapack.h"
   #define test_gelqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_gelqf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#endif
#ifdef TimeF77
   #define test_gerqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
   ATL_assert(!Mjoin(PATL,f77gerqf)(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_))
#elif defined(TimeC)
   #define test_gerqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_gerqf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#else
   #include "clapack.h"
   #define test_gerqf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_gerqf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#endif
#ifdef TimeF77
   #define test_geqlf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
   ATL_assert(!Mjoin(PATL,f77geqlf)(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_))
#elif defined(TimeC)
   #define test_geqlf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_geqlf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#else
   #include "clapack.h"
   #define test_geqlf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_geqlf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#endif
#ifdef TimeF77
   #define test_geqrf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
   ATL_assert(!Mjoin(PATL,f77geqrf)(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_))
#elif defined(TimeC)
   #define test_geqrf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_geqrf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#else
   #include "clapack.h"
   #define test_geqrf(Major_, M_, N_, A_, lda_, tau_, wrk_, lw_) \
      ATL_assert(!ATL_geqrf(M_, N_, A_, lda_, tau_, wrk_, lw_))
#endif


TYPE *GetGE(int M, int N, int lda)
{
   TYPE *A;
   A = malloc(ATL_MulBySize(lda)*N);
   if (A)
   {
      #if defined(ATL_USEPTHREADS) && !defined(ATL_NONUMATOUCH)
         ATL_NumaTouchSpread(ATL_MulBySize(lda)*N, A);
      #endif
      Mjoin(PATL,gegen)(M, N, A, lda, M*N+lda);
   }
   return(A);
}

static void CrapUpTri
   (enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, int N, TYPE *A, int lda)
/*
 * Puts crap on opposite triangle to Uplo, so as to ensure error on use
 */
{
   const int lda2=(lda SHIFT), ldap1=((lda+1)SHIFT);
   int j;

   if (Order == CblasRowMajor)
   {
      if (Uplo == CblasLower) Uplo = CblasUpper;
      else Uplo = CblasLower;
   }
   if (Uplo == CblasLower)
   {
      A += lda2;
      for (j=1; j < N; j++, A += lda2)
         Mjoin(PATLU,set)(j SHIFT, -50000000.0, A, 1);
   }
   else
   {
      for (j=0; j < N; j++, A += ldap1)
         Mjoin(PATLU,set)((N-j-1)SHIFT, -5500000.0, A+(1 SHIFT), 1);
   }
}

static TYPE *DupMat(enum ATLAS_ORDER Order, int M, int N, TYPE *A, int lda,
                    int ldc)
/*
 * returns a duplicate of the A matrix, with new leading dimension
 */
{
   int i, j, M2;
   const int ldc2 = (ldc SHIFT), lda2 = (lda SHIFT);
   TYPE *C;
   if (Order == CblasRowMajor)
   {
      i = M;
      M = N;
      N = i;
   }
   M2 = M SHIFT;
   ATL_assert(ldc >= M);
   C = malloc(ATL_MulBySize(ldc)*N);
   ATL_assert(C);
   #if defined(ATL_USEPTHREADS) && !defined(ATL_NONUMATOUCH)
      ATL_NumaTouchSpread(ATL_MulBySize(ldc)*N, C);
   #endif
   for (j=0; j != N; j++)
   {
      for (i=0; i != M2; i++) C[i] = A[i];
      C += ldc2;
      A += lda2;
   }
   return(C-N*ldc2);
}

#include <math.h>
static void PosDefGen
   (enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, int N, TYPE *A, int lda)
/*
 * Generates a reasonably conditioned positive definite matrix
 */
{
   TYPE *aa, *L;
   TYPE val, bias, sign;
   int j;
   const int lda2=(lda SHIFT), ldap1=((lda+1)SHIFT);

   Mjoin(PATL,gegen)(N, N, A, lda, N*N+lda);
   if (Order == CblasRowMajor)
   {
      if (Uplo == CblasLower) Uplo = CblasUpper;
      else Uplo = CblasLower;
   }
/*
 * It should be enough to make diagonal non-zero, but small numbers are very
 * ill-conditioned, and therefore may not be solvable in practice.  Therefore,
 * scale the diagonal by log(N).
 */
   bias = log(N);
   bias = (bias < 1.0) ? 1.0 : bias;
   for (aa=A,j=0; j < N; j++, aa += ldap1)
   {
      val = *aa;
      sign = (val < 0.0) ? -1.0 : 1.0;
      val = (val < 0.0) ? -val : val;
      val = (val+bias)*sign;
      *aa = val;
   }
/*
 * For imaginary numbers, force zero imaginary component on diagonal
 */
   #ifdef TCPLX
      Mjoin(Mjoin(ATL_,UPR),set)(N, 0.0, A+1, ldap1);
   #endif
/*
 * Zero non-active portion of matrix
 */
   if (Uplo == CblasLower)
   {
      for (j=0, aa=A; j < N; j++, aa += lda2)
         Mjoin(PATL,zero)(j, aa, 1);
   }
   else
   {
      for (j=0, aa=A+(1 SHIFT); j < N; j++, aa += ldap1)
         Mjoin(PATL,zero)(N-j-1, aa, 1);
   }
/*
 * Force A = L * L', where L is Lower or Upper as requested, to make pos def
 */
   L = DupMat(CblasColMajor, N, N, A, lda, N);
   #ifdef TCPLX
      Mjoin(CBP,herk)(CblasColMajor, Uplo, CblasNoTrans, N, N, ATL_rone, L, N,
                      ATL_rzero, A, lda);
   #else
      Mjoin(CBP,syrk)(CblasColMajor, Uplo, CblasNoTrans, N, N, ATL_rone, L, N,
                      ATL_rzero, A, lda);
#endif
   free(L);
/*
 * Make sure non-triangular elements are bad for error detection
 */
   CrapUpTri(CblasColMajor, Uplo, N, A, lda);
}

static void MakeHEDiagDom
   (enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, int N, TYPE *A, int lda)
/*
 * Makes hermitian matrix diagonally dominant
 */
{
   int j;
   const int lda2=(lda SHIFT), ldap1=((lda+1)SHIFT);

   if (Order == CblasRowMajor)
   {
      if (Uplo == CblasLower) Uplo = CblasUpper;
      else Uplo = CblasLower;
   }
   if (Uplo == CblasLower)
   {
      for (j=0; j < N; j++, A += ldap1)
      {
         #ifdef TREAL
            *A = 1.0 + cblas_asum(N-j, A, 1);
            *A += cblas_asum(j, A-lda*j, lda);
         #elif defined(SCPLX)
            *A = 1.0 + cblas_scasum(N-j, A, 1);
            *A += cblas_scasum(j, A-lda2*j, lda);
         #else
            *A = 1.0 + cblas_dzasum(N-j, A, 1);
            *A += cblas_dzasum(j, A-lda2*j, lda);
         #endif
         #ifdef TCPLX
            A[1] = ATL_rzero;
         #endif
      }
   }
   else /* Upper */
   {
      for (j=0; j < N; j++, A += ldap1)
      {
         #ifdef TREAL
            *A = 1.0 + cblas_asum(N-j, A, lda);
            *A += cblas_asum(j, A-j, 1);
         #else
            #ifdef SCPLX
               *A = 1.0 + cblas_scasum(N-j, A, lda);
               *A += cblas_scasum(j, A-j*2, 1);
            #else
               *A = 1.0 + cblas_dzasum(N-j, A, lda);
               *A += cblas_dzasum(j, A-j*2, 1);
            #endif
            A[1] = ATL_rzero;
         #endif
      }
   }
}

static void hegen
   (enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, int N, TYPE *A, int lda)
{
#ifdef POSDEFGEN
   PosDefGen(Order, Uplo, N, A, lda);
#else
   MakeHEDiagDom(Order, Uplo, N, A, lda);
   CrapUpTri(Order, Uplo, N, A, lda);
#endif
}

static TYPE *GetHE(enum CBLAS_ORDER Order, enum CBLAS_UPLO Uplo, int N, int lda)
/*
 * Gets symm/hemm matrix, and puts a bunch of crap in other side to make
 * sure factorization doesn't use it, and makes pos def by making it
 * diag dominant
 */
{
   TYPE *A;

   A = GetGE(N, N, lda);
   if (!A) return(NULL);
   hegen(Order, Uplo, N, A, lda);
   return(A);
}
static int ONB=0;       /* optimal NB to return in ILAENV */
#if defined(ATL_USEPTHREADS) && defined(ATL_USE_ATL_ILAENV) &&  \
    defined(ATL_USER_ILAENV)
   #define ATL_ilaenv ATL_itlaenv
#endif
/*
 * These are ISPEC values, which control what question is being asked ILAENV
 */
#define LA_OPT_NB 1     /* best NB if workspace allows */
#define LA_MIN_NB 2     /* minimal nb that will give speedup */
#define LA_NBXOVER 3    /* when N < this, use unblocked code */

#if defined(Add_) || defined(Add__)
   #define ILAENV ilaenv_
#elif defined(NoChange)
   #define ILAENV ilaenv
#elif !defined(UpCase)  /* no #define necessary for this UpCase! */
   #error "Unknown name decoration!"
#endif

#ifndef ATL_USER_ILAENV /* ignore this ILAENV if user provides */
#ifdef ATL_USE_ATL_ILAENV
   #include "atlas_lapack.h"
#endif
F77_INTEGER ILAENV
#ifdef StringStructVal
   (F77_INTEGER *ISPEC, F77_CHAR NAME, F77_CHAR OPTS,
#elif defined(StringStructPtr)
   (F77_INTEGER *ISPEC, F77_CHAR *NAME, F77_CHAR *OPTS,
#elif defined(StringSunStyle)
   (F77_INTEGER *ISPEC, char *NAME, char *OPTS,
#else
    #error "No supported string handling for this architecture!"
#endif
    F77_INTEGER *N1, F77_INTEGER *N2, F77_INTEGER *N3, F77_INTEGER *N4
#if defined(StringSunStyle)
    , F77_INTEGER namelen, F77_INTEGER optslen
#endif
   )
{
   char *name, *opts, *fname, *fopts;
   int ispec, n1, n2, n3, n4;
#ifndef StringSunStyle
   int namelen, optslen;
#endif
   int nb, i;
#ifdef ATL_USE_ATL_ILAENV
   int irout, iopt;
#endif

/*
 * Translate F77 inputs into C; strings upcased for ease of comparison
 */
#ifdef StringStructVal
   namelen = NAME.len;
   fname = NAME.cp;
   optslen = OPTS.len;
   fopts = OPTS.cp;
#elif defined(StringStructPtr)
   namelen = NAME->len;
   fname = NAME->cp;
   optslen = OPTS->len;
   fopts = OPTS->cp;
#elif defined(StringSunStyle)
   fname = NAME;
   fopts = OPTS;
#endif
   name = malloc(sizeof(char)*(namelen+1));
   opts = malloc(sizeof(char)*(optslen+1));
   ATL_assert(name && opts);
   for (i=0; i < namelen; i++)
      name[i] = toupper(fname[i]);
   name[i] = '\0';
   for (i=0; i < optslen; i++)
      opts[i] = toupper(fopts[i]);
   opts[i] = '\0';
   ispec = *ISPEC; n1 = *N1; n2 = *N2; n3 = *N3; n4 = *N4;
/*
 * Should ask for min blk factor only when workspace inadequate, which
 * should never occur since we malloc the space
 */
   ATL_assert(ispec != LA_MIN_NB);
   if (ispec == 3)                      /* ask for 8 cols before blocking */
      return((F77_INTEGER) 8);
   if (ispec != LA_OPT_NB)              /* only query left in tuning calls */
   {
      fprintf(stderr, "unexpected ISPEC : %d\n", ispec);
      exit(-1);
   }
   ATL_assert(ispec == LA_OPT_NB);

/*
 * If we use ATLAS's tuned ATL_ilaenv, just report selected NB in ONB
 */
#ifdef ATL_USE_ATL_ILAENV
   irout = -1;
   iopt = 0;
   if (!strcmp(name+1, "ORMQR"))
      irout = LArorgen;
   else if (!strcmp(name+1, "GETRF"))
      irout = LAgetrf;
   else if (!strcmp(name+1, "POTRF"))
      irout = LApotrf;
   else if (!strcmp(name+1, "GEQRF"))
   {
      irout = LAgeqrf;
      iopt += LARight+LAUpper;
   }
   else if (!strcmp(name+1, "GELQF"))
   {
      irout = LAgeqrf;
      iopt += LALeft+LALower;
   }
   else if (!strcmp(name+1, "GERQF"))
   {
      irout = LAgeqrf;
      iopt += LALeft+LAUpper;
   }
   else if (!strcmp(name+1, "GEQLF"))
   {
      irout = LAgeqrf;
      iopt += LARight+LALower;
   }
   ATL_assert(irout != -1)
   switch (*name)
   {
   case 'D':
      iopt += LADreal;
      break;
   case 'S':
      iopt += LASreal;
      break;
   case 'Z':
      iopt += LADcplx;
      break;
   case 'C':
      iopt += LAScplx;
      break;
   default :
      ATL_assert(0 == 1);
   }
   ONB = ATL_ilaenv(ispec, irout, iopt, n1, n2, n3, n4);
#else
/*
 * NOTE: need to scope what DORMQR should return!
 */
   if ( strcmp(name+1, "GEQRF") &&       /* QR, M >= N */
        strcmp(name+1, "GELQF") &&       /* LQ  */
        strcmp(name+1, "GERQF") &&       /* RQ  */
        strcmp(name+1, "GEQLF") &&       /* QL  */
        strcmp(name+1, "ORMQR") &&       /* QR/LQ wrk rout */
        strcmp(name+1, "GETRF") &&       /* LU */
        strcmp(name+1, "POTRF") )        /* Cholesky */
   {
      fprintf(stderr, "Unexpected name request = '%s'!\n", name);
      exit(-2);
   }
#endif
   return(ONB);
}
/*
 * same as before, but overriding the C interface to ilaenv
 */
int clapack_ilaenv(enum ATL_ISPEC ispec, enum ATL_LAROUT rout,
                   unsigned int opts, int N1, int N2, int N3, int N4)
{
   if (ispec == 3)
      return(8);                     /* ask for 8 cols before blocking */
   ATL_assert(ispec == LA_OPT_NB);   /* only query left in tuning calls */
   #ifdef ATL_USE_ATL_ILAENV
      ONB = ATL_ilaenv(ispec, rout, opts, N1, N2, N3, N4);
   #endif
   return(ONB);
}
/*
 * Override ATL_ilaenv & threaded version as well
 */
#ifndef ATL_USE_ATL_ILAENV
int ATL_ilaenv(enum ATL_ISPEC ispec, enum ATL_LAROUT rout,
               unsigned int opts, int N1, int N2, int N3, int N4)
{
   return(clapack_ilaenv(ispec, rout, opts, N1, N2, N3, N4));
}
int ATL_itlaenv(enum ATL_ISPEC ispec, enum ATL_LAROUT rout,
                unsigned int opts, int N1, int N2, int N3, int N4)
{
   return(clapack_ilaenv(ispec, rout, opts, N1, N2, N3, N4));
}
#endif
#endif   /* end #ifndef ATL_USER_ILAENV */

static FILE *my_fopen(char *name, char *acc)
{
   if (!name)
      return(*acc == 'r' ? stdin : stdout);
   if (!strncmp(name, "stderr", 6))
      return(stderr);
   if (!strncmp(name, "stdout", 6))
      return(stdout);
   if (!strncmp(name, "stdin", 5))
      return(stdin);
   return(fopen(name, acc));
}

static void my_fclose(FILE *fp)
{
   if (fp && fp != stderr && fp != stdout && fp != stdin)
      fclose(fp);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -R <#> <rout1> ... <rout#>\n");
   fprintf(stderr, "      routs: getrf, potrf, geqrf\n");
   fprintf(stderr, "   -F <mflop> : force <mflops> of timed computation\n");
   fprintf(stderr, "   -r m/n # : restrict given dimension to #\n");
   fprintf(stderr, "   -N # : largest problem to tune for is #\n");
   fprintf(stderr, "   -# <#> : repeat each timing # times\n");
   fprintf(stderr,
           "   -s # : stop repeating timing when # seconds are exceeded\n");
   fprintf(stderr, "   -a <ldagap> : lda = M + <ldagap> foreach M\n");
   fprintf(stderr, "   -f <flushKB> : flush at least this mem in LRU timers\n");
   fprintf(stderr, "   -S <#> <side1> ... <side#>\n");
   fprintf(stderr, "   -U <nuplo> <up1> ... <upN> : Vals: [u,l,q,g]\n");
   fprintf(stderr,
           "   -o[f,c,b] <file>: output nb selection as F77/C/both file\n");
   exit(ierr ? ierr : -1);
}

int *GetIntList1(int ival)
/*
 * returns integer array with iarr[0] = 1, iarr[1] = ival
 */
{
   int *iarr;
   iarr = malloc(2*sizeof(int));
   ATL_assert(iarr);
   iarr[0] = 1;
   iarr[1] = ival;
   return(iarr);
}

int *GetIntList2(int ival1, int ival2)
/*
 * returns integer array with iarr[0] = 1, iarr[1] = ival1, ival[2] = ival2
 */
{
   int *iarr;
   iarr = malloc(3*sizeof(int));
   ATL_assert(iarr);
   iarr[0] = 1;
   iarr[1] = ival1;
   iarr[2] = ival2;
   return(iarr);
}

int *GetIntList(int nargs, char **args, int i, int nmul)
/*
 * Gets a list of integers, whose length is given by atoi(args[i])*nmul
 * list is this length+1, since 0'th location gets atoi(args[i])
 */
{
   int n, *iarr, k;

   if (++i >= nargs)
      PrintUsage(args[0], i, NULL);
   n = atoi(args[i]) * nmul;
   ATL_assert(n > 0);
   iarr = malloc(sizeof(int)*(n+1));
   ATL_assert(iarr);

   iarr[0] = n / nmul;
   for (k=0; k < n; k++)
   {
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      iarr[k+1] = atoi(args[i]);
   }
   return(iarr);
}

int *RoutNames2IntList(int nargs, char **args, int i)
{
   int n, *iarr, k;

   if (++i >= nargs)
      PrintUsage(args[0], i, NULL);
   n = atoi(args[i]);
   ATL_assert(n > 0);
   iarr = malloc(sizeof(int)*(n+1));
   ATL_assert(iarr);

   iarr[0] = n;
   for (k=0; k < n; k++)
   {
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      if (!strcmp(args[i], "getrf") || !strcmp(args[i], "GETRF"))
         iarr[k+1] = LAgetrf;
      else if (!strcmp(args[i], "potrf") || !strcmp(args[i], "POTRF"))
         iarr[k+1] = LApotrf;
      else if (!strcmp(args[i], "geqrf") || !strcmp(args[i], "GEQRF"))
         iarr[k+1] = LAgeqrf;
      else if (!strcmp(args[i], "geqlf") || !strcmp(args[i], "GEQLF"))
         iarr[k+1] = LAgeqrf;
      else if (!strcmp(args[i], "gerqf") || !strcmp(args[i], "GERQF"))
         iarr[k+1] = LAgeqrf;
      else if (!strcmp(args[i], "gelqf") || !strcmp(args[i], "GELQF"))
         iarr[k+1] = LAgeqrf;
      else
         PrintUsage(args[0], i, args[i]);
   }
   return(iarr);
}

int *IntRange2IntList(int N0, int NN, int incN)
{
   int i, n;
   int *iarr;

   for (i=N0, n=0; i <= NN; i += incN) n++;
   iarr = malloc(sizeof(int)*(n+1));
   ATL_assert(iarr);
   iarr[0] = n;
   for (i=N0, n=1 ; i <= NN; i += incN, n++)
      iarr[n] = i;
   return(iarr);
}

#include "atlas_mmparse.h"
int GetGoodNForFlopRate(int rout, double nsec, int restM, int restN)
/*
 * Reads speed of kernel code, and computes the problem size that can be
 * completed in nsec seconds
 */
{
   ATL_mmnode_t *mmp;
   double floprate, time;
   int m, n, k;
   char pre;
   double GetFlopCount(enum ATL_LAROUT rout, enum ATL_LAFLG flags,
                       int M, int N, int KL, int KU, int nb);

   #ifdef SCPLX
      pre = 'c';
   #elif defined (DCPLX)
      pre = 'z';
   #elif defined(SREAL)
      pre = 's';
   #else
      pre = 'd';
   #endif
   #ifdef ATL_USEPTHREADS
      #ifdef ATL_PTR64
         k = 20000;
      #else
         k = 8000;
      #endif
   #else
      #ifdef ATL_PTR64
         k = 8000;
      #else
         k = 2000;
      #endif
   #endif
   mmp = ReadMMFileWithPath(pre, "../blas/gemm/res", "ipgen.sum");
   if (!mmp)
   {
      printf("lansrch goes with default size due to lack of timing!\n");
      return(k);
   }
   while (mmp->next)
      mmp = KillMMNode(mmp);
   floprate = mmp->mflop[0];
   KillAllMMNodes(mmp);
   #ifdef ATL_USEPTHREADS
      floprate *= ATL_NCPU;
      printf("Parallel flop rate at %g over %u cores\n", floprate, ATL_NCPU);
   #else
      printf("Serial flop rate at %g\n", floprate);
   #endif
   if (floprate < 100)  /* nonsensical result */
      return(k);
   k = 600;
   do
   {
       k += 80;
       m = (restM) ? restM : k;
       n = (restN) ? restN : k;
       time = GetFlopCount(rout, 0, m, n, 0, 0, 1) / (floprate*1000000.0);
   }
   while (time < nsec);
   k -= 80;
/*
 * Don't let K be too huge, or extremely tiny regardless of computation
 */
   if (restM || restN)
   {
      if (k > 50000)
         k = 50000;
   }
   else if (k > 10000)
      k = 10000;
   if (k < 500)
      k = 500;
   return(k);
}
void GetFlags(int nargs, char **args, int *nreps, int *flsizeKB, int *mflop,
              int *secs, int *rout, int *ldagap, int *maxN,
              int **UPLOs, int **SDs, int *F77out, char **outnam,
              int *restM, int *restN)
{
   int *NBs=NULL, *ns=NULL, *ms=NULL, *ups=NULL, *sds=NULL, *ip;
   int i, k, n;

   *ldagap = 0;
   *flsizeKB = L2SIZE/1024;
   #ifdef ATL_USEPTHREADS
      *maxN = 4000;
   #else
      *maxN = 2000;
   #endif
   *maxN = 0;
   *rout = LAgeqrf;
   *nreps = 3;
   *secs = 5;
   *restN = *restM = 0;
   *F77out = 0;
   *outnam = NULL;
   *mflop = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'R':        /* -R # <rout1> ... <routN#>  */
         NBs = RoutNames2IntList(nargs, args, i);
         i += NBs[0] + 1;
         *rout = NBs[1];
         break;
      case 'r' : /* one dimension is restricted */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (args[i-1][0] == 'm' || args[i-1][0] == 'M')
            *restM = atoi(args[i]);
         else
            *restN = atoi(args[i]);
         break;
      case 'N':
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *maxN = atoi(args[i]);
         break;
      case 'o':   /* -o[f,c,b] <file> */
         *F77out = (args[i][2] == 'f' || args[i][2] == 'F');
         if (args[i][2] == 'f' || args[i][2] == 'F')
            *F77out = 1;
         else if (args[i][2] == 'c' || args[i][2] == 'C')
            *F77out = 0;
         else
            *F77out = 2;
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *outnam = args[i];
         break;
      case 'f':                         /* -f <flushKB> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *flsizeKB = atoi(args[i]);
         break;
      case 's':                         /* -s <secs> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *secs = atoi(args[i]);
         break;
      case 'F':                         /* -F <mflop> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *mflop = atoi(args[i]);
         break;
      case 'U':                         /* -U <nup> <u1> ... <uN>;[u,l,q,g] */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         ups = malloc(sizeof(int)*(n+1));
         ups[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'U':
            case 'u':
               ups[k+1] = LAUpper;
               break;
            case 'l':
            case 'L':
            default:
               ups[k+1] = LALower;
               break;
            }
         }
         break;
      case 'S':                         /* -S <#> <side1> ... <sideN> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         n = atoi(args[i]);
         ATL_assert(n > 0);
         sds = malloc(sizeof(int)*(n+1));
         sds[0] = n;
         for (k=0; k < n; k++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, NULL);
            switch(args[i][0])
            {
            case 'L':
            case 'l':
               sds[k+1] = LALeft;
               break;
            default:
               sds[k+1] = LARight;
               break;
            }
         }
         break;
      case 'a':                         /* -a <ldagap> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *ldagap = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
/*
 * Take default values
 */
   if (*outnam == NULL)
   {
      *F77out = 0;
      *outnam = "lanbsrch_out";
   }
   if (!(*maxN))
   {
      #ifdef ATL_USEPTHREADS
         *maxN = GetGoodNForFlopRate(*rout, 3.25, *restM, *restN);
      #else
         *maxN = GetGoodNForFlopRate(*rout, 4.50, *restM, *restN);
      #endif
      fprintf(stderr, "maxN = %d\n", *maxN);
   }
   if (NBs)
      free(NBs);
   if (!ups)
      ups = GetIntList1(LAUpper);
   if (!sds)
      sds = GetIntList1(LARight);
   *UPLOs = ups;
   *SDs = sds;
}

double GetFlopCount(enum ATL_LAROUT rout, enum ATL_LAFLG flags,
                    int M, int N, int KL, int KU, int nb)
/*
 * These numbers copied from LAPACK timer routines TIMING/[EIG,LIN]/dopla[2].f
 */
{
    double m=(M?M:N), n=N, kl=KL, ku=ku, mn;
    double adds=0.0, muls=0.0;

    if (rout & LAgetrf)
    {
       mn = (m >= n) ? n : m;   /* mn = MIN(M,N) */
       adds = mn * ( m*n-(m+n)*(mn+1.0)/2.0 + (mn+1.0)*(2.0*mn+1.0)/6.0 );
       muls = adds + mn * ( m-(mn+1.0)/2.0 );
    }
    else if (rout & LAgeqrf)
    {
       if (flags & LARight)  /* LAgeqrf || LAgeqlf */
       {
          if (M >= N)
          {
             muls = n*( ((23.0/6.0)+m+(n/2.0)) + n*(m-(n/3.0)) );
             adds = n*( (5.0/6.0) + n*((1.0/2.0) + (m-(n/3.0))) );
          }
          else
          {
             muls = m*( ((23.0/6.0) + 2.0*n - (m/2.0)) + m*(n-(m/3.0)) );
             adds = m*( (5.0/6.0) + n - (m/2.0) + m*(n-(m/3.0)) );
          }
       }
       else  /* LAgerqf || LAgelqf */
       {
          if (M >= N)
          {
             muls = n*( ((29.0/6.0) + m + n/2.0) + n*(m-n/3.0) );
             adds = n*( (5.0/6.0) + m + n*(-0.5 + (m - n/3.0)) );
          }
          else
          {
             muls = m*( ((29.0/6.0) + 2.0*n - 0.5*m) +m*(n - m/3.0) );
             adds = m*( (5.0/6.0) + 0.5*m + m*(n - m/3.0) );
          }
       }
    }
    else if (rout & LApotrf)
    {
       muls = m*( (1.0/3.0) + m*((1.0/2.0) + (m/6.0)) );
       adds = (1.0/6.0)*m*(-1.0+m*m);
    }
    #if defined (SCPLX) || defined(DCPLX)
       return(2.0*adds + 6.0*muls);
    #else
       return(adds+muls);
    #endif
}


double Time2Flops(int rout, int UPLO, int M, int N, double time)
{
   double mflop;
   mflop = GetFlopCount(rout, UPLO, M, N, 0, 0, CAN_NB);
   if (mflop > 0)
      mflop /= time*1e6;  /* translate flops & time to MFLOPS */
   return(mflop);
}

double GetTimeWithReps_LU
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int i;
/*
 * Keep setsz a multiple of TYPE size for alignment reasons.  LU only accesses
 * M*N of matrix and all of IPIV.
 */
   setsz = lda*N*ATL_sizeof +
           ((M*sizeof(int)+ATL_sizeof-1)/ATL_sizeof)*ATL_sizeof;
   setszT = M*N*ATL_sizeof + M*sizeof(int);
   mflop = GetFlopCount(LAgetrf, 0, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      Mjoin(PATL,gegen)(M, N, (TYPE*)(wrksets+i*setsz), lda, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_getrf(CblasColMajor, M, N, (TYPE*)(wrksets+i*setsz), lda,
                 (int*)(wrksets+i*setsz+lda*N*ATL_sizeof));
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}

double GetTimeWithReps_LLT
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int i;
   setsz=lda*N*ATL_sizeof;   /* matrix is entire working set of LLt */
   setszT=N*N*ATL_sizeof;    /* only touch N*N portion */
   mflop = GetFlopCount(LApotrf, Uplo, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      PosDefGen(CblasColMajor, Uplo_LA2ATL(Uplo), N,
                (TYPE*)(wrksets+i*setsz), lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_potrf(Uplo, N, (TYPE*)(wrksets+i*setsz), lda);
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}

double GetTimeWithReps_LQ
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   TYPE dtmp, dtmp1;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int wlen;            /* length of QR's workspace */
   int i;
/*
 * Figure out how much workspace is required, and allocate it
 */
   test_gelqf(CblasColMajor, M, N, &dtmp1, lda, &dtmp1, &dtmp, -1);
   wlen = dtmp;
/*
 * QR accesses matrix, Min(M,N)-length tau & workspace, but for flush purposes
 * be conservative and say it only accesses A
 */
   setsz = (lda*N + wlen + Mmin(M,N)) * ATL_sizeof;
   setszT = M*N*ATL_sizeof;

   mflop = GetFlopCount(LAgeqrf, LALeft+LALower, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      Mjoin(PATL,gegen)(M, N, (TYPE*)(wrksets+i*setsz), lda, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_gelqf(CblasColMajor, M, N, (TYPE*)(wrksets+i*setsz), lda,
                 (TYPE*)(wrksets+i*setsz+(N*lda+wlen)*ATL_sizeof),
                 (TYPE*)(wrksets+i*setsz+N*lda*ATL_sizeof), wlen);
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}

double GetTimeWithReps_RQ
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   TYPE dtmp, dtmp1;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int wlen;            /* length of QR's workspace */
   int i;
/*
 * Figure out how much workspace is required, and allocate it
 */
   test_gerqf(CblasColMajor, M, N, &dtmp1, lda, &dtmp1, &dtmp, -1);
   wlen = dtmp;
/*
 * QR accesses matrix, Min(M,N)-length tau & workspace, but for flush purposes
 * be conservative and say it only accesses A
 */
   setsz = (lda*N + wlen + Mmin(M,N)) * ATL_sizeof;
   setszT = M*N*ATL_sizeof;

   mflop = GetFlopCount(LAgeqrf, LALeft+LAUpper, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      Mjoin(PATL,gegen)(M, N, (TYPE*)(wrksets+i*setsz), lda, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_gerqf(CblasColMajor, M, N, (TYPE*)(wrksets+i*setsz), lda,
                 (TYPE*)(wrksets+i*setsz+(N*lda+wlen)*ATL_sizeof),
                 (TYPE*)(wrksets+i*setsz+N*lda*ATL_sizeof), wlen);
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}

double GetTimeWithReps_QL
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   TYPE dtmp, dtmp1;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int wlen;            /* length of QR's workspace */
   int i;
/*
 * Figure out how much workspace is required, and allocate it
 */
   test_geqlf(CblasColMajor, M, N, &dtmp1, lda, &dtmp1, &dtmp, -1);
   wlen = dtmp;
/*
 * QR accesses matrix, Min(M,N)-length tau & workspace, but for flush purposes
 * be conservative and say it only accesses A
 */
   setsz = (lda*N + wlen + Mmin(M,N)) * ATL_sizeof;
   setszT = M*N*ATL_sizeof;

   mflop = GetFlopCount(LAgeqrf, LARight+LALower, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      Mjoin(PATL,gegen)(M, N, (TYPE*)(wrksets+i*setsz), lda, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_geqlf(CblasColMajor, M, N, (TYPE*)(wrksets+i*setsz), lda,
                 (TYPE*)(wrksets+i*setsz+(N*lda+wlen)*ATL_sizeof),
                 (TYPE*)(wrksets+i*setsz+N*lda*ATL_sizeof), wlen);
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}

double GetTimeWithReps_QR
   (int mflopF, int lda, int M, int N, int nb, int Uplo, int Side, int flsizeKB)
{
   double mflop, t0, t1, drep;
   TYPE dtmp, dtmp1;
   char *wrksets;       /* working sets for kernel calls */
#ifdef TCPLX
   const int lda2 = lda+lda;
#else
   const int lda2 = lda;
#endif
   size_t setsz, setszT;   /* work set size in memory, and amnt of it touched */
   size_t nrep;            /* # of reps required to force mflopF flops */
   size_t nset;            /* # of working sets allocated */
   int wlen;            /* length of QR's workspace */
   int i;
/*
 * Figure out how much workspace is required, and allocate it
 */
   test_geqrf(CblasColMajor, M, N, &dtmp1, lda, &dtmp1, &dtmp, -1);
   wlen = dtmp;
/*
 * QR accesses matrix, Min(M,N)-length tau & workspace, but for flush purposes
 * be conservative and say it only accesses A
 */
   setsz = (lda*N + wlen + Mmin(M,N)) * ATL_sizeof;
   setszT = M*N*ATL_sizeof;

   mflop = GetFlopCount(LAgeqrf, LARight+LAUpper, M, N, 0, 0, CAN_NB);
/*
 * Cannot reuse matrices (bogus to factor an already factored matrix), so we
 * must take as our total memspace MAX(nrep,nset)*setsz
 */
   ATL_assert(mflop > 0.0);
   drep = (mflopF*1.0e6) / mflop;
   nrep = (int)(drep+0.999999);
/*
 * If cacheline flush doesn't work, then we must use this method
 */
   #if ATL_LINEFLUSH
      if (nrep < 2)
         return(-1.0);                                /* do wt normal timer */
   #else
      nrep = (nrep >= 1) ? nrep : 1;
   #endif

   nset = (flsizeKB*1024+setszT-1)/setszT;
   if (nset < nrep)
      nset = nrep;
   wrksets = malloc(nset * setsz);
   ATL_assert(wrksets);

   for (i=0; i < nset; i++)
      Mjoin(PATL,gegen)(M, N, (TYPE*)(wrksets+i*setsz), lda, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      test_geqrf(CblasColMajor, M, N, (TYPE*)(wrksets+i*setsz), lda,
                 (TYPE*)(wrksets+i*setsz+(N*lda+wlen)*ATL_sizeof),
                 (TYPE*)(wrksets+i*setsz+N*lda*ATL_sizeof), wlen);
   }
   t1 = time00();
   free(wrksets);

   return((t1-t0)/((double)nrep));
}


double GetTime(int rout, int mflopF, int lda, int M, int N, int nb, int Uplo,
               int Side, int flsizeKB)
{
#if ATL_LINEFLUSH
   FLSTRUCT *flp=NULL;
#endif
   TYPE *A, *wrk=NULL, dtmp, dtmp1, *tau=NULL;
   int *ipiv=NULL, itmp, wlen;
/*
 * If the matrix alone is larger than four times the flush size, then
 * the matrix init should be self-flushing, so we avoid the time & memory
 * waste of doing the flush.  This is critical for large problems.
 */
   int DOFLUSH = (((size_t)M)*N*ATL_sizeof < ((size_t)flsizeKB)*1024*4);
   double t0, t1;
/*
 * Call routs that force particular flop count if requested; they return -1.0
 * if one invocation will suffice to force mflopF, in which case do the timing
 * in this routine, which is simpler & doesn't require LRU & as much workspace
 * If we don't have the ability to do cacheline flushing, must use LRU rout!
 */
#if ATL_LINEFLUSH
   if (mflopF > 0)
   {
#endif
      if (rout == LApotrf)
         t1 = GetTimeWithReps_LLT(mflopF, lda, M, N, nb, Uplo, Side, flsizeKB);
      else if (rout == LAgeqrf)
      {
         if (Side == LARight)
         {
            if (Uplo == LAUpper)
               t1 = GetTimeWithReps_QR(mflopF, lda, M, N, nb, Uplo, Side,
                                       flsizeKB);
            else
               t1 = GetTimeWithReps_QL(mflopF, lda, M, N, nb, Uplo, Side,
                                       flsizeKB);
         }
         else if (Uplo == LAUpper)
            t1 = GetTimeWithReps_RQ(mflopF, lda, M, N, nb, Uplo, Side,
                                    flsizeKB);
         else
            t1 = GetTimeWithReps_LQ(mflopF, lda, M, N, nb, Uplo, Side,
                                    flsizeKB);
      }
      else
         t1 = GetTimeWithReps_LU(mflopF, lda, M, N, nb, Uplo, Side, flsizeKB);
#if ATL_LINEFLUSH == 0
      return(t1);
#else
      if (t1 >= 0.0)
         return(t1);
   }
#endif
#if ATL_LINEFLUSH != 0
/*
 * Generate operands
 */
   A = GetGE(M, N, lda);
   ATL_assert(A);
   if (DOFLUSH)
      flp = ATL_GetFlushStruct(A, N*((size_t)lda)*ATL_sizeof, NULL);
   if (rout == LApotrf)
      PosDefGen(CblasColMajor, Uplo_LA2ATL(Uplo), N, A, lda);
   else if (rout & LAgeqrf)
   {                            /* QR must allocate workspace */
      if (Side == LARight)
      {
         if (Uplo == LAUpper)
         {
            test_geqrf(CblasColMajor, M, N, A, lda, &dtmp1, &dtmp, -1);
         }
         else
         {
            test_geqlf(CblasColMajor, M, N, A, lda, &dtmp1, &dtmp, -1);
         }
      }
      else if (Uplo == LAUpper)
      {
         test_gerqf(CblasColMajor, M, N, A, lda, &dtmp1, &dtmp, -1);
      }
      else
      {
         test_gelqf(CblasColMajor, M, N, A, lda, &dtmp1, &dtmp, -1);
      }
      wlen = dtmp;
      wrk = calloc(wlen, ATL_sizeof);
      ATL_assert(wrk);
      if (DOFLUSH)
         flp = ATL_GetFlushStruct(wrk, wlen*ATL_sizeof, flp);
      itmp = (M >= N) ? M : N;
      tau = calloc(itmp, ATL_sizeof);
      if (DOFLUSH)
         flp = ATL_GetFlushStruct(tau, itmp*ATL_sizeof, flp);
   }
   else
   {
      ipiv = calloc(M, sizeof(int));
      ATL_assert(ipiv);
      if (DOFLUSH)
         flp = ATL_GetFlushStruct(ipiv, M*sizeof(int), flp);
   }
/*
 * Flush cache, and do timing
 */
   if (DOFLUSH)
      ATL_FlushAreasByCL(flp);
   if (rout == LApotrf)
   {
      t0 = time00();
      test_potrf(Uplo, N, A, lda);
      t1 = time00();
   }
   else if (rout == LAgeqrf)
   {
      if (Side == LARight)
      {
         if (Uplo == LAUpper)
         {
            t0 = time00();
            test_geqrf(CblasColMajor, M, N, A, lda, tau, wrk, wlen);
            t1 = time00();
         }
         else
         {
            t0 = time00();
            test_geqlf(CblasColMajor, M, N, A, lda, tau, wrk, wlen);
            t1 = time00();
         }
      }
      else if (Uplo == LAUpper)
      {
         t0 = time00();
         test_gerqf(CblasColMajor, M, N, A, lda, tau, wrk, wlen);
         t1 = time00();
      }
      else
      {
         t0 = time00();
         test_gelqf(CblasColMajor, M, N, A, lda, tau, wrk, wlen);
         t1 = time00();
      }
   }
   else
   {
      t0 = time00();
      test_getrf(CblasColMajor, M, N, A, lda, ipiv);
      t1 = time00();
   }
   if (tau)
      free(tau);
   if (wrk)
      free(wrk);
   if (ipiv)
      free(ipiv);
   free(A);
   if (DOFLUSH)
      ATL_KillAllFlushStructs(flp);
   return(t1 - t0);
#endif
}

static char Uplo2Char(enum ATL_LAROUT rout, enum ATL_LAFLG flags)
{
   switch (rout)
   {
   case LAgeqrf:
      if (flags & LARight)
         return('Q');
      return( (flags & LAUpper) ? 'R' : 'L');
   case LApotrf:
      return( (flags & LAUpper) ? 'U' : 'L' );
   default:
      return('G');
   }
}

static char Side2Char(enum ATL_LAROUT rout, enum ATL_LAFLG i)
{
   if (rout == LAgeqrf)
   {
      if (i & LARight) return( (i & LAUpper) ? 'R' : 'L');
      else return('Q');
   }
   return( (i & LARight) ? 'R' : 'L');
}

static int Flags2Bitmap(int Rout, int Upper, int Side)
{
   return((Rout<<8)|Upper|Side);
}
enum FLGSEL {FSUplo, FSSide, FSRout};

char *Bitmap2Char(enum FLGSEL wch, int flags)
{
   int rout;

   switch(wch)
   {
   case FSUplo:
      if (flags & LAUpper)
         return("UPPER");
      else if (flags & LALower)
         return("LOWER");
      else
         return("GEREC");
   case FSSide:
      if (flags & LARight)
         return("RIGHT");
      else if (flags & LALeft)
         return("LEFT");
      else
         return("N/A");
   case FSRout:
      rout = flags >> 8;
      if (rout & LApotrf)
         return("POTRF");
      else if (rout & LAgetrf)
         return("GETRF");
      else if (rout & LAgeqrf)
      {
         if (flags & LARight)
         {
            if (flags & LALower)
               return("GEQLF");
            else
               return("GEQRF");
         }
         else if (flags & LAUpper)
            return("GERQF");
         else
            return("GELQF");
      }
   }
   return("UNKNOWN");
}

char *NameExcludingPath(char *path)
{
   char *sp=path;
   int i;
   for (i=0; path[i]; i++)
      if (path[i] == '/')
         sp = path+i+1;
   return(sp);
}

int UseN(int N, int *Ms, int *Ns)
/*
 * RETURNS: 0 if all Ns are same value, 1 else
 */
{
   int i;
   for (i=0; i < N-1 && Ns[i] == Ns[i+1]; i++);
   return(Ns[i] != Ns[i+1]);
}

void WriteCFile(char *outnam, int N, int *flgs, int *ms, int *ns, int *nbs,
                int nc, int *NBc, int *Xover)
{
   int i, nhalf, ONM, n, nn;
   char *els="", *ln, *nameNoPath;
   FILE *nbout;

   ln = malloc(strlen(outnam)+3);
   ATL_assert(ln);
   sprintf(ln, "%s.h", outnam);
   nbout = my_fopen(ln, "w");
   free(ln);
   nameNoPath = NameExcludingPath(outnam);
   fprintf(nbout, "#ifndef %s\n\n", nameNoPath);
   fprintf(nbout, "/*\n * NB selection for %s: Side='%s', Uplo='%s'\n",
           Bitmap2Char(FSRout, flgs[0]), Bitmap2Char(FSSide, flgs[0]),
           Bitmap2Char(FSUplo, flgs[0]));
   fprintf(nbout, " * M : %d", ms[0]);
   for (i=1; i < N; i++)
      fprintf(nbout, ",%d", ms[i]);
   fprintf(nbout, "\n * N : %d", ns[0]);
   for (i=1; i < N; i++)
      fprintf(nbout, ",%d", ns[i]);
   fprintf(nbout, "\n * NB : %d", nbs[0]);
   for (i=1; i < N; i++)
      fprintf(nbout, ",%d", nbs[i]);
   fprintf(nbout, "\n */\n");

   fprintf(nbout, "#define %s(n_, nb_) \\\n", nameNoPath);
   if (nc && NBc && Xover)      /* if we've got compressed arrays, */
   {                            /* just use predefined crossover points */
      if (nc == 1)
         fprintf(nbout, "   (nb_) = %d\n", NBc[0]);
      else
      {
         fprintf(nbout, "{ \\\n");
         fprintf(nbout, "   if ((n_) < %d) (nb_) = %d; \\\n", Xover[0], NBc[0]);
         for (i=1; i < nc-1; i++)
            fprintf(nbout, "   else if ((n_) < %d) (nb_) = %d; \\\n",
                    Xover[i], NBc[i]);
         fprintf(nbout, "   else (nb_) = %d; \\\n", NBc[nc-1]);
         fprintf(nbout, "}\n\n");
      }
   }
   else
   {
      for (i=0; i < N-1 && ns[i] == ns[i+1]; i++);
      ONM = (ns[i] == ns[i+1]);  /* all N the same, base choice on M */
      if (ONM) fprintf(nbout, "{ /* n_ of this func is M of matrix! */ \\\n");
      else fprintf(nbout, "{ \\\n");
      for (i=0; i < N-1; i++)
      {
         if (!ONM) { n = ns[i]; nn = ns[i+1]; }
         else { n = ms[i]; nn = ms[i+1]; }
         nhalf = (nn - n)>>1;
         fprintf(nbout, "   %sif ((n_) < %d) (nb_) = %d; \\\n",
                 els, n+nhalf, nbs[i]);
         els = "else ";
      }
      fprintf(nbout, "   else (nb_) = %d; \\\n", nbs[i]);
      fprintf(nbout, "}\n\n");
   }
   fprintf(nbout, "\n#endif    /* end ifndef %s */\n", nameNoPath);
   my_fclose(nbout);
}

void PrintIntArr(char *name, int N, int *iarr)
{
   int i;
   if (N > 0)
   {
      fprintf(stdout, "%s = %d", name, iarr[0]);
      for (i=1; i < N; i++)
         fprintf(stdout, ", %d", iarr[i]);
      fprintf(stdout, "\n");
   }
   else
      fprintf(stderr, "%s has no entries!\n", name);
}
int *CompressDecisionTree(int N, int *Ns, int *NBs)
/*
 * Takes raw arrays from timing, and removes any duplicate NBs, and returns
 * N & NB arrays that can be directly printed
 */
{
   int *ns, *nbs, n, i, gap;
   ATL_assert(N > 1);
   ns = malloc((N+N+1)*sizeof(int));
   ATL_assert(ns);
   ns++;
   nbs = ns + N;

   for (i=0; i < N-1; i++) /* Translate Ns into boundary checks */
   {
      gap = (Ns[i+1] - Ns[i])>>1;
      ns[i] = Ns[i] + gap;
      nbs[i] = NBs[i];
   }
   ns[i] = 0;
   nbs[i] = NBs[i];
   n = 0;
   i = 1;
   while (i < N)   /* remove redundant NBs */
   {
      while (i < N && nbs[n] == NBs[i])
      {
         ns[n] = ns[i];
         i++;
      }
      if (i < N)
      {
         nbs[++n] = NBs[i];
         ns[n] = ns[i++];
      }
   }
   n++;
/*
 * Copy nbs into location ns+n, rather than ns+N
 */
   for (i=0; i < n; i++)
      ns[n+i] = nbs[i];
   ns--;
   *ns = n;
   return(ns);
}

void WriteNbFile(int F77out, char *outnam, int N, int *flgs, int *ms, int *ns,
                 int *nbs)
{
   int i;
   int *ip;
   for (i=0; i < N; i++)
      fprintf(stdout, "%s: Side='%s', Uplo='%s', M=%d, N=%d, bestNB=%d\n",
              Bitmap2Char(FSRout, flgs[i]), Bitmap2Char(FSSide, flgs[i]),
              Bitmap2Char(FSUplo, flgs[i]), ms[i], ns[i], nbs[i]);

   if (UseN(N, ms, ns))
      ip = CompressDecisionTree(N, ns, nbs);
   else
      ip = CompressDecisionTree(N, ms, nbs);
   WriteCFile(outnam, N, flgs, ms, ns, nbs, ip[0], ip+1+ip[0], ip+1);
   free(ip);
}

int GetMyReps(int N, int *nreps)
/*
 * Finds the correct nreps for this N
 */
{
   int n, i;

   n = *nreps++;
   for (i=n+n-2; i>=0; i -= 2)
   {
      if (N >= nreps[i])
         return(nreps[i+1]);
   }
   return(nreps[1]);
}


#ifdef SCPLX
   #define C_PRE 'c'
#elif defined(DCPLX)
   #define C_PRE 'z'
#elif defined(SREAL)
   #define C_PRE 's'
#else
   #define C_PRE 'd'
#endif
char *latimeName(int rout, int uplo, int side, int M, int N, int lda, int kb,
                 int nsample)
/*
 * Get the name for a lapack block timing file of form:
 *   res/<pre><rout><Sd><Up>_<M>x<N>_<lda>_<kb>_<nsample>.tim
 */
{
   char *rn, *on;
   int len, i;
   rn = Bitmap2Char(FSRout, (rout<<8)|uplo|side);
   len = 17 + strlen(rn) + NumDecDigits(M) + NumDecDigits(N) + NumDecDigits(lda)
         + NumDecDigits(kb) + NumDecDigits(nsample);
   on = malloc(len*sizeof(char));
   assert(on);
   if (rout == LApotrf)
      i = sprintf(on, "res/%c%s%c_%ux%u_%u_%u_%u.tim", C_PRE, rn,
                  (uplo==LAUpper)?'U':'L', M, N, lda, kb, nsample);
   else
      i = sprintf(on, "res/%c%s_%ux%u_%u_%u_%u.tim", C_PRE, rn, M, N, lda, kb,
                  nsample);
   assert(i < len);
   return(on);
}
double GetMultSampleTimes
(
   int verb,                    /* 0 : no output, else print mflop on line */
   int retavg,                  /* 0: return min time, else return average */
   int nsample,                 /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int rout,                    /* routine to time */
   int uplo,
   int side,
   int M, int N, int lda,       /* prob dims */
   int nb                       /* nb to use */
)
/*
 * Returns the minimum of NSAMP timing samples, where NSAMP is either nsample,
 * or nsecs seconds of execution, whichever happens first.
 */
{
   double t0, t1, time, mintime, avgtime;
   int i;
   char *fn;
   extern int ONB;
   FILE *fp;

   fn = latimeName(rout, uplo, side, M, N, lda, nb, nsample);
   fp = fopen(fn, "r");
   if (fp)
   {
      avgtime = mintime = 0.0;
      for (i=0; i < nsample; i++)
      {
         if (fscanf(fp, " %le", &time) != 1)
         {
            fclose(fp);
            if (!i)
               fp = NULL;
            break;
         }
         if (verb)
            printf("%.2f ", Time2Flops(rout, uplo, M, N, time));
         avgtime += time;
         mintime = (mintime <= time && i) ? mintime : time;
      }
      if (i == nsample)          /* don't nullify so that */
         fclose(fp);   /* timing won't be rerun */
   }
/*
 * If file didn't exist, or wasn't intelligible, run timing
 */
   if (!fp)
   {
      fp = fopen(fn, "w");
      ONB = nb;
      avgtime = mintime = 0.0;
      t0 = time00();
      for (i=0; i < nsample; i++)
      {
         time = GetTime(rout, mflopF, lda, M, N, nb, uplo, side, flshszKB);
         if (verb)
            printf("%.2f ", Time2Flops(rout, uplo, M, N, time));
         avgtime += time;
         if (mintime)
            mintime = (mintime <= time) ? mintime : time;
         else
            mintime = time;
         t1 = time00() - t0;
         fprintf(fp, "%le\n", time);
/*
 *       Time runs too long, skip extra samples
 */
         if (nsecs && t1 > nsecs)
            break;
      }
      fclose(fp);
   }
   free(fn);
   return(retavg ? avgtime/i : mintime);
}


#include Mstr(Mjoin(Mjoin(atlas_,PRE),opgen_view.h))
#include Mstr(Mjoin(Mjoin(atlas_,PRE),ipgen_view.h))

ATL_UINT fastopIdxInRange(ATL_UINT minIdx, ATL_UINT maxIdx)
/*
 * RETURNS: index of smallest SPF in opgen_view within idx range [min,max]
 */
{
   unsigned int i, idx, imin=minIdx;
   float minSPF=2e16;

   for (i=minIdx; i <= maxIdx; i++)
   {
      const unsigned int *ip;
      const unsigned int idx = ATL_VWopgen_IDXMUL(i);
      float spf;
      ip = ATL_VIEW_opgen + idx;
      spf = *((float*)ip);
      if (spf < minSPF)
      {
         minSPF = spf;
         imin = i;
      }
   }
   assert(minSPF != 2e16);
   return(imin);
}

int ipfindNB  /* RETURNS: ipgen KB to use for this problem */
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int uplo, int side,
   int M, int N, int lda,       /* prob dims */
   int minKB, int maxKB,        /* min/max KB to try with opgen range */
   double *TIMB                 /* if non-NULL OUTPUT: min time */
)
{
   ATL_UINT i, kbB = Mmin(M,N)>>1;
   double tim, timB=0.0, tim0=0.0, tim1=0.0;
   maxKB = Mmin(maxKB, kbB);
   maxKB = Mmin(maxKB, 1024);
   kbB = 0;
   if (minKB < ATL_VWipgen_BEST_KB)
   {
      for (i=0; i < ATL_VWipgen_NCASES; i++)
      {
         ATL_UINT kb;
         kb = ATL_GetVWipgenKB(i);
         if (kb+kb >= minKB)
            break;
      }
      for (; i < ATL_VWipgen_NCASES; i++)
      {
         ATL_UINT kb;
         kb = ATL_GetVWipgenKB(i);
         kb += kb;
         if (kb > maxKB)
            break;
         tim = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                                   uplo, side, M, N, lda, kb);
         printf("      nb=%d: %.2f MFLOP.\n", kb,Time2Flops(rout,uplo,M,N,tim));
         if (!kbB || tim < timB)
         {
            kbB = kb;
            timB = tim;
         }
         tim0 = tim1;
         tim1 = tim;
      }
   }
   if (3*ATL_VWipgen_BEST_KB <= maxKB)
   {
      ATL_UINT kb;
      kb=ATL_VWipgen_BEST_KB*3;
      do
      {
         kb += ATL_VWipgen_BEST_KB;
         tim = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                                   uplo, side, M, N, lda, kb);
         printf("      nb=%d: %.2f MFLOP.\n", kb,Time2Flops(rout,uplo,M,N,tim));
         if (!kbB || tim < timB)
         {
            kbB = kb;
            timB = tim;
         }
         else if (tim > tim1 && tim > tim0 && tim > timB*1.03 && tim0 != 0.0)
             break;
         tim0 = tim1;
         tim1 = tim;
      }
      while (kb <= maxKB);
   }
   if (TIMB)
      *TIMB = timB;
   return(kbB);
}

int opfindNB  /* RETURNS: opgen KB to use for this problem */
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int uplo, int side,
   int M, int N, int lda,       /* prob dims */
   int minKB, int maxKB,        /* min/max KB to try with opgen range */
   double *TIMB                 /* if non-NULL OUTPUT: min time */
)
/*
 * RETURNS: best KB for this problem size in KBs within [minIdx,maxIdx] range.
 * IDEA: scope view
 */
{
   unsigned int kb = Mmin(M,N)>>1, kbB=0, inc, i;
   double tim, tim0=0.0, tim1=0.0, timB = 0.0;
   if (!kb)
      return(1);
   maxKB = Mmin(kb, maxKB);
   maxKB = Mmin(maxKB, ATL_VWopgen_MAX_KB);
   printf("   FINDING KB [%u,%u] FOR M=%d, N=%d, LDA=%d:\n", minKB, maxKB,
          M, N, lda);

/*
 * Check KB=1,2, which are GER & GER2 respectively
 */
   if (minKB < 3)
   {
      kbB = 1;
      timB = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                                uplo, side, M, N, lda, 1);
      printf("      nb=%d: %.2f MFLOP.\n", 1, Time2Flops(rout,uplo,M,N,timB));
      tim0 = timB;
      tim1 = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                                uplo, side, M, N, lda, 2);
      printf("      nb=%d: %.2f MFLOP.\n", 2, Time2Flops(rout,uplo,M,N,tim1));
      if (tim1 < timB)
      {
         timB = tim1;
         kbB = 2;
      }
   }
   assert(maxKB > minKB);
   minKB = (minKB <= 3) ? 0 : minKB-3;
   maxKB = (maxKB > 3) ? maxKB-3 : 0;
/*
 * See how course to do initial search
 */
   inc = maxKB - minKB;
   if (inc <= 4)
      inc = 1;
   else if (inc <= 8)
      inc = 2;
   else if (inc <= 32)
      inc = 4;
   else if (inc <= 64)
      inc = 8;
   else if (inc <= 128)
      inc = 16;
   else
      inc = 32;

   for (i=minKB; i < maxKB; i += inc)
   {
      unsigned int idx;
      idx = fastopIdxInRange(i, i+inc);
      kb = i+3;
      tim = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                               uplo, side, M, N, lda, kb);
      printf("      nb=%d: %.2f MFLOP.\n", kb, Time2Flops(rout,uplo,M,N,tim));
      if (!kbB || tim < timB)
      {
         timB = tim;
         kbB = kb;
      }
/*
 *    Bail if present time noticably slower than best, and we have been going
 *    down in performance for last 3 tries
 */
      if (tim > 1.03*timB && tim > tim1 && tim > tim0)
         if (timB > 0.0 && tim1 > 0.0 && tim0 > 0.0)
            break;
      tim0 = tim1;
      tim1 = tim;
   }
/*
 * Refine inc
 */
   while (inc > 1)
   {
      maxKB = minKB + inc;
      minKB = (minKB >= inc) ? minKB-inc : 0;
      inc = (inc >= 16) ? (inc>>2) : (inc>>1);
      for (i=minKB; i < maxKB; i += inc)
      {
         unsigned int idx;
         idx = fastopIdxInRange(i, i+inc);
         kb = idx+3;
         tim = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                                  uplo, side, M, N, lda, kb);
         printf("      nb=%d: %.2f MFLOP.\n", kb,Time2Flops(rout,uplo,M,N,tim));
         if (!kbB || tim < timB)
         {
            timB = tim;
            kbB = kb;
         }
      }
   }
   if (TIMB)
      *TIMB = timB;
   assert(kbB);
   return(kbB);
}

/*
 * RETURNS: best-performing block factor, starting with best-performing
 *          opgen, and then trying any 2*KB larger than that found in ipgen
 */
ATL_UINT findMaxKB              /* RETURNS: best performing KB */
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int uplo, int side,
   int M, int N, int lda        /* prob dims */
)
{
   double timB, tim, tim0, tim1;
   ATL_UINT kbB=ATL_VWopgen_BEST_KB, kb, KB, i;
   ATL_CUINT maxKB=Mmin(M,N)>>1;

   if (kbB > (N>>1))  /* if N too small for large blockings */
      return(kbB);    /* just say search with opgen_view */
/*
 * Try best-performing outer-product KB
 */
   printf("   FINDING MAX KB FOR M=%d, N=%d, LDA=%d:\n", M, N, lda);
   timB = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                             uplo, side, M, N, lda, kbB);
   printf("      nb=%d: %.2f MFLOP.\n", kbB, Time2Flops(rout,uplo,M,N,timB));
   tim0 = tim1 = timB;
/*
 * If best inner-product KB is larger than outer-product, try it even though
 * we'll call ipgemm with outer-product formulation
 */
   if (ATL_VIEW_BEST_KB > kbB)
   {
      tim = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                               uplo, side, M, N, lda, ATL_VIEW_BEST_KB);
      printf("      nb=%d: %.2f MFLOP.\n", ATL_VIEW_BEST_KB,
             Time2Flops(rout, uplo, M, N, tim));
      if (tim < timB)
      {
         timB = tim;
         kbB = ATL_VIEW_BEST_KB;
      }
      tim1 = tim;
   }
   kb = ipfindNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side, M, N, lda,
                 ATL_VWopgen_BEST_KB, 1024, &tim);
   if (tim < timB)
   {
      kbB = kb;
      timB = tim;
   }
   if (kbB <= ATL_VWopgen_BEST_KB)
   {
      kb = opfindNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side, M, N, lda,
                    1, kbB, &tim);
      if (tim < timB)
      {
         kbB = kb;
         timB = tim;
      }
   }
   printf("   MAX KB=%u, MFLOP=%.2f\n", kbB,Time2Flops(rout, uplo, M, N, timB));
   return(kb);
}

int findNB
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int uplo, int side,
   int M, int N, int lda,       /* prob dims */
   int minKB, int maxKB,        /* smallest & largest KB to try */
   double *spdMin,              /* time(minKB) / time(bestKB) */
   double *spdMax               /* time(maxKB) / time(bestKB) */
)
/*
 * Find best KB for this exact problem size, assuming in range [minKB,maxKB]
 */
{
   int kb=0;
   double tmin, timeL=0.0, timMax, timMin;
   printf("SEARCHING KB RANGE [%u,%u] for D=(%u,%u), LD=%u\n", minKB, maxKB,
          M, N, lda);
   timMin = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                               uplo, side, M, N, lda, minKB);
   timMax = GetMultSampleTimes(0, 0, nreps, nsecs, flshszKB, mflopF, rout,
                               uplo, side, M, N, lda, minKB);
   if (minKB <= ATL_VWopgen_MAX_KB)
      kb = opfindNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side, M, N, lda,
                    minKB, maxKB, &timeL);
   if (maxKB > ATL_VWopgen_MAX_KB)
   {
      int kbR;
      kbR = ipfindNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side, M, N,
                     lda, minKB, maxKB, &tmin);
      if (kb)
      {
         if (!kbR || tmin > timeL)
            tmin = timeL;
         else
            kb = kbR;
      }
      else
         kb = kbR;
   }
   else
      tmin = timeL;
   assert(kb);
   if (spdMin)
      *spdMin = timMin / tmin;
   if (spdMax)
      *spdMax = timMax / tmin;
   printf("BEST KB=%u (%.2fMFLOP) for D=(%u,%u), LD=%u\n", kb,
          Time2Flops(rout, uplo, M, N, tmin), M, N, lda);
   return(kb);
}

int *FindAllNBs
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int uplo,
   int side,
   int restM,
   int restN,
   int ldagap,
   int iL,                      /* index in nnb to find left case */
   int iR,                      /* index in nnb to find right case */
   int *nnb                     /* array of all N/NB combos found so far */
)
/*
 * Recursively searches the dimension space for good NBs given two known
 * boundary points
 */
{
   int m, n, k, lda, i, iM, kb;
   int *tp;
   double Lperc, Rperc;

   k = (nnb[iR] - nnb[iL])>>1;  /* gap between endpoints */
   if (k < 32)                  /* if within 32, close enough */
      return(nnb);              /* return done */
   k += nnb[iL];
   m = (restM) ? restM : k;
   n = (restN) ? restN : k;
   lda = m + ldagap;
/*
 * If nnb is too short, get roughly a double-length array and continue
 */
   i = nnb[1];
   if (i+2 > nnb[0])
   {
      tp = malloc(sizeof(int)*(i+i));
      ATL_assert(tp);
      memcpy(tp, nnb, (i)*sizeof(int));
      free(nnb);
      nnb = tp;
      tp[0] = i+i;
   }
/*
 * Find best NB for a problem midway between left and right N's
 */
   iM = i;
   nnb[iM] = k;
   kb = findNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,
               m, n, lda, nnb[iL+1], nnb[iR+1], &Lperc, &Rperc);
   nnb[iM+1] = kb;
   nnb[1] = i + 2;
   printf("   M=%d, N=%d, lda=%d, NB=%d\n", m, n, lda, kb);
/*
 * If middle NB greater than left NB, and loss of using different NB is more
 * than 2%, must search space between
 */
   if (nnb[iM+1] > nnb[iL+1] && Lperc > 1.02)
      nnb = FindAllNBs(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,
                       restM, restN, ldagap, iL, iM, nnb);
/*
 * If middle NB less than right NB, and loss of using different NB is more
 * than 2%, must search space between
 */
   if (nnb[iM+1] < nnb[iR+1] && Rperc > 1.02)
      nnb = FindAllNBs(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,
                       restM, restN, ldagap, iM, iR, nnb);
   return(nnb);
}

int *FindNBByN
(
   int rout,                    /* routine to time */
   int nreps,                   /* # of times to repeat timing */
   int nsecs,                   /* max # secs to spend repeating */
   int flshszKB,                /* KB to flush */
   int mflopF,                  /* mflops to force in one timing*/
   int side,
   int uplo,
   int maxN,                    /* largest dimension to try */
   int restM,                   /* 0 : vary M wt tuning dim, else fixed size */
   int restN,                   /* 0 : vary N wt tuning dim, else fixed size */
   int ldagap                   /* lda = ldagap+M */
)
/*
 * RETURNS: integer array nnb: nnb[0] holds the array length, nnb[1] holds
 *          the number of used entries.  For all other entries, even entries
 *          hold a unique N, and the next entry (odd) is the best NB to use
 *          The Ns are not sorted (will be in recursive order)
 */
{
   int *nnb;
   int i, m, n, lda, kbmax;
   double t0, tN, percL, percR;
   m = (restM) ? restM : maxN;
   n = (restN) ? restN : maxN;
   lda = m + ldagap;
   kbmax = findMaxKB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,m,n,lda);
/*
 * Find best NB for minimum problem (dimension 4)
 */
   nnb = malloc(sizeof(int)*4000);
   nnb[0] = 4000;   /* length of array in pos 0 */
   nnb[1] = 6;      /* number of entries used so far in pos 1 */
   nnb[2] = 25;
   nnb[4] = maxN;
   nnb[5] = kbmax;
   m = (restM) ? restM : 24;
   n = (restN) ? restN : 24;
   lda = m + ldagap;
   nnb[3] = findNB(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,
                   m, n, lda, 1, kbmax, NULL, NULL);
   nnb = FindAllNBs(rout, nreps, nsecs, flshszKB, mflopF, uplo, side,
                    restM, restN, ldagap, 2, 4, nnb);
   #if 0
/*
 * Translate indices into actual NBs
 */
   for (n=nnb[1], i=3; i < n; i += 2)
   {
      int j = nnb[i];
      assert(j >= 0 && j < ATL_AMM_NCASES);
      nnb[i] = ATL_AMM_KBs[j];
   }
#endif
   return(nnb);
}

void SortNNB(int *nnb)
/*
 * Sorts the array nnb on length (N, even entries).  nnb has the form:
 * nnb[0] holds the array length, nnb[1] holds the number of used entries.
 * For all other entries, even entries hold a unique N, and the next entry
 * (odd) is the best NB to use.
 */
{
   int i, j, n, nmin, imin;

   n = nnb[1];
   for (i=2; i < n-2; i += 2)  /* simple selection sort on N */
   {
      imin = i;
      nmin = nnb[i];
      for (j=i+2; j < n; j += 2)
      {
         if (nnb[j] < nmin)
         {
            nmin = nnb[j];
            imin = j;
         }
      }
      if (imin != i)
      {
         nnb[imin] = nnb[i];
         j = nnb[imin+1];
         nnb[imin+1] = nnb[i+1];
         nnb[i] = nmin;
         nnb[i+1] = j;
      }
   }
}

void PrintNNB(int *nnb)
{
   int i, n;

   printf("         N     NB\n");
   printf("==========  =====\n");
   n = nnb[1];
   for (i=2; i < n; i += 2)
      printf("%10d %5d\n", nnb[i], nnb[i+1]);
}

int TransNNB(int rout, int restM, int restN, int side, int uplo, int *nnb,
             int **FLAGs, int **Ms, int **Ns, int **NBs)
/*
 * Translates lanbsrch's nnb array into lanbtst's style in order to call
 * common WriteNBFile
 */
{
   int n, nn, i, j, iflag;
   int *ms, *ns, *nbs, *flags;
   nn = nnb[1];
   n = (nn - 2)>>1;
   ms = malloc(sizeof(int)*n);
   ATL_assert(ms);
   ns = malloc(sizeof(int)*n);
   ATL_assert(ns);
   nbs = malloc(sizeof(int)*n);
   ATL_assert(nbs);
   flags = malloc(sizeof(int)*n);
   ATL_assert(flags);

   iflag = Flags2Bitmap(rout, uplo, side);
   for (j=0, i=2; i < nn; j++, i += 2)
   {
      ms[j] = (restM) ? restM : nnb[i];
      ns[j] = (restN) ? restN : nnb[i];
      nbs[j] = nnb[i+1];
      flags[j] = iflag;
   }

   *Ms = ms;
   *Ns = ns;
   *NBs = nbs;
   *FLAGs = flags;
   return(n);
}

int main(int nargs, char **args)
/*
 * lanbsrch is a specialization of lanbtst which tries to find good LAPACK
 * blocking factors much faster.  It does this by assuming good NBs rise
 * monotonicly with N (this is not strictly true, but shouldn't cost us
 * too much performance when wrong; lack of monotonic increase probably
 * down to particular problems sizes and cleanup), and by using recursive
 * halving to quickly search the space
 */
{
   int nsample, flsizeKB, mflopF, rout, ldagap, maxN, F77out, restM, restN;
   int nsecs, n;
   int *UPLOs, *SIDEs;
   int *nnb, *Ms, *Ns, *NBs, *flags;
   char *outnam;

   GetFlags(nargs, args, &nsample, &flsizeKB, &mflopF, &nsecs, &rout, &ldagap,
            &maxN, &UPLOs, &SIDEs, &F77out, &outnam, &restM, &restN);
   nnb = FindNBByN(rout, nsample, nsecs, flsizeKB, mflopF, SIDEs[1],
                   UPLOs[1], maxN, restM, restN, ldagap);
   if (!nnb)
      return(-1);
   SortNNB(nnb);
   PrintNNB(nnb);
   n = TransNNB(rout, restM, restN, SIDEs[1], UPLOs[1], nnb,
                &flags, &Ms, &Ns, &NBs);
   free(nnb);
   free(UPLOs);
   free(SIDEs);
   WriteNbFile(F77out, outnam, n, flags, Ms, Ns, NBs);
   free(Ms);
   free(Ns);
   free(NBs);
   free(flags);
   return(0);
}
