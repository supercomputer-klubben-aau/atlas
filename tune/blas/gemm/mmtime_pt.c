#ifdef ATL_NCPU
   #if ATL_NCPU < 2
      #undef ATL_NCPU
   #endif
#endif
#ifdef ATL_NCPU
   #define _GNU_SOURCE 1 /* what manpage says you need to get CPU_SET */
   #define __USE_GNU   1 /* what you actually need set on linuxes I've seen */
   #define __USE_XOPEN2K8/* needed to avoid undef locale_t on some systems */
   #include <sched.h>    /* must include this before pthreads */
   #include <pthread.h>
   #define dumb_prand(is_) ( 0.5 - ((double)rand_r(is_))/((double)RAND_MAX) )
#else
   #define dumb_prand(is_) ( 0.5 - ((double)rand())/((double)RAND_MAX) )
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define uint unsigned int
#ifdef ATL_SZA
   #define ATL_getszA(m_, n_, mu_, nu_, vl_) ATL_SZA
#else
   uint ATL_getszA(uint M, uint K, uint mu, uint ku, uint vlen)
   {
      uint nmu=(M+mu-1)/mu;
      return(nmu*mu*K);
   }
#endif
#ifdef ATL_SZB
   #define ATL_getszB(m_, n_, mu_, nu_, vl_) ATL_SZB
#else
   uint ATL_getszB(uint K, uint N, uint ku, uint nu, uint vlen)
   {
      uint nnu=(N+nu-1)/nu;
      return(nnu*nu*K);
   }
#endif
#ifdef ATL_SZC
   #define ATL_getszC(m_, n_, mu_, nu_, vl_) ATL_SZC
#else
   uint ATL_getszC(uint M, uint N, uint mu, uint nu, uint vlen)
   {
      uint nmu=(M+mu-1)/mu, nnu=(N+nu-1)/nu, blksz=((mu*nu+vlen-1)/vlen)*vlen;
      return(nmu*nnu*blksz);
   }
#endif
#if defined(TIME_TRMVK) || defined(TIME_TRMM)
   #include "atlas_misc.h"
#else
   double ATL_walltime(void);
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
   #define ATL_Cachelen 128
   #define ATL_AlignPtr(vp) \
      (void*) ( ATL_Cachelen + ((((size_t) (vp))>>7)<<7) )
#endif
#define CINT const int

#ifdef TIME_TRMVK
   #ifdef TIME_AMM_SM
      #include "atlas_amm.h"
      #ifdef SIDE_R_
         #define ATL_usm_alloc Mjoin(PATL,utrsmR_alloc)
      #else
         #define ATL_usm_alloc Mjoin(PATL,utrsmL_alloc)
      #endif
      #ifdef SIDE_R_
         #define SD_Right 1
      #else
         #define SD_Right 0
      #endif
      sminfo_t IPINFO;
      void ATL_UTRSM(sminfo_t*, const enum ATLAS_DIAG, ATL_CINT N, ATL_CINT R,
                     const SCALAR, const TYPE*, ATL_CSZT, TYPE*, ATL_CSZT,
                     TYPE*, TYPE*, TYPE*, TYPE*);
      #ifdef TCPLX
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         void AMM_bn(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                     TYPE*, const TYPE*, const TYPE*, const TYPE*);
         void BLK2C(const size_t, const size_t, const TYPE*, const TYPE*,
                       const TYPE*, const TYPE*, TYPE*, const size_t);
      #else
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void BLK2C(const size_t, const size_t, const TYPE*, const TYPE*,
                       const TYPE*, const TYPE*, TYPE*, const size_t);
      #endif
      void AMM_b0(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                  TYPE*, const TYPE*, const TYPE*, const TYPE*);
      void AMM_b1(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                  TYPE*, const TYPE*, const TYPE*, const TYPE*);
   #else
      #ifdef TCPLX
         #error "Complex TRSM kernel timing not yet supported!"
      #else
         #ifdef SD_Right
            #error "Side=Right not yet supported!"
         #else
            #ifdef UDIAG
               #define my_trsmk Mjoin(PATL,ktrsmLLNU_rk4)
            #else
               #define my_trsmk Mjoin(PATL,ktrsmLLN_rk4)
            #endif
         void my_trsmk(ATL_CINT M, ATL_CINT N, const SCALAR alpha,
             const TYPE *T, TYPE *B, ATL_CINT ldb, TYPE *W);
         #endif
      #endif
   #endif
#elif defined(TIME_TRMM)
      #include "atlas_amm.h"
      #ifdef SIDE_R_
         #define SD_Right 1
      #else
         #define SD_Right 0
      #endif
      #ifdef TCPLX
         #if SD_Right == 1
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         void B2BLK(const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         #else
         void A2BLK(const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*, TYPE*);
         #endif
         void BLK2C(const size_t, const size_t, const TYPE*, const TYPE*,
                       const TYPE*, const TYPE*, TYPE*, const size_t);
         void AMM_bn(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                     TYPE*, const TYPE*, const TYPE*, const TYPE*);
      #else
         #if SD_Right == 1
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void B2BLK(const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         #else
         void A2BLK(const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         #endif
         void BLK2C(const size_t, const size_t, const SCALAR, const TYPE*,
                       const SCALAR, const TYPE*, const size_t);
      #endif
         void AMM_b0(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                     TYPE*, const TYPE*, const TYPE*, const TYPE*);
         void AMM_b1(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                     TYPE*, const TYPE*, const TYPE*, const TYPE*);
         #define KMM AMM_b0
#elif defined(TIME_COPY)
   #ifndef TO_BLK
      #define TO_BLK 0
   #endif
   #ifndef COPYK
      #define COPYK ATL_USERCPMM
   #endif
   #ifdef TREAL
      #ifdef COPY_C
         #if TO_BLK
            void COPYK(const size_t, const size_t, const TYPE, const TYPE*,
                       const size_t, const TYPE, TYPE*);
         #else
            void COPYK(const size_t, const size_t, const TYPE, const TYPE*,
                       const TYPE, TYPE*, const size_t);
         #endif
      #else
         #if TO_BLK
            void COPYK(const size_t, const size_t, const TYPE,
                       const TYPE*, const size_t, TYPE*);
         #else
            void COPYK(const size_t, const size_t, const TYPE,
                       const TYPE*, TYPE*, const size_t);
         #endif
      #endif
   #else
      #ifdef COPY_C
         #if TO_BLK
            void COPYK(const size_t, const size_t, const TYPE*, const TYPE*,
                       const size_t, TYPE*, TYPE*, TYPE*);
         #else
            void COPYK(const size_t, const size_t, const TYPE*, const TYPE*,
                       const TYPE*, const TYPE*, TYPE *, const size_t);
         #endif
      #else
         #if TO_BLK
            void COPYK(const size_t, const size_t, const TYPE*,
                       const TYPE*, const size_t, TYPE*, TYPE*);
         #else
            void COPYK(const size_t, const size_t, const TYPE*,
                       const TYPE*, const TYPE*, TYPE*, const size_t);
         #endif
      #endif
   #endif
#elif defined(TCPLX)
   size_t rszA, rszB, rszC;
   void CAMM_b0(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                const TYPE *pCn);
   void CAMM_b1(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                const TYPE *pCn);
   void CAMM_bn(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
                 const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
                 const TYPE *pCn);
   #ifdef BETA0
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         #ifdef CORDER
            CAMM_b0(mblks, nblks, K, Ai, Bi, Cr, Ar, Br, Cr);
            CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_b0(mblks, nblks, K, Ar, Bi, Ci, Ai, Br, Ci);
            CAMM_b1(mblks, nblks, K, Ai, Br, Ci, pAn, pBn, pCn);
         #else
            CAMM_b0(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
            CAMM_b0(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
            CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_b1(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
         #endif
      }
   #elif defined(BETA1)
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         #ifdef CORDER
            CAMM_bn(mblks, nblks, K, Ai, Bi, Cr, Ar, Br, Cr);
            CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_b1(mblks, nblks, K, Ar, Bi, Ci, Ai, Br, Ci);
            CAMM_b1(mblks, nblks, K, Ai, Br, Ci, pAn, pBn, pCn);
         #else
            CAMM_bn(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
            CAMM_b1(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
            CAMM_bn(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_b1(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
         #endif
      }
   #else
      void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *Ai,
               const TYPE *Bi, TYPE *Ci, const TYPE *pAn, const TYPE *pBn,
               const TYPE *pCn)
      {
         extern size_t rszA, rszB, rszC;
         const TYPE *Ar=Ai+rszA, *Br=Bi+rszB;
         TYPE *Cr = Ci + rszC;
         #ifdef CORDER
            CAMM_b1(mblks, nblks, K, Ai, Bi, Cr, Ar, Br, Cr);
            CAMM_b1(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_bn(mblks, nblks, K, Ar, Bi, Ci, Ai, Br, Ci);
            CAMM_bn(mblks, nblks, K, Ai, Br, Ci, pAn, pBn, pCn);
         #else
            CAMM_b1(mblks, nblks, K, Ai, Bi, Cr, Ai, Br, Ci);
            CAMM_bn(mblks, nblks, K, Ai, Br, Ci, Ar, Br, Cr);
            CAMM_b1(mblks, nblks, K, Ar, Br, Cr, Ar, Bi, Ci);
            CAMM_bn(mblks, nblks, K, Ar, Bi, Ci, pAn, pBn, pCn);
         #endif
      }
   #endif
#else
   #ifndef KMM
      #define KMM ATL_USERMM
   #endif
   void KMM(ATL_CSZT mblks, ATL_CSZT nblks, ATL_CSZT K, const TYPE *A,
            const TYPE *B, TYPE *C, const TYPE *pAn, const TYPE *pBn,
            const TYPE *pCn);
#endif

struct kmm_struct{
   #ifdef PUTC
      TYPE *C;
      size_t nmblks, ldc;
   #endif
   #ifdef TIME_COPY
      TYPE *A;
      size_t lda;
      size_t nmblks, nnblks;
      int COLWISE, BM, BN;
   #endif
   int mb, nb, kb;                      /* C: mbxnb, At: kbxmb, B: kbXnb */
   int mu, nu, ku;                      /* needed to compute mblks/nblks */
   int movA, movB, movC;                /* which mat move in flush array? */
   size_t FLSIZE;                       /* min area to move in in bytes */
   int reps;                            /* # calls to kmm in one timing */
   int vlen;                            /* vector length */
   int iam;                             /* thread rank */
   int p;                               /* total number of threads */
   int *pids;                           /* IDs of processors for affinity */
   double mf;                           /* mflop returned by timing */
   volatile unsigned char *chkin;       /* P-len array to signal start/done */
};

/*
 * Generic copy used to generate bus traffic for outer-product timings
 * It won't actually copy things in right order, but we don't care.
 * It will also not be exact same speed as final copy, but the only way
 * to know that speed is to time all copies for all kerns, which is too expens.
 * Copy time should only matter for small problems, and should normally be
 * dominated by cache affects, which are only minorly affected by
 * implementation details.  This copies best-possible-ordering should make
 * it look like a fast copy despite lack of vectorization.
 */
#ifdef PUTC
   #ifdef TCPLX
static void gecpy(const int M, const int N, const TYPE *rC, const TYPE *iC,
                  TYPE *C, size_t ldc)
{
   register int j;

   ldc += ldc;
   for (j=0; j < N; j++, C += ldc, iC += M, rC += M)
   {
      register int i;
      for (i=0; i < M; i++)
      {
         const register int i2=i+i;
         C[i2] = rC[i];
         C[i2+1] = iC[i];
      }
   }
}
   #else
static void gecpy(const int M, const int N, const TYPE *bC, TYPE *C, size_t ldc)
{
   register int i, j;
   for (j=0; j < N; j++, C += ldc, bC += M)
      for (i=0; i < M; i++)
         C[i] = bC[i];
}
   #endif
#endif
#ifdef TIME_COPY
static double getMflops(double m, double n)
#else
static double getMflops(double m, double n, double k)
#endif
{
   #if defined(TIME_TRMVK) || defined(TIME_TRMM)
      double muls, adds;
      if (m == k)
      {
         muls = n*m*(m+1.0) / 2.0;
         adds = n*m*(m-1.0) / 2.0;
      }
      else if (n == k)
      {
         muls = m*n*(n+1.0) / 2.0;
         adds = m*n*(n-1.0) / 2.0;
      }
      else
         assert(0);
      #ifdef TCPLX
         return(1e-6*(6.0*muls + 2.0*adds));
      #else
         return(1e-6*(muls + adds));
      #endif
   #elif defined(TIME_SYRKK)
      return(1e-6*n*(n+1)*k);
   #elif defined(TIME_COPY)
      return(1e-6*m*n);
   #else
      return(1e-6*2.0*m*n*k);
   #endif
}
static size_t getSetSz(int mvA, int mvB, int mvC,
                       size_t szA, size_t szB, size_t szC, size_t *NOMVSZ)
{
   size_t nomvsz=0, setsz=0;
   if (mvA)
      setsz += szA;
   else
      nomvsz += szA;
   if (mvB)
      setsz += szB;
   else
      nomvsz += szB;
   if (mvC)
      setsz += szC;
   else
      nomvsz += szC;
   if (NOMVSZ)
      *NOMVSZ = nomvsz;
   return(setsz);
}

double GetKmmMflop
(
   CINT mb, CINT nb, CINT kb,           /* C: mbxnb, At: kbxmb, B: kbXnb */
   CINT mu, CINT nu, CINT ku, CINT vlen,
   int movA, int movB, int movC,        /* which mat move in flush array? */
   struct kmm_struct *pd,               /* problem definition */
   size_t FLSIZE,                       /* min area to move in in bytes */
   CINT reps                            /* # calls to kmm in one timing */
)
/*
 * Returns MFLOP rate of matmul kernel KMM
 */
{
   #if defined(TIME_AMM_SM) || defined(TIME_TRMM)
      #ifdef TCPLX
         TYPE one[2] = {ATL_rone, ATL_rzero};
         TYPE zero[2] = {ATL_rzero, ATL_rzero};
      #else
         #define one ATL_rone
         #define zero ATL_rzero
      #endif
   #endif
   void *vp=NULL;
   TYPE *C, *A, *B, *a, *b, *c;
   TYPE *ep, *sp;  /* extra & set ptrs */
   double t0, t1, mf;
   size_t szA, szB, szC, extra, setsz, nsets, tsz, i, j, incA, incB, incC, n;
   #ifdef PUTC
      const size_t ldc=pd->ldc, NMB=pd->nmblks;
      size_t mbcnt=0;
      TYPE *CC=pd->C+(pd->ldc SHIFT)*pd->iam, *cc=CC;
   #endif
   #ifdef TIME_COPY
      const size_t lda=pd->lda;
      size_t INCBLK, INCPAN, II, JJ, NN, MM;
      unsigned int B0, B1;
      TYPE *AA=pd->A, *aa;
      #ifdef TCPLX
         TYPE alpha[2] = {0.0, 0.0}, *beta = alpha;
      #else
         TYPE alpha=0.0, beta=0.0;
      #endif
   #else
   const TYPE alpha=1.0;
   TYPE beta=1.0;
   #endif
   CINT mblks = mb/mu, nblks = nb/nu;
   const int NOMOVE = !(movA|movB|movC);
   unsigned int seed = mb*kb + (nb<<14);

#if defined(TIME_TRMVK) || defined(TIME_TRMM)
   size_t lda, ldc;
   #ifdef TIME_AMM_SM
      void *wp;
      TYPE *diag, *L, *R, *w;
      void* ATL_usm_alloc(sminfo_t*, const int, TYPE**, TYPE**, TYPE**, TYPE**);
      wp = ATL_usm_alloc(&IPINFO, kb, &diag, &L, &R, &w);
   #endif
/*
 * special declaration for TRMM
 */
   #ifdef TIME_TRMM
      void *wp=NULL;
      TYPE *pt, *pr, *pc;
      size_t sz, szT, szR;
      int nmu = (mb+mu-1)/mu, nnu=(nb+nu-1)/nu;
      int tmb=nmu*mu, tnb=nnu*nu;
      #if SD_Right == 1
         CINT K = ((nb+ku-1)/ku)*ku;
      #else
         CINT K = ((mb+ku-1)/ku)*ku;
      #endif
      #ifdef TCPLX
         TYPE alpT[2] = {ATL_rone, ATL_rzero};
      #else
         #define alpT ATL_rone
      #endif
   #endif
   setsz = mb*kb + kb*nb;
   if (!NOMOVE && FLSIZE)
      nsets = (ATL_DivBySize(FLSIZE)+setsz-1) / setsz;
   else
      movA = movB = movC = nsets = 0;
   if (nsets < 1) nsets=1;
   i = Mmax(nsets*kb, 2000) / kb;
   i = (i > 1) ? i : 2;
   nsets = i;
#if SD_Right == 1
   lda = (i*nb)|1;
   ldc = (i*mb)|1;
   if (movA) incA = nb SHIFT; else incA = 0;
   if (movC) incC = mb SHIFT; else incC = 0;
   incB = 0;
   vp = malloc(ATL_MulBySize(lda*nb + ldc*nb + mb*nb) + ATL_Cachelen);
   A = ATL_AlignPtr(vp);
   B = A + ((nb*lda) SHIFT);
   C = B + ((mb*nb) SHIFT);
   /* initialization */
   for (i=0; i<mb*(nb SHIFT); i++) B[i] = dumb_prand(&seed);
   for (i=0; i<lda*(nb SHIFT); i++) A[i] = ATL_rzero;
   for (i=0; i<ldc*(nb SHIFT); i++) C[i] = dumb_prand(&seed);
#else
   lda = (i*mb)|1;
   ldc = (i*mb)|1;
   if (movA) incA = mb SHIFT; else incA = 0;
   if (movC) incC = mb SHIFT; else incC = 0;
   incB = 0;
   vp = malloc(ATL_MulBySize(lda*mb + ldc*nb + mb*nb) + ATL_Cachelen);
   A = ATL_AlignPtr(vp);
   B = A + ((mb*lda) SHIFT);
   C = B + ((mb*nb) SHIFT);
   /* initialization */
   for (i=0; i<mb*(nb SHIFT); i++) B[i] = dumb_prand(&seed);
   for (i=0; i<lda*(mb SHIFT); i++) A[i] = ATL_rzero;
   for (i=0; i<ldc*(nb SHIFT); i++) C[i] = dumb_prand(&seed);
#endif

   #ifndef TIME_AMM_SM
      /*incA *= mb; */
      #ifdef TIME_TRMM
/*
 *       allocate workspace for TRMM
 */
         #if SD_Right == 1
            szT = ((tnb*K+vlen-1)/vlen)*vlen;
            szR = ((tmb*K+vlen-1)/vlen)*vlen;
            szC = (((mu*nu+vlen-1)/vlen)*vlen)*nmu*nnu;
            sz = ATL_MulBySize(szT + nu*ku + szR + mu*ku + szC + (mu+mu)*nu +
                               3*ATL_Cachelen);
            if (sz < ATL_MaxMalloc)
               wp = malloc(sz);
            if (!wp)
            {
               fprintf(stderr, "ERROR: can't allocate memory!!!");
               exit(-1);
            }
            pt = ATL_AlignPtr(vp);
            pr = pt + (szT SHIFT);
            pr = ATL_AlignPtr(pr);
            pc = pr + (szR SHIFT);
            pc = ATL_AlignPtr(pc);
         #else
            szT = ((tmb*K+vlen-1)/vlen)*vlen;
            szR = ((tnb*K+vlen-1)/vlen)*vlen;
            szC = (((mu*nu+vlen-1)/vlen)*vlen)*nmu*nnu;
            sz = ATL_MulBySize(szT + mu*ku + szR + nu*ku + szC + (mu+mu)*nu +
                         3*ATL_Cachelen);
            if (sz < ATL_MaxMalloc)
            wp = malloc(sz);
            if (!wp)
            {
               fprintf(stderr, "ERROR: can't allocate memory!!!");
               exit(-1);
            }
            pt = ATL_AlignPtr(vp);
            pr = pt + (szT SHIFT);
            pr = ATL_AlignPtr(pr);
            pc = pr + (szR SHIFT);
            pc = ATL_AlignPtr(pc);
         #endif
      #endif
   #else
      {
         #if SD_Right == 1
/*            int Mjoin(PATL, GetTRSMkRNU)();
            const int TU = Mjoin(PATL,GetTRSMkRNU)(); */
            const int TU = IPINFO.nu;
            const int TUTU = TU*TU;
         #else
/*            int Mjoin(PATL, GetTRSMkLMU)(); */
            const int TU = IPINFO.mu;  /* Mjoin(PATL,GetTRSMkLMU)(); */
            const int TUTU = TU*TU;
         #endif
         TYPE *dd;
         #ifdef SIDE_R_
            #ifdef TRANSA_T_
               IPINFO.incA = lda SHIFT;  /* LT,UT */
            #else
               IPINFO.incA = 1 SHIFT;         /* UN,LN */
            #endif
         #else
            #ifdef TRANSA_T_
               IPINFO.incA = 1 SHIFT;
            #else
               IPINFO.incA = lda SHIFT;
            #endif
         #endif
         for (a=L; a!=R; a++) *a = ATL_rzero; /* zero everything */
         /* now make diag blks identity */
         for (a=L, dd=diag, i=0; a<R; i++, a+=i*(TUTU SHIFT), dd+=(TU SHIFT))
         {
            int ii;
            for (ii=0; ii<TU; ii++)
            {
               a[(ii*TU+ii) SHIFT] = ATL_rone;
               dd[ii SHIFT] = ATL_rone;
            }
         }
      }
   #endif
   for (a=A, i=0; i<nsets; i++, a+=(incA) SHIFT)
   {
      int ii;
      for (ii=0; ii < kb; ii++) a[(ii*lda+ii) SHIFT] = ATL_rone;
   }
#else
   #ifdef TIME_COPY
      if (pd->COLWISE)
      {
         NN = pd->nnblks;
         MM = pd->nmblks;
         INCBLK = mb SHIFT;
         INCPAN = lda*(nb SHIFT);
         B0 = mb;
         B1 = nb;
      }
      else
      {
         NN = pd->nmblks;
         MM = pd->nnblks;
         INCBLK = lda * (nb SHIFT);
         INCPAN = (mb SHIFT);
         B0 = nb;
         B1 = mb;
      }
      AA += pd->iam * INCPAN;
      INCPAN -= INCBLK*MM;
      II = JJ = 0;
   #endif
/*
 * Get size for each matrix, and round up to ensure we keep alignment
 */
   #ifdef TIME_COPY
      szA = szB = 0;
   #else
      szA = ATL_getszA(mb, kb, mu, ku, vlen);
      szB = ATL_getszB(kb, nb, ku, nu, vlen);
   #endif
   szC = ATL_getszC(mb, nb, mu, nu, vlen);
/*
 * Compute the setsz & extra
 * setsz: size of moving mats, extra: size of stationary mats
 */
   setsz = getSetSz(movA, movB, movC, szA, szB, szC, &extra);
   if (!NOMOVE && FLSIZE)
      nsets = (ATL_DivBySize(FLSIZE)+setsz-1) / setsz;
   else
      movA = movB = movC = nsets = 0;
   extra += (mu+mu)*nu;
   tsz = ATL_Cachelen << 1;
   tsz = ATL_DivBySize(tsz);
   tsz += nsets*setsz+extra;
   vp = malloc(2*ATL_Cachelen + ATL_MulBySize(tsz));
   assert(vp);
   sp = ATL_AlignPtr(vp);
   ep = sp + nsets*setsz;
   ep = ATL_AlignPtr(ep);
   if (movA)
   {
      A = sp;
      incA = szA SHIFT;
      sp += incA * nsets;
      sp = ATL_AlignPtr(sp);
   }
   else
   {
      A = ep;
      ep += szA SHIFT;
      ep = ATL_AlignPtr(ep);
      incA = 0;
   }
   if (movB)
   {
      B = sp;
      incB = szB SHIFT;
      sp += incB * nsets;
      sp = ATL_AlignPtr(sp);
   }
   else
   {
      B = ep;
      ep += szB SHIFT;
      ep = ATL_AlignPtr(ep);
      incB = 0;
   }
   if (movC)
   {
      C = sp;
      incC = szC SHIFT;
   }
   else
   {
      C = ep;
      incC = 0;
   }
   sp = ATL_AlignPtr(vp);
   for (i=0; i < tsz; i++)
      sp[i] = dumb_prand(&seed);
#endif

   a = A; b = B; c = C;
   t0 = ATL_walltime();
   for (j=0,i=reps; i; i--)
   {
      TYPE *an, *bn, *cn;
      if (++j != nsets)
      {
         an = a+incA;
         bn = b+incB;
         cn = c+incC;
      }
      else
      {
         j = 0;
         an = A; bn = B; cn = C;
      }
      #ifdef TIME_TRMVK
         #ifdef TIME_AMM_SM
            ATL_UTRSM(&IPINFO, AtlasNonUnit, mb, nb, one, NULL, lda, c, ldc,
                      diag, L, R, w);
         #else
            my_trsmk(mb, nb, 1.0, a, c, mb, B);
         #endif
      #elif defined(TIME_TRMM)
         #if SD_Right == 1
            #ifdef TCPLX
               A2BLK(nb, mb, alpT, b, ldc, pr, pr+szR);
               B2BLK(nb, alpT, a, lda, pt, pt+szT);
               AMM_b0(nmu, nnu, K, pr+szR, pt+szT, pc, pr,pt+szT, pc+szC);
               AMM_b0(nmu, nnu, K, pr, pt+szT, pc+szC, pr, pt, pc);
               AMM_bn(nmu, nnu, K, pr, pt, pc, pr+szR, pt, pc+szC);
               AMM_b1(nmu, nnu, K, pr+szR, pt, pc+szC, pr+szR, pt+szT, pc);
               BLK2C(mb, nb, one, pc, pc+szC, zero, b, ldc);
            #else
               A2BLK(nb, mb, alpT, b, ldc, pr);
               B2BLK(nb, alpT, a, lda, pt);
               AMM_b0(nmu, nnu, K, pr, pt, pc, pr, pt, pc);
               BLK2C(mb, nb, one, pc, zero, b, ldc);
            #endif
         #else
            #ifdef TCPLX
               A2BLK(mb, alpT, a, lda, pt, pt+szT);
               B2BLK(mb, nb, alpT, b, ldc, pr, pr+szR);
               AMM_b0(nmu, nnu, K, pt+szT, pr+szR, pc, pt, pr+szR, pc+szC);
               AMM_b0(nmu, nnu, K, pt, pr+szR, pc+szC, pt, pr, pc);
               AMM_bn(nmu, nnu, K, pt, pr, pc, pt+szT, pr, pc+szC);
               AMM_b1(nmu, nnu, K, pt+szT, pr, pc+szC, pt+szT, pr+szR, pc);
               BLK2C(mb, nb, one, pc, pc+szC, zero, b, ldc);
            #else
               A2BLK(mb, alpT, a, lda, pt);
               B2BLK(mb, nb, alpT, b, ldc, pr);
               AMM_b0(nmu, nnu, K, pt, pr, pc, pt, pr, pc);
               BLK2C(mb, nb, one, pc, zero, b, ldc);
            #endif
         #endif
      #elif defined(TIME_SYRKK)
         KMM(mblks, nblks, kb, a, a, c, an, an, cn);
      #elif defined(TIME_COPY)
         #ifdef COPY_C
            #ifdef TCPLX
               #if TO_BLK
                  COPYK(mb, nb, alpha, AA, lda, beta, c+szC, c);
               #else
                  COPYK(mb, nb, alpha, c+szC, c, beta, AA, lda);
               #endif
            #else
               #if TO_BLK
                  COPYK(mb, nb, alpha, AA, lda, beta, c);
               #else
                  COPYK(mb, nb, alpha, c, beta, AA, lda);
               #endif
            #endif
         #else  /* timing A or B copy */
            #ifdef TREAL
               #if TO_BLK
                  COPYK(B0, B1, alpha, AA, lda, c);
               #else
                  COPYK(B0, B1, alpha, c, AA, lda);
               #endif
            #else
               #if TO_BLK
                  COPYK(B0, B1, alpha, AA, lda, c+szC, c);
               #else
                  COPYK(B0, B1, alpha, c+szC, c, AA, lda);
               #endif
            #endif
         #endif
         aa += INCBLK;
         if (++II == MM)
         {
            II = 0;
            aa += INCPAN;
            if (++JJ == NN)
            {
               JJ = 0;
               aa = AA;
            }
         }
      #else
         #ifdef ALLTRANS_
            KMM(nblks, mblks, kb, b, a, c, bn, an, cn);
         #else
            KMM(mblks, nblks, kb, a, b, c, an, bn, cn);
         #endif
         #ifdef PUTC
            #ifdef TCPLX
               gecpy(mb, nb, c+szC, c, cc, ldc);
            #else
               gecpy(mb, nb, c, cc, ldc);
            #endif
            if (++mbcnt == NMB)
            {
               cc = CC;
               mbcnt = 0;
            }
            else
               cc += mb SHIFT;
         #endif
      #endif
      a = an;
      b = bn;
      c = cn;
   }
   t1 = ATL_walltime() - t0;
   #ifdef TIME_COPY
      #ifdef COPY_C
         mf = reps*getMflops(mb, nb) / t1;
      #else
         mf = reps*getMflops(nb, kb) / t1;  /* need to dist mb/nb eventually */
      #endif
      #ifdef TCPLX
         mf *= 2.0;
      #endif
   #else
      mf = reps*getMflops(mb, nb, kb) / t1;
      #if defined(TCPLX) && !defined(TIME_TRMVK)
         mf *= 4.0;
      #endif
   #endif
   #if (defined(TIME_TRMVK) && defined(TIME_AMM_SM)) || defined(TIME_TRMM)
      free(wp);
   #endif
   free(vp);
   return(mf);
}

#ifdef PRINT_COREID
   #include <utmpx.h>
#endif
/*#define PRINT_NUMAIDS */
#ifdef PRINT_NUMAIDS
   #define _GNU_SOURCE 1
   #include <unistd.h>
   #include <sys/syscall.h>
#endif
void *TimeOnCore(void *vp)
{
   struct kmm_struct *kp = vp;
   const int P = kp->p;
   int i;
   volatile unsigned char *chkin = kp->chkin;
   #ifdef PRINT_COREID
      printf("core=%d\n", sched_getcpu());
   #endif
#ifdef PRINT_NUMAIDS
    unsigned cpu, node;
    syscall(SYS_getcpu, &cpu, &node, NULL);
    printf("cpu=%u, node=%u\n", cpu, node);
#endif
/*
 * First we barrier, so that all cores are active.  Otherwise, 1st core to
 * start may run with bus to himself, and not give us a measure of true
 * parallel performance.  Want full-on contention as in perfect parallel code.
 * We wait on chkin array to have all non-zero entries.  Even on weakly-ordered
 * caches this should work, though the delay may be long.
 */
   #ifdef ATL_NCPU
      chkin[kp->iam] = 1;
      for (i=0; i < P; i++)
         while(!chkin[i]);
   #endif
   kp->mf = GetKmmMflop(kp->mb, kp->nb, kp->kb, kp->mu, kp->nu, kp->ku,
                        kp->vlen, kp->movA, kp->movB, kp->movC, kp,
                        kp->FLSIZE, kp->reps);
   return(NULL);
}


#ifndef ATL_NCPU
double *TimeOnCores(struct kmm_struct *kb)
{
   double *mflops;
   int i, p;

   mflops = malloc(sizeof(double));
   assert(mflops);
   p = kb->p;
   kb->iam = 0;
   TimeOnCore(kb);
   free(kb->pids);
   mflops[0] = kb->mf;
   return(mflops);
}
#else
double *TimeOnCores(struct kmm_struct *kb)
{
   struct kmm_struct *kp;
   pthread_t *threads;
   pthread_attr_t *attr;
   cpu_set_t *cpuset;
   double *mflops;
   int i, p;
   unsigned char *chkin;

   p = kb->p;
/*
 * On PowerPC Linux, pthread_attr_setaffinity_np sometimes attempts to
 * realloc the cpuset variable, so it must be malloced not taken from
 * the stack.  This is crazy behaviour, but it is what happens.
 */
   cpuset = malloc(sizeof(cpu_set_t));
   kp = malloc(sizeof(struct kmm_struct)*p);
   threads = malloc(sizeof(pthread_t)*p);
   attr = malloc(sizeof(pthread_attr_t)*p);
   mflops = malloc(sizeof(double)*p);
   chkin = malloc(sizeof(char)*p);
   assert(cpuset && kp && threads && attr && mflops && chkin);
   for (i=0; i < p; i++)   /* init chkin to 0 before starting any threads */
      chkin[i] = 0;        /* when all entries non-zero, all thrds started */
   for (i=0; i < p; i++)
   {
      memcpy(kp+i, kb, sizeof(struct kmm_struct));
      kp[i].chkin = (volatile char*)chkin;
      kp[i].iam = i;
      CPU_ZERO(cpuset);
      CPU_SET(kp->pids[i], cpuset);
      assert(!pthread_attr_init(attr+i));
      assert(!pthread_attr_setaffinity_np(attr+i, sizeof(cpu_set_t), cpuset));
      assert(!pthread_attr_setdetachstate(attr+i, PTHREAD_CREATE_JOINABLE));
      assert(!pthread_create(threads+i, attr+i, TimeOnCore, kp+i));
   }
   for (i=0; i < p; i++)
   {
      pthread_join(threads[i], NULL);
      pthread_attr_destroy(attr+i);
      mflops[i] = kp[i].mf;
   }
   free(kp->pids);
   free(kp);
   free(threads);
   free(attr);
   free(cpuset);
   free(chkin);
   return(mflops);
}
#endif

void GetStat(int n, double *d, double *min, double *max, double *avg)
{
   int i;
   double dmin, dmax, dsum;

   dmin = dmax = dsum = d[0];
   for (i=1; i < n; i++)
   {
      dmax = (dmax >= d[i]) ? dmax : d[i];
      dmin = (dmin <= d[i]) ? dmin : d[i];
      dsum += d[i];
   }
   *min = dmin;
   *max = dmax;
   *avg = dsum / (double)n;
}

void PrintUsage(char *name, int iarg, char *arg)
{
   fprintf(stderr, "\nERROR around arg %d (%s).\n", iarg, arg ? arg:"unknown");
   fprintf(stderr, "USAGE: %s [flags], where flags are:\n", name);
   fprintf(stderr, "   -p <#> : use # threads (with affinity)\n");
   fprintf(stderr, "   -tl <#> id1 ... id# : spawn # threads to given IDs\n");
   fprintf(stderr, "   -B <#> : mb = nb = kb = #\n");
   fprintf(stderr, "   -m <#> : mb = #\n");
   fprintf(stderr, "   -n <#> : nb = #\n");
   fprintf(stderr, "   -k <#> : kb = #\n");
   fprintf(stderr, "   -V <veclen>\n");
   #ifdef TIME_COPY
      fprintf(stderr, "   -D <M> <lda> <N> <COLWISE>: copy only\n");
      #ifndef COPY_C
      fprintf(stderr, "   -A [a,b]: for A/B copy, set which you are copying\n");
      #endif
   #endif
   fprintf(stderr, "   -u[mnk] <#> : M/N/K loop unrolling is #\n");
   fprintf(stderr, "   -r <#> : set the # of times to call KMM\n");
   fprintf(stderr, "   -R <mf>: set # reps to force <mf> MFLOPs\n");
   fprintf(stderr, "   -F <kb> : set flush size in kilobytes\n");
   fprintf(stderr, "   -M[a,b,c] <#> : mov[A,B,C] = #\n");
   exit(iarg ? iarg : -1);
}

struct kmm_struct *GetFlags(int nargs, char **args, FILE **fpout)
{
   FILE *fp;
   struct kmm_struct *kp;
   double mflops=750.0;
   int i, j, ACPY=1, IGMF=0;  /* ignore mflops.frc file? */

   *fpout = NULL;
   kp = malloc(sizeof(struct kmm_struct));
   assert(kp);
   #ifdef TIME_COPY
      kp->nmblks = kp->nnblks = 2000;
      kp->lda = kp->nmblks | 1;  /* make odd to force align probs */
      #ifdef COPY_C
         kp->COLWISE = 1;
      #else
         kp->COLWISE = 0;
      #endif
   #endif
   kp->pids = NULL;
   kp->p = 1;
   kp->mb = kp->nb = kp->kb = 40;
   kp->mu = kp->nu = 4;
   kp->ku = 1;
   kp->movA = kp->movB = kp->movC = 0;
   kp->FLSIZE = L2SIZE;
   kp->reps = 0;
   kp->vlen = 1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      #if defined(TIME_COPY) && !defined(COPY_C)
      case 'A':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         ACPY = (args[i][0] != 'b' && args[i][0] != 'B');
         break;
      #endif
      case 'f':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         *fpout = fopen(args[i], "w");
         break;
      case 'V':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->vlen = atoi(args[i]);
         break;
      #ifdef TIME_COPY
      case 'D':  /* -D <M> <lda> <N> <COL> */
         if (i+4 >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->nmblks = atol(args[i+1]);
         kp->lda = atol(args[i+2]);
         kp->nnblks = atol(args[i+3]);
         kp->COLWISE = atoi(args[i+4]);
         i += 4;
         break;
      #endif
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->FLSIZE = atol(args[i]) * 1024;
         break;
      case 'u':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         j = atoi(args[i]);
         if (args[i-1][2] == 'k')
            kp->ku = j;
         else if (args[i-1][2] == 'n')
            kp->nu = j;
         else
            kp->mu = j;
         break;
      case 'R':
         if (args[i][2] == 'f')
            IGMF=1;
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         mflops = atof(args[i]);
         break;
      case 'r':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->reps = atoi(args[i]);
         break;
      case 't':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->p = atoi(args[i]);
         kp->pids = malloc(sizeof(int)*kp->p);
         assert(kp->pids);
         for (j=0; j < kp->p; j++)
         {
            if (++i >= nargs)
               PrintUsage(args[0], i, "out of arguments");
            kp->pids[j] = atoi(args[i]);
         }
         break;
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->p = atoi(args[i]);
         break;
      case 'm':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->mb = atoi(args[i]);
         break;
      case 'n':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->nb = atoi(args[i]);
         break;
      case 'k':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->kb = atoi(args[i]);
         break;
      case 'B':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         kp->mb = kp->nb = kp->kb = atoi(args[i]);
         break;
      case 'M':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         switch(args[i-1][2])
         {
         case 'c':
         case 'C':
            kp->movC = atoi(args[i]);
            break;
         case 'b':
         case 'B':
            kp->movB = atoi(args[i]);
            break;
         case 'a':
         case 'A':
            kp->movA = atoi(args[i]);
            break;
         default:
            PrintUsage(args[0], i-1, "unknown mov matrix");
         }
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
/*
 * If there is a tuned <upr>mflops.frc file, use this in preference to the
 * commandline argument unless the -Rf override flag was used
 */
   if (!IGMF)
   {
      char fn[16] = {'s','m','f','l','o','p','s','.','f','r','c','\0'};
      #if defined(DREAL) || defined(DCPLX)
         fn[0] = 'd';
      #endif
      fp = fopen(fn, "r");
      if (fp)
      {
         assert(fscanf(fp, "%le", &mflops) == 1);
         fclose(fp);
      }
      #if defined(TIME_TRMVK) && defined(TIME_AMM_SM)
         mflops *= 0.25;
      #endif
   }
   if (!kp->reps)
   {
      kp->reps = (mflops*1000000.0/((2.0*kp->mb)*kp->nb*kp->kb));
      if (kp->reps < 1)
         kp->reps = 1;
   }
   if (!kp->pids)
   {
      kp->pids = malloc(sizeof(int)*kp->p);
      assert(kp->pids);
      #ifdef ATL_ARCH_XeonPHI
      {
         int n4 = ((kp->p)>>2)<<2, nr = kp->p - n4;
         for (j=0; j < n4; j += 2)
         {
            kp->pids[j] = 2*j;
            kp->pids[j+1] = 2*j+1;
         }
         switch(nr)
         {
         case 3:
            kp->pids[j+2] = 2*(j+2);
         case 2:
            kp->pids[j+1] = 2*j+1;
         case 1:
            kp->pids[j] = 2*j;
            break;
         case 0:;
         }
      }
      #else
         for (j=0; j < kp->p; j++)
             kp->pids[j] = j;
      #endif
   }
   #if defined(TIME_TRMVK) && defined(TIME_AMM_SM)
      IPINFO.mu = kp->mu;
      IPINFO.nu = kp->nu;
   #endif
   #ifdef TIME_COPY  /* Reduce M/N to nmblks, nnblks, handle P */
   {
      const size_t gap = kp->lda - kp->nmblks;
      unsigned int nb;

      kp->BM = kp->mb;
      kp->BN = kp->nb;
      nb = kp->BM;
      kp->nmblks = (kp->nmblks+nb-1)/nb;
      if (!kp->COLWISE)
         kp->nmblks = (kp->nmblks >= kp->p) ? kp->nmblks : kp->p;
      kp->lda = kp->nmblks*nb + gap;
      nb = kp->BN;
      kp->nnblks = (kp->nnblks+nb-1)/nb;
      if (kp->COLWISE)
         kp->nnblks = (kp->nnblks >= kp->p) ? kp->nnblks : kp->p;
      kp->movA = kp->movB = 0;
      kp->movC = 1;
   }
   #endif
   #ifdef TIME_SYRKK
      kp->movB = 0;
   #endif
   #ifdef PUTC
      i = 2000 / kp->mb;
      i = (i > 1) ? i : 2;
      kp->nmblks = i;
      kp->ldc = (i*kp->mb)|1;
      i = kp->p * kp->nb;
      kp->C = malloc(sizeof(TYPE)*kp->ldc*(i SHIFT));
      assert(kp->C);
   #endif
   return(kp);
}

int main(int nargs, char **args)
{
   struct kmm_struct *kp;
   int i, p;
   double *dp;
   double min, max, avg;
   FILE *fpout = stdout;

   kp = GetFlags(nargs, args, &fpout);
   p = kp->p;
   #ifndef ATL_NCPU
      assert(p < 2);
   #endif
   #if defined(TIME_TRMVK) && defined(TIME_AMM_SM)
      IPINFO.amm_b0 = AMM_b0;
      IPINFO.amm_b1 = AMM_b1;
      #ifdef TCPLX
         IPINFO.amm_bn = AMM_bn;
      #endif
      IPINFO.b2blk = B2BLK;
      IPINFO.a2blk = A2BLK;
   #endif
   #ifdef TIME_COPY
   {
      TYPE *A;
      size_t k, N;
      unsigned int seed;

      N = kp->lda * kp->nnblks * kp->nb;
      seed = kp->lda+(kp->nnblks<<8);
      kp->A = A = malloc(ATL_MulBySize(N));
      assert(A);
      for (k=0; k < N; k++)
         A[k] = dumb_prand(&seed);
   }
   #endif
   dp = TimeOnCores(kp);
   #ifdef PUTC
      free(kp->C);
   #endif
   free(kp);
   GetStat(p, dp, &min, &max, &avg);
   printf("PER-CORE: %le", dp[0]);
   for (i=1; i < p; i++)
      printf(", %le", dp[i]);
   printf("\nALL CORES: min=%.2f, max=%.2f, avg=%.2f\n", min, max, avg);
   if (fpout)
   {
      fprintf(fpout, "%d 1\n", p);
      for (i=0; i < p; i++)
         fprintf(fpout, "%e\n", dp[i]);
   }
   free(dp);
   if (fpout && fpout != stdout && fpout != stderr)
      fclose(fpout);
   exit(0);
}
