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
#ifndef ATL_SZA
   uint ATL_getszA(uint M, uint K, uint mu, uint ku, uint vlen)
   {
      uint nmu=(M+mu-1)/mu;
      return(nmu*mu*K);
   }
#endif
#ifndef ATL_SZB
   uint ATL_getszB(uint K, uint N, uint ku, uint nu, uint vlen)
   {
      uint nnu=(N+nu-1)/nu;
      return(nnu*nu*K);
   }
#endif
#ifndef ATL_SZC
   uint ATL_getszC(uint M, uint N, uint mu, uint nu, uint vlen)
   {
      uint nmu=(M+mu-1)/mu, nnu=(N+nu-1)/nu, blksz=((mu*nu+vlen-1)/vlen)*vlen;
      return(nmu*nnu*blksz);
   }
#endif
#ifdef TIME_TRMVK
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
      void ATL_UTRSM(sminfo_t *ip, const enum ATLAS_DIAG Diag, ATL_CINT N,
                     ATL_CINT R, const SCALAR alpha, const TYPE *A,
                     ATL_CINT lda, TYPE *X, ATL_CINT ldx, TYPE *diag,
                     TYPE *L, TYPE *RW, TYPE *w);
      #ifdef TCPLX
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
      #else
         void A2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
         void B2BLK(const size_t, const size_t, const SCALAR,
                    const TYPE*, const size_t, TYPE*);
      #endif
      void AMM_b0(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                  TYPE*, const TYPE*, const TYPE*, const TYPE*);
      void AMM_b1(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                  TYPE*, const TYPE*, const TYPE*, const TYPE*);
      void AMM_bn(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                  TYPE*, const TYPE*, const TYPE*, const TYPE*);
      #if defined(SDN) && (SDN == 1)
         #define SD_Right 1
      #else
         #define SD_Right 0
      #endif
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
   int mb, nb, kb;                      /* C: mbxnb, At: kbxmb, B: kbXnb */
   int mu, nu, ku;                      /* needed to compute mblks/nblks */
   int movA, movB, movC;                /* which mat move in flush array? */
   long FLSIZE;                       /* min area to move in in bytes */
   long flA, flB, flC;                /* individual flush sizes in bytes*/
   int reps;                            /* # calls to kmm in one timing */
   int vlen;                            /* vector length */
   int iam;                             /* thread rank */
   int p;                               /* total number of threads */
   int *pids;                           /* IDs of processors for affinity */
   double mf;                           /* mflop returned by timing */
   volatile unsigned char *chkin;       /* P-len array to signal start/done */
};

static double getMflops(double m, double n, double k)
{
   #ifdef TIME_TRMVK
      #if SD_Right == 1
         double muls = n*m*(m+1.0) / 2.0;
         double adds = n*m*(m-1.0) / 2.0;
      #else
         double muls = m*n*(n+1.0) / 2.0;
         double adds = m*n*(n-1.0) / 2.0;
      #endif
      #ifdef TCPLX
         return(1e-6*(6.0*muls + 2.0*adds));
      #else
         return(1e-6*(muls + adds));
      #endif
   #elif defined(TIME_SYRKK)
      return(1e-6*n*(n+1)*k);
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

#if 1

double GetKmmMflop
(
   CINT mb, CINT nb, CINT kb,           /* C: mbxnb, At: kbxmb, B: kbXnb */
   CINT mu, CINT nu, CINT ku, CINT vlen,
   int movA, int movB, int movC,        /* which mat move in flush array? */
   long FLSIZE,                         /* min area to move in in bytes */
   long flA, long flB, long flC,        /* min area to move in in bytes */
   CINT reps                            /* # calls to kmm in one timing */
)
/*
 * Returns MFLOP rate of matmul kernel KMM
 */
{
   #ifdef TCPLX
      TYPE one[2] = {ATL_rone, ATL_rzero};
   #else
      #define one ATL_rone
   #endif
   CINT mblks = mb/mu, nblks = nb/nu;
   const int NOMOVE = !(movA|movB|movC);
   size_t szA, szB, szC, extra, setsz, nsets, tsz, i, j, incA, incB, incC, n;
   TYPE *C, *A, *B, *a, *b, *c;
   TYPE *ep, *sp;  /* extra & set ptrs */
   double t0, t1, mf;
   const TYPE alpha=1.0;
   TYPE beta=1.0;
   void *vp=NULL;
   unsigned int seed = mb*kb + (nb<<14);
   size_t nseta, nsetb, nsetc, ja, jb, jc;

#ifdef TIME_TRMVK
   size_t lda, ldc;
   #ifdef TIME_AMM_SM
      void *wp;
      TYPE *diag, *L, *R, *w;
      void* ATL_usm_alloc
         (sminfo_t *ip, const int, TYPE**, TYPE**, TYPE**, TYPE**);
      wp = ATL_usm_alloc(&IPINFO, kb, &diag, &L, &R, &w);
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
   #else
      {
         #if SD_Right == 1
            int Mjoin(PATL, GetTRSMkRNU)();
            const int TU = Mjoin(PATL,GetTRSMkRNU)();
            const int TUTU = TU*TU;
         #else
            int Mjoin(PATL, GetTRSMkLMU)();
            const int TU = Mjoin(PATL,GetTRSMkLMU)();
            const int TUTU = TU*TU;
         #endif
         TYPE *dd;
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
/*
 * Get size for each matrix, and round up to ensure we keep alignment
 */
   szA = ATL_getszA(mb, kb, mu, ku, vlen);
   szB = ATL_getszB(kb, nb, ku, nu, vlen);
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
   if (flA > 0 || flB > 0 || flC > 0)
   {
      nseta = (ATL_DivBySize(flA)+szA-1) / szA;
      nsetb = (ATL_DivBySize(flA)+szB-1) / szB;
      nsetc = (ATL_DivBySize(flA)+szC-1) / szC;
      if (!movA || nseta < 1) nseta = 1;
      if (!movB || nsetb < 1) nsetb = 1;
      if (!movC || nsetc < 1) nsetc = 1;

      tsz += nseta*szA + nsetb*szB + nsetc*szC + extra;
      vp = malloc(6*ATL_Cachelen + ATL_MulBySize(tsz));
   }
   else
   {
      nseta = nsetb = nsetc = nsets;
      if (!movA || nseta < 1) nseta = 1;
      if (!movB || nsetb < 1) nsetb = 1;
      if (!movC || nsetc < 1) nsetc = 1;
      tsz += nsets*setsz+extra;
      vp = malloc(8*ATL_Cachelen + ATL_MulBySize(tsz));
   }
   assert(vp);
   sp = ATL_AlignPtr(vp);
   if (flA > 0 || flB > 0 || flC > 0)
      ep = sp + nseta*szA + nsetb*szB + nsetc*szC;
   else
      ep = sp + nsets*setsz;
   ep += ATL_DivBySize(3*ATL_Cachelen);
   ep = ATL_AlignPtr(ep);
   if (movA)
   {
      A = sp;
      incA = szA SHIFT;
      sp += incA * nseta;
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
      sp += incB * nsetb;
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
   /* OLD Init *//*
   sp = ATL_AlignPtr(vp);
   for (i=0; i < (tsz+ATL_DivBySize(5*ATL_Cachelen)); i++)
      sp[i] = dumb_prand(&seed);
   */
   for (i=0; i<nseta*szA; i++)
      A[i] = dumb_prand(&seed);
   for (i=0; i<nsetb*szB; i++)
      B[i] = dumb_prand(&seed);
   for (i=0; i<nsetc*szC; i++)
      C[i] = dumb_prand(&seed);
#endif

   a = A; b = B; c = C;
   t0 = ATL_walltime();
   for (ja=jb=jc=0,i=reps; i; i--)
   {
      TYPE *an, *bn, *cn;
      if (++ja != nseta)
         an = a+incA;
      else
      {
         ja = 0; an = A;
      }
      if (++jb != nsetb)
         bn = b+incB;
      else
      {
         jb = 0; bn = B;
      }
      if (++jc != nsetc)
         cn = c+incC;
      else
      {
         jc = 0; cn = C;
      }
      #ifdef TIME_TRMVK
         #ifdef TIME_AMM_SM
            my_trsmk(AtlasNonUnit, mb, nb, one, NULL, lda, c, ldc,
                     diag, L, R, w);
         #else
            my_trsmk(mb, nb, 1.0, a, c, mb, b);
         #endif
      #elif defined(TIME_SYRKK)
         KMM(mblks, nblks, kb, a, a, c, an, an, cn);
      #else
         KMM(mblks, nblks, kb, a, b, c, an, bn, cn);
      #endif
      a = an;
      b = bn;
      c = cn;
   }
   t1 = ATL_walltime() - t0;
   mf = reps*getMflops(mb, nb, kb) / t1;
   #if defined(TCPLX) && !defined(TIME_TRMVK)
      mf *= 4.0;
   #endif
   #if defined(TIME_TRMVK) && defined(TIME_AMM_SM)
      free(wp);
   #endif
   free(vp);
   return(mf);
}

#else

double GetKmmMflop
(
   CINT mb, CINT nb, CINT kb,           /* C: mbxnb, At: kbxmb, B: kbXnb */
   CINT mu, CINT nu, CINT ku, CINT vlen,
   int movA, int movB, int movC,        /* which mat move in flush array? */
   long FLSIZE,                       /* min area to move in in bytes */
   CINT reps                            /* # calls to kmm in one timing */
)
/*
 * Returns MFLOP rate of matmul kernel KMM
 */
{
   CINT mblks = mb/mu, nblks = nb/nu;
   const int NOMOVE = !(movA|movB|movC);
   size_t szA, szB, szC, extra, setsz, nsets, tsz, i, j, incA, incB, incC, n;
   TYPE *C, *A, *B, *a, *b, *c;
   TYPE *ep, *sp;  /* extra & set ptrs */
   double t0, t1, mf;
   const TYPE alpha=1.0;
   TYPE beta=1.0;
   void *vp=NULL;
   unsigned int seed = mb*kb + (nb<<14);
/*
 * Get size for each matrix, and round up to ensure we keep alignment
 */
   szA = ATL_getszA(mb, kb, mu, ku, vlen);
   szB = ATL_getszB(kb, nb, ku, nu, vlen);
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
/*
 * For copy timing, barrier here so above init traffic not charged to cpytime
 * Kernels usually have some copying to be done, so its OK to overlap a little
 * init mem traffic!
 */
   #if defined(ATL_NCPU) && defined(TIME_COPY)
      chkin[kp->iam] = 2;
      for (i=0; i < P; i++)
         while(chkin[i] < 2);
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
         my_trsmk(mb, nb, 1.0, a, c, mb, b);
      #elif defined(TIME_SYRKK)
         KMM(mblks, nblks, kb, a, a, c, an, an, cn);
      #else
         KMM(mblks, nblks, kb, a, b, c, an, bn, cn);
      #endif
      a = an;
      b = bn;
      c = cn;
   }
   t1 = ATL_walltime() - t0;
   mf = reps*getMflops(mb, nb, kb) / t1;
   #ifdef TCPLX
      mf *= 4.0;
   #endif
   free(vp);
   return(mf);
}

#endif

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
                        kp->vlen, kp->movA, kp->movB, kp->movC,
                        kp->FLSIZE, kp->flA, kp->flB, kp->flC, kp->reps);
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
      mflops[i] = kp[i].mf;
   }
   free(cpuset);
   free(kp->pids);
   free(kp);
   free(threads);
   free(attr);
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
   fprintf(stderr, "   -u[mnk] <#> : M/N/K loop unrolling is #\n");
   fprintf(stderr, "   -r <#> : set the # of times to call KMM\n");
   fprintf(stderr, "   -R <mf>: set # reps to force <mf> MFLOPs\n");
   fprintf(stderr, "   -F <kb> : set flush size in kilobytes\n");
   fprintf(stderr, "   -M[a,b,c] <#> : mov[A,B,C] = #\n");
   fprintf(stderr, "   -F[a,b,c] <kb> : set flush size for [A,B,C]\n");
   exit(iarg ? iarg : -1);
}

struct kmm_struct *GetFlags(int nargs, char **args, FILE **fpout)
{
   FILE *fp;
   struct kmm_struct *kp;
   double mflops=750.0;
   int i, j, IGMF=0;  /* ignore mflops.frc file? */

   *fpout = NULL;
   kp = malloc(sizeof(struct kmm_struct));
   assert(kp);
   kp->pids = NULL;
   kp->p = 1;
   kp->mb = kp->nb = kp->kb = 40;
   kp->mu = kp->nu = 4;
   kp->ku = 1;
   kp->movA = kp->movB = kp->movC = 0;
   kp->FLSIZE = L2SIZE;
   kp->flA = kp->flB = kp->flC = -1;
   kp->reps = 0;
   kp->vlen = 1;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
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
      case 'F':
         if (++i >= nargs)
            PrintUsage(args[0], i, "out of arguments");
         switch(args[i-1][2])
         {
         case 'c':
         case 'C':
            kp->flC = atoi(args[i]) * 1024;
            break;
         case 'b':
         case 'B':
            kp->flB = atoi(args[i]) * 1024;
            break;
         case 'a':
         case 'A':
            kp->flA = atoi(args[i]) * 1024;
            break;
         case 0:
            kp->FLSIZE = atoi(args[i]) * 1024;
            break;
         default:
            PrintUsage(args[0], i-1, "unknown flush operand");
         }
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
   #ifndef TIME_COPY
/*
 * If there is a tuned <upr>mflops.frc file, use this in preference to the
 * commandline argument unless the -Rf override flag was used
 */
   if (!IGMF)
   {
      char fn[16] = {'s', 'm', 'f', 'o', 'p', 's', '.', 'f', 'r', 'c', '\0'};
      #if defined(DREAL) || defined(DCPLX)
         fn[0] = 'd';
      #endif
      fp = fopen(fn, "r");
      if (fp)
      {
         assert(fscanf(fp, "%le", &mflops) == 1);
         fclose(fp);
      }
   }
   #endif
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
   #ifdef TIME_SYRKK
      kp->movB = 0;
   #endif
   if (kp->FLSIZE <= 0)
      kp->FLSIZE = L2SIZE;
   if (kp->flA > 0 || kp->flB > 0 || kp->flC > 0)
   {
      if (kp->flA < 0) kp->flA = kp->FLSIZE;
      if (kp->flB < 0) kp->flB = kp->FLSIZE;
      if (kp->flC < 0) kp->flC = kp->FLSIZE;
   }
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
   dp = TimeOnCores(kp);
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
