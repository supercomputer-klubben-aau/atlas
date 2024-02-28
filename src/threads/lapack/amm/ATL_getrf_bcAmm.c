/*
 * This file uses several macros to run with different configuration.
 *
 * #define DISABLE_LOOKAHEAD - For disabling infinite lookahead.
 * #define SERIAL_TRSM - For running USWP and TRSM in serial.
 * #define USE_SGETF2 -  For running getrf in serial for getf2.
 * #define NO_RECURSIVE_GETF2 - For iterative blocked getf2.
 * #define NO_GETF2 - For debugging, to see the impact of getf2.
 * #define ATL_LAPROF - For timings used my model driver.
 * #define SHOW_PROF - Used with ATL_LAPROF, to show the timings.
 * #define ATL_COLAFF - Not verified. To use column-major affinity.
 * #define DEBUG - Just adds some assert that helps debug race conditions.
 */
#define ATL_CBC_NAME /* to use the tgetrf_CBC name to be included in ATLAS. */
#define USE_BCSWP    /* to use ATLAS's block-cyclic amm swap routines. */


/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"


#include "atlas_cbc2d.h"

#include "atlas_threads.h"
#include "atlas_tlapack.h"
#include "atlas_pthreads.h"
#include "atlas_lvl3.h"
#include "atlas_lapack.h"
#include "atlas_lvl2.h"
#include "atlas_lamch.h"

#define UAMM_VER u0

#define UAMM_PRE Mjoin(Mjoin(atlas_,PRE),UAMM_VER)
#define UAMM_UPR Mjoin(Mjoin(atlas_,UPR),UAMM_VER)

#include Mstr(Mjoin(UAMM_PRE,amm_blk.h))
#include Mstr(Mjoin(UAMM_PRE,amm_flag.h))
#include Mstr(Mjoin(UAMM_PRE,amm_kern.h))
#include Mstr(Mjoin(UAMM_PRE,amm_perf.h))

#include Mstr(Mjoin(UAMM_PRE,amm_ablk2cmat.h))
#include Mstr(Mjoin(UAMM_PRE,amm_cmat2ablk.h))
#include Mstr(Mjoin(UAMM_PRE,amm_cm2am_an.h))
#include Mstr(Mjoin(UAMM_PRE,amm_cm2am_a1.h))
#include Mstr(Mjoin(UAMM_PRE,amm_am2cm_an.h))
#include Mstr(Mjoin(UAMM_PRE,amm_am2cm_a1.h))

enum error_codes {NOT_ENOUGH_MEMORY=-1, CASE_NOT_HANDLED=-2};

#define Mmod(x, y) ( (x) - ( ((x) / (y)) * (y) ) )

#define NUMROC(mn, id, MN, T) \
{ \
   int mn_p, mn_r; \
   mn = (MN) / (T); \
   mn_r = (MN) - (mn * (T)); \
   if ((id) < mn_r) mn += 1; \
}

#define RID(ID, PT, QT) ((ID) / (QT))
#define CID(ID, PT, QT) Mmod((ID), (QT))

#define GID(RID, CID, PT, QT) ((RID)*(QT) + (CID))

#define IDXL2G(IL, ID, NB_, NT) \
   (((IL)/(NB_))*(NB_)*(NT) + (ID)*(NB_) + Mmod((IL), (NB_)))

#define IDXG2L(IG, NB_, NT) ((NB_)*((IG)/((NB_)*(NT))) + Mmod((IG), (NB_)))
#define IDXG2P(IG, NB_, NT) Mmod(((IG)/(NB_)), (NT))

#define IDXL2BL(IL, NB_) (((IL)/(NB_))*(NB_)*(NB_) + Mmod((IL), (NB_)))

#define BLKG2L(BL, NT) ((BL)/(NT))
#define BLKL2G(BL, ID, NT) ((BL)*(NT) + (ID))

#define BLKG2P(BL, NT) Mmod((BL), (NT))

#define TRSM_BLK(KB, RID, PT) (((KB) - (RID)) / (PT))

#define NEXT_GEMM_BLK(KB, RID, PT) \
   (((KB)/(PT)) + ((Mmod((KB), (PT)) >= RID) ? 1 : 0))

#if 0
#else
   #define my_trsm Mjoin(PATL, trsm)
   #define my_swap Mjoin(PATL, swap)
#endif

#define myRBS (nbnb SHIFT) /* stride of real blks */
#define myIBS (myRBS)      /* stride of imag blks */

#ifdef TCPLX
   #if !defined(CAMM_USEA) && !defined(CAMM_USEB) && !defined(CAMM_USEC)
      #define CAMM_USEA /* provides better performance than USEC on opt32 */
   #endif
   #define ca2cc(m, n, alp, W, beta, A, lda) \
      ca2cc((m), (n), (alp), (W), (W)+nbnb, (beta), (A), (lda))
   #define cc2ca(m, n, alp, A, lda, beta, W) \
      cc2ca((m), (n), (alp), (A), (lda), (beta), (W), (W)+nbnb)

   #ifdef REALbIMAG    /* if defined, stores real before imag, no big effect */
      #define myROFF 0
      #define myIOFF nbnb
      #define c2a1(m, n, alp, A, lda, W) \
         c2a1((m), (n), (alp), (A), (lda), (W), (W)+myIOFF)
      #define a2c1(m, n, alp, A, lda, W) \
         a2c1((m), (n), (alp), (A), (lda), (W), (W)+myIOFF)

      #ifdef CAMM_USEA
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            rB = B0; \
            iB = B0 + nbnb; \
            rA = A0; \
            iA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #elif defined (CAMM_USEB)
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            rB = B0; \
            iB = B0 + nbnb; \
            rA = A0; \
            iA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #elif defined(CAMM_USEC)
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            rB = B0; \
            iB = B0 + nbnb; \
            rA = A0; \
            iA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #endif
   #else /* using imag, then real storage */
      #define myROFF nbnb
      #define myIOFF (-nbnb)
      #define c2a1(m, n, alp, A, lda, W) \
         c2a1((m), (n), (alp), (A), (lda), (W)+myROFF, (W))
      #define a2c1(m, n, alp, A, lda, W) \
         a2c1((m), (n), (alp), (A), (lda), (W)+myROFF, (W))

      #ifdef CAMM_USEA
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            iB = B0; \
            rB = B0 + nbnb; \
            iA = A0; \
            rA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #elif defined (CAMM_USEB)
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            iB = B0; \
            rB = B0 + nbnb; \
            iA = A0; \
            rA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #elif defined(CAMM_USEC)
         #define camm(nmu, nnu, kb, A0, B0, C0, amm_b1, amm_bn, nbnb, pfd) \
         { \
            TYPE *rA, *iA, *rB, *iB, *rC, *iC; \
            iB = B0; \
            rB = B0 + nbnb; \
            iA = A0; \
            rA = A0 + nbnb; \
            rC = C0; \
            iC = C0 + nbnb; \
            amm_bn(nmu, nnu, nb, iA, iB, rC, rA, rB, rC); \
            amm_bn(nmu, nnu, nb, rA, rB, rC, rA, iB, iC); \
            amm_b1(nmu, nnu, nb, rA, iB, iC, iA, rB, iC); \
            amm_b1(nmu, nnu, nb, iA, rB, iC, iA+((pfd) SHIFT), \
                   iB, rC+((pfd) SHIFT)); \
         }
      #endif
   #endif
   #define bcAm2rm(N, M, rW, rBS, A, lda, nb, nt, a2r) \
      Mjoin(PATL, bcAm2rm)((N), (M), (rW), (rBS), (rW)+(myIOFF), myIBS, \
                           (A), (lda), (nb), (nt), (a2r))
   #define bcRm2am(N, M, A, lda, rW, rBS, nb, nt, r2a) \
      Mjoin(PATL, bcRm2am)((N), (M), (A), (lda), (rW), (rBS), \
                           (rW)+(myIOFF), myIBS, (nb), (nt), (r2a))
   #define bcAblk2cmat(M, N, rW, rBS, A, lda, nb, nt, ca2cc) \
      Mjoin(PATL, bcAblk2cmat)((M), (N), (rW), (rBS), (rW)+nbnb, myIBS, \
                               (A), (lda), (nb), (nt), (ca2cc))
#else
   #define myROFF 0
   #define bcAm2rm Mjoin(PATL, bcAm2rm)
   #define bcRm2am Mjoin(PATL, bcRm2am)
   #define bcAblk2cmat Mjoin(PATL, bcAblk2cmat)
#endif

typedef struct
{
   int id, rid, cid;
   int pt, qt, nt;

   int nb, nbi;
   int mu, nu;
   int NMu, NNu, RMu, RNu;
   int rmu;

   cm2am_t r2a;
   cm2am_t c2a;
   am2cm_t a2r;
   am2cm_t a2c;
   ablk2cmat_t ca2cc;
   cmat2ablk_t cc2ca;
   ammkern_t amm_b1;
   ammkern_t amm_bn; /* needed for complex */

   const int major;
   TYPE *A;
   int lda;
   int *ipiv;
   int M, N, MN;
   int Rm, Rn;
   int oM, oN, oMb, oNb;
   int m, n;

   TYPE **WRKS;
   int *ldws;
   int *m_bs;
   int n_b;

   int *PROG;
   int *JB;
   ATL_CBC2D_t *cbc;
   int info;

   ATL_TGETF2_M_t *lups;
   ATL_INT *maxindx, *stage;
   TYPE **wrks2;
   ATL_INT *ldws2;
   ATL_INT ldw2;
   int info2;
   #ifdef USE_BCSWP
      ATL_bcpiv_t *bp;
   #endif

   char space[256 - (28)*sizeof(int) - 14*sizeof(void*)];
} ATL_GETRF_2D_t;

#ifdef ATL_LAPROF
   #define DEFINE_PROFVARS
   #include "atlas_laprof.h"
#endif

#ifdef SHOW_PROF
void PrintAverageProfiling(int nthreads, int ntimers, int tstart, int tend,
                           long long timers[ATL_NTHREADS][ntimers],
                           char** timerNames, int ignore)
{
   long long ltsc, ctsc, ttsc;
   double ttime, ctime;
   int k, t;
   fprintf(stdout, "GEMM MFLOPs : ");
   for (k=0; k<nthreads; k++)
   {
      if (timers[k][GEMM] <= 0 || timers[k][MMFLOP] <= 0) continue;
      fprintf(stdout, "%.4lf ",
            timers[k][MMFLOP]*1e-9/Tick2Sec(timers[k][GEMM]));
      timers[k][MMFLOP] = 0;
   }
   fprintf(stdout, "\n");
   ttsc = 0;
   for (t=tstart; t<tend; t++)
   {
      ctsc = 0;
      for (k=0; k<nthreads; k++)
      {
         ctsc += Timers[k][t];
      }
      ctsc /= nthreads;
      ttsc += ctsc;
   }
   ttime = Tick2Sec(ttsc)*1e3;
   fprintf(stdout, "Total time : %lf ms\n", ttime);
   fprintf(stdout, "--------------------------------------\n");
   for (t=tstart; t<tend; t++)
   {
      ctsc = 0;
      for (k=0; k<nthreads; k++)
      {
         ctsc += Timers[k][t];
      }
      ctsc /= nthreads;
      if (ctsc <= 0) continue;
      ctime = Tick2Sec(ctsc)*1e3;
      fprintf(stdout, "%15s : %lf ms (%.2lf%%)\n", timerNames[t],
                  ctime, ctime*100.0/ttime);
   }
}

void PrintProfilingInfo_pt(int nthreads, int ntimers,
                           long long timers[ATL_NTHREADS][ntimers],
                           char** timerNames, int ignore)
{
   long long LastTimeStamp, CurrentTimeStamp;
   double time, ctime;
   int k, t;
   LastTimeStamp = 0;
   for (k=0; k<nthreads; k++)
   {
      CurrentTimeStamp = 0;
      for (t=0; t<ntimers; t++)
      {
         CurrentTimeStamp += timers[k][t];
      }
      if (ignore >= 0)
         CurrentTimeStamp -= timers[k][ignore];
      ctime = Tick2Sec(CurrentTimeStamp)*1e3;
      fprintf(stdout, "\nThread %2d : %lf ms\n", k, ctime);
      fprintf(stdout, "--------------------------------------\n");
      for (t=0; t<ntimers; t++)
      {
         if (timers[k][t] <= 0) continue;
         time = Tick2Sec(timers[k][t])*1e3;
         fprintf(stdout, "%15s : %lf ms (%.2lf%%)\n", timerNames[t],
                  time, time*100.0/ctime);
      }
      if (CurrentTimeStamp > LastTimeStamp)
         LastTimeStamp = CurrentTimeStamp;
   }
   fprintf(stdout, "\nMax total time : %e secs\n", Tick2Sec(LastTimeStamp));
}
#endif

/*
 * RCW: added IAMM param
 */
static void SelectGenKernels(int IDXAMM, int *nb, int *mu, int *nu,
                        cm2am_t *r2a_an, cm2am_t *c2a_a1,
                        am2cm_t *a2r_an, am2cm_t *a2c_a1,
                        ablk2cmat_t *ca2cc, cmat2ablk_t *cc2ca,
                        ammkern_t *ammkrnl_b1, ammkern_t *ammkrnl_bn
                        )
{
   int i, index = -1, i2, ik1=-1;
   if (IDXAMM < 0 || IDXAMM > ATL_UAMM_NCASES)
   {
      for (i=0; i<ATL_UAMM_NCASES; i++)
      {
         if (ATL_UAMM_KBs[i] == *nb)
         {
            index = i;
            break;
         }
         else if (*nb < ATL_UAMM_KBs[i] && i>0 && *nb > ATL_UAMM_KBs[i-1])
         {
            index = ik1;
            break;
         }
         if (ATL_UAMM_KFLAG[i]==1)
            ik1 = i; /* store the largest KB kernel with KU=1 */
      }
   }
   else
      index = IDXAMM;
   /* if (*nb > ATL_UAMM_KBs[ATL_UAMM_NCASES-1]) index = ATL_UAMM_NCASES-1; */
   ATL_assert(index >= 0);
   #ifdef ATL_LAPROF
      printf("Choosing KB:%d Kernel with MU=%d NU=%d KU=%d\n",
         ATL_UAMM_KBs[index], ATL_UAMM_MUs[index], ATL_UAMM_NUs[index],
         ATL_UAMM_KUs[index]);
   #endif
   *nb = ATL_UAMM_KBs[index];
   *mu = ATL_UAMM_MUs[index];
   *nu = ATL_UAMM_NUs[index];
   *r2a_an = ATL_UAMM_AT2BLK_an[index];
   *c2a_a1 = ATL_UAMM_B2BLK_a1[index];
   *a2r_an = ATL_UAMM_BLK2AT_an[index];
   *a2c_a1 = ATL_UAMM_BLK2B_a1[index];
   *ca2cc = ATL_UAMM_BLK2C_a1_b0[index];
   *cc2ca = ATL_UAMM_C2BLK_a1_b0[index];
   *ammkrnl_b1 = ATL_UAMM_KERN_b1[index];
   *ammkrnl_bn = ATL_UAMM_KERN_bn[index];
}

void SelectBlockAndGrid(ATL_CINT M, ATL_CINT N, int *nb, int *nbi,
                        int *pt, int *qt)
{

   #if defined(ATL_LAPROF) && !defined(SHOW_PROF)
      Mjoin(PATL, ilaenv_CBC)(M, N, nb, pt, qt);
      return;
   #endif
   #ifdef AUTO_DIST
      const int L3_per_chip = 5<<20; /* 5MB for opt32 */
      const int tnt = ATL_NTHREADS;
      int ipt = 4;
      int iqt = 8;
      int ikb = 12; /* ATL_UAMM_66KB */
      int gkb = 72; /* ATL_UAMM_98KB */
      int gkb_idx = 11; /* ATL_UAMM_66IDX */
      int mkb = 168; /* ATL_AMM_MAXKB */
      int mkb_idx = 19; /* IDX_OF(ATL_AMM_MAXKB) */
      int inb = Mmin(Mmin(M/ipt, N/iqt), ikb);
      int fnb=0, fpt=0, fqt=0;
      int small_range_limit = 1000; /* need to change for precision */
      int size_ratio = M / small_range_limit;
      int min_panels_per_thread = 4;
      int tnb;
      if (size_ratio == 0)
      {
         fnb = inb;
         fpt = 1;
         fqt =  Mmin(N / (fnb*min_panels_per_thread), tnt);
      }
      else
      {
         tnb = ATL_UAMM_KBs[gkb_idx +
            (int)((mkb_idx-gkb_idx)*0.5*Mmin(1.0, size_ratio/40.0))];
         fnb = Mmin(size_ratio*inb, tnb);
         fpt = ipt;
         fqt =  Mmin(N / (fnb*min_panels_per_thread), tnt/ipt);
      }
      *nb = fnb;
      *pt = fpt;
      *qt = fqt;
      if (*nb > 12) *nbi = 12;
      else *nbi = 4;
   #else
      #ifndef PT
         #define PT 2
      #endif
      #ifndef QT
         #define QT 2
      #endif

      *pt = PT;
      *qt = QT;
      #ifdef NB_d
         *nb = NB_d;
      #else
         *nb = ATL_UAMM_80KB;
      #endif
      #ifdef NB_i
         *nbi = NB_i;
      #else
         if (*nb > 12) *nbi = 12;
         else *nbi = 4;
      #endif
   #endif
}


/*
 * Swap one row that is in col-major with one that is in C-major (access).
 *
 * n: number of elements to swap.
 * start: start index of the row.
 * an: total number of elements of that row.
 * A: points to the row that needs to be swapped.
 * W: points to the start of the C-major block.
 * pi: index of the pivot row in W.
 */
static void ATL_amswap(int n, int start, int an, int mu, int nu,
                        TYPE *A, int nb, TYPE *W, int pi)
{
   #define FindDest(i, j, np, mu) \
      (((i)/(mu))*(mu)*(np) +(j)*(mu) +Mmod((i),(mu)))

   int i, j;
   TYPE t, *Ac = A + nb*(start SHIFT);
   #ifdef TCPLX
      TYPE *iW = W + nb*nb;
   #endif

   /* find the destination of the first element */
   j = FindDest(pi, start, an, mu);
   for (i=0; i<n; i++, Ac+=(nb SHIFT), j+=mu)
   {
      /* swap Ac[0] with W[j] */
      t = Ac[0];
      Ac[0] = W[j];
      W[j] = t;
      #ifdef TCPLX
         t = Ac[1];
         Ac[1] = iW[j];
         iW[j] = t;
      #endif
   }
}

/*
 * Most complicated routine to apply pivoting on a block-major block
 * with other access-major blocks.
 *
 * W  : points to the start of the block.
 * jb : local panel no.
 * kb : block no. for which pivoting will be applied.
 * start/end : needed for parallel swap.
 */
static void ATL_amlaswp(int id, int start, int an, int n,
                     int jb, int kb, TYPE *W, ATL_GETRF_2D_t *pd)
{
   const int pt = pd->pt, qt = pd->qt;
   const int cid = CID(id, pt, qt);
   const int nb = pd->nb;
   const int mu = pd->mu;
   const int nu = pd->nu;
   TYPE **WRKS = pd->WRKS;
   int *ldws = pd->ldws;
   int *ipiv = pd->ipiv;
   int i, I, K = kb*nb;
   int partner, cbb, irem;
   int nbnb2 = nb*nb*2;
   TYPE *A0;
   ipiv += K;
   cbb = K+nb; /* current block bound */
   for (i=0; i<nb; i++)
   {
      I = ipiv[i]; /* I is the global index of pivot row */
      partner = GID(BLKG2P(I/nb, pt), cid, pt, qt);
      if (partner != id || I >= cbb) /* swap between block and access major */
      {
         I = IDXG2L(I, nb, pt); /* I is now the local index for partner */
         irem = Mmod(I, nb); /* irem is index within the block */
         I = IDXL2BL(I-irem, nb); /* correct the offset for block-major style */

         A0 = WRKS[partner]+ ((jb*nb*((size_t)(ldws[partner])) /* correct pnl */
            + I + nbnb2) SHIFT); /* offset of two blocks for C major blocks */

         ATL_amswap(n, start, an, mu, nu, W+(i SHIFT), nb, A0, irem);
      }
      else /* swap within current block */
      {
         I = IDXG2L(I-K, nb, pt); /* I is now the local index ??? */
         if (i != I) /* we need to swap */
            my_swap(n, W+((start*nb+i)SHIFT), nb, W+((start*nb+I) SHIFT), nb);
      }
   }
}

/*
 * Applies Swap and TRSM to one block. It has to find its data pointers.
 */
static void ApplySwapTRSM(int cjb, int kb, int nnu, int cn, ATL_GETRF_2D_t *pd)
{
   const int id = pd->id;
   const int rid = pd->rid, cid = pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   const int mu = pd->mu, nu = pd->nu;
   const int nbnb = nb * nb;
   const int active = Mmod(kb, pt);
   TYPE **WRKS = pd->WRKS;
   int *ldws = pd->ldws;
   TYPE *A0, *B0, *Bu;
   cm2am_t c2a1 = pd->c2a;
   #ifdef TREAL
      const TYPE one = ATL_rone, zero = ATL_rzero;
   #elif defined(TCPLX)
      const TYPE one[2]={ATL_rone,ATL_rzero}, zero[2]={ATL_rzero,ATL_rzero};
   #endif

   #ifdef SERIAL_TRSM
      if (rid == active)
      {
         int partner = GID(active, Mmod(kb,qt), pt, qt);
         A0 = WRKS[partner] + (((kb/qt)*nb*((size_t)(ldws[partner]))) SHIFT);
         A0 += (((kb/pt)*nbnb) SHIFT);
         B0 = WRKS[id] + ((cjb*nb*((size_t)(ldws[id]))) SHIFT);
         B0 += ((TRSM_BLK(kb, active, pt)*nbnb + nbnb) SHIFT);
         #ifdef USE_BCSWP
            Mjoin(PATL,bcLaswp_amm)(pd->bp, nnu*nu, 0, B0, nb, kb*nb, kb*nb+nb,
                                    cjb, nnu, nbnb);
         #else
            ATL_amlaswp(id, 0, nnu*nu, nnu*nu, cjb, kb, B0, pd);
         #endif
         my_trsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                  nb, cn, one, A0, nb, B0, nb);

         c2a1(nb, cn, one, B0, nb, B0-(nbnb SHIFT));
      }
   #else /* Parallel TRSM */
      int ts = rid*nnu / pt;
      int tn = ((rid+1)*nnu/pt) - ts;
      int partner;
      ts *= nu;
      tn *= nu;
      /* find the A block for TRSM */
      partner = GID(active, Mmod(kb,qt), pt, qt);
      A0 = WRKS[partner] + (((kb/qt)*nb*((size_t)(ldws[partner]))) SHIFT);
      A0 += (((kb/pt)*nbnb) SHIFT);
      /* find the block for SWAP and TRSM */
      partner = GID(active, cid, pt, qt);
      B0 = WRKS[partner] + ((cjb*nb*((size_t)(ldws[partner]))) SHIFT);
      B0 += ((TRSM_BLK(kb, active, pt)*nbnb + nbnb) SHIFT);

      #ifdef USE_BCSWP
         Mjoin(PATL,bcLaswp_amm)(pd->bp, tn, ts, B0, nb, kb*nb, kb*nb+nb,
                                 cjb, nnu, nbnb);
      #else
         ATL_amlaswp(partner, ts, nnu*nu, tn, cjb, kb, B0, pd);
      #endif

      Bu = B0 + ts*nb; /* needed for complex B-major block */
      B0 = B0 + ((ts*nb) SHIFT);
      my_trsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
               nb, tn, one, A0, nb, B0, nb);

      c2a1(nb, tn, one, B0, nb, Bu-(nbnb SHIFT));
   #endif
}

/*
 * Applies GEMM, if provided the correct pointers.
 * sb: starting block number for GEMM.
 */
static void ApplyGEMM(int kb, int sb, int m_b, int nmu, int rmu, int nnu,
               TYPE* A0, TYPE* B0, TYPE* C0, const int nb,
               ammkern_t amm_b1, ammkern_t amm_bn)
{
   #ifndef FWD_GEMM
      #define FWD_GEMM 0
   #endif
   const nbnb = nb*nb;
   int ib, mb = m_b - sb;
   if (FWD_GEMM || (kb&1)) /* do top-down gemm */
   {
      if (nmu)
      {
         for (ib=0; ib<mb; ib++, A0+=(nbnb SHIFT), C0+=(nbnb SHIFT))
         {
            #ifdef TCPLX
               camm(nmu, nnu, nb, A0, B0, C0, amm_b1, amm_bn, nbnb, nbnb);
            #else
               amm_b1(nmu, nnu, nb, A0, B0, C0, A0+nbnb, B0, C0+nbnb);
            #endif
         }
      }
      if (rmu)
      {
         #ifdef TCPLX
            camm(rmu, nnu, nb, A0, B0, C0, amm_b1, amm_bn, nbnb, nbnb);
         #else
            amm_b1(rmu, nnu, nb, A0, B0, C0, A0+nbnb, B0, C0+nbnb);
         #endif
      }
   }
   else /* do bottom-up gemm */
   {
      A0 += mb*(nbnb SHIFT);
      C0 += mb*(nbnb SHIFT);
      if (rmu)
      {
         #ifdef TCPLX
            camm(rmu, nnu, nb, A0, B0, C0, amm_b1, amm_bn, nbnb, -nbnb);
         #else
            amm_b1(rmu, nnu, nb, A0, B0, C0, A0-nbnb, B0, C0-nbnb);
         #endif
      }
      A0-=(nbnb SHIFT);
      C0-=(nbnb SHIFT);
      if (nmu)
      {
         for (ib=mb; ib; ib--, A0-=(nbnb SHIFT), C0-=(nbnb SHIFT))
         {
            #ifdef TCPLX
               camm(nmu, nnu, nb, A0, B0, C0, amm_b1, amm_bn, nbnb, -nbnb);
            #else
               amm_b1(nmu, nnu, nb, A0, B0, C0, A0-nbnb, B0, C0-nbnb);
            #endif
         }
      }
   }
}

/*
 * This function will apply one phase of update on a panel.
 * jb: local panel no.
 * Jb: global panel no.
 * kb: which panel the update is being applied from.
 * Wj: points to first TRSM block
 */
static void ApplyOneUpdate(int cjb, int cJb, int kb, int nnu, int cn,
                     TYPE *Wj, TYPE *Aj, ATL_GETRF_2D_t *pd)
{
   const int id = pd->id;
   const int rid = pd->rid, cid = pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   const int nu = pd->nu;
   const int nbnb = nb * nb;
   ATL_CBC2D_t *cbc = pd->cbc;
   TYPE **WRKS = pd->WRKS;
   int *ldws = pd->ldws;
   int sb, partner;
   int active = Mmod(kb, pt);
   TYPE *A0, *B0, *C0, *Wk;
   #ifdef TREAL
      TYPE one = ATL_rone, zero = ATL_rzero;
   #elif defined(TCPLX)
      const TYPE one[2]={ATL_rone,ATL_rzero}, zero[2]={ATL_rzero,ATL_rzero};
   #endif
   #ifdef ATL_LAPROF
      long long ltsc=0, ctsc=0, ttsc=0;
      rdtsc(ltsc);
   #endif

   ApplySwapTRSM(cjb, kb, nnu, cn, pd);

   /* sync to confirm that everyone is done with TRSM */
   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][SWSM], ttsc);
   #endif

   /* find first GEMM block */
   sb = NEXT_GEMM_BLK(kb, rid, pt);
   C0 = Wj + ((sb*nbnb + nbnb) SHIFT);
   /* find first A block for GEMM */
   partner = GID(rid, Mmod(kb,qt), pt, qt);
   A0 = WRKS[partner] + (((kb/qt)*nb*((size_t)(ldws[partner]))) SHIFT);
   A0 += ((sb*nbnb) SHIFT);
   /* find first B block for GEMM */
   partner = GID(active, cid, pt, qt);
   B0 = WRKS[partner] + ((cjb*nb*((size_t)(ldws[partner]))) SHIFT);
   B0 += ((TRSM_BLK(kb, active, pt)*nbnb) SHIFT);

   ApplyGEMM(kb, sb, pd->m_bs[id], pd->NMu, pd->RMu, nnu,
               A0, B0, C0, nb, pd->amm_b1, pd->amm_bn);
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][GEMM], ttsc);
   #endif

   if (rid == active)
   {
      #ifndef EARLY_B_COPY
         if ((kb+1) < cJb) /* if i am not active  do the copy, else postpone */
      #endif
         {
            /* copy back U block to original storage */
            A0 = Aj + ((kb*nb) SHIFT);
            B0 = B0 + (nbnb SHIFT);
            Mjoin(PATL, gecopy)(nb, cn, B0, nb, A0, pd->lda);
         }
   }
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][Copy2GU], ttsc);
   #endif
   kb++;
   if (++active == pt) active = 0;
   if (rid == active && kb<cJb) /* more updates needed */
   {
      /* shift top block for swap and trsm */
      pd->ca2cc(nb, nnu*nu, one, C0, zero, C0-(nbnb SHIFT), nb);
   }
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][CopyS], ttsc);
   #endif
   /* sync before next phase of update */
   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][CSync2M], ttsc);
   #endif
   if (rid == active)
      pd->PROG[cJb] = kb; /* kb already incremented */
   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][CSync2M], ttsc);
   #endif
}

static void ApplyLeftPivots(ATL_GETRF_2D_t *pd)
{
   const int id = pd->id, rid=pd->rid, cid=pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   int *JB = pd->JB;
   const int MN = pd->MN;
   TYPE *A = pd->A;
   const int lda = pd->lda;
   int *ipiv = pd->ipiv;
   volatile int *PROG = (volatile int *)pd->PROG;
   const int oNb = pd->oNb;
   const int aNb = oNb-1;
   const int nbnt = nb * nt;
   const int nbpt = nb * pt;
   const int nbqt = nb * qt;
   const size_t nbntlda = nbnt * lda;
   const size_t nbptlda = nbpt * lda;
   const size_t nbqtlda = nbqt * lda;
   int cJB, mn;
   const int ts = rid*nb/pt;
   const int tn = (rid+1)*nb/pt - ts;
   const int tslda = ts * lda;
   int j, Jb;
   TYPE *A0;
   #ifdef ATL_LAPROF
      long long ltsc=0, ctsc=0, ttsc=0;
      rdtsc(ltsc);
   #endif
   #if 0
      cJB = *JB;
      mn = cJB * nb;
      mn = Mmin(mn, MN);
      j = cid*nb;
      if (j+nb > mn) mn = j+nb;
      A0 = A + ((j*((size_t)(lda))) SHIFT);
      A0 += (tslda SHIFT);
      for (Jb=cid; Jb<aNb; Jb+=qt, j+=nbqt, A0+=(nbqtlda SHIFT))
      {
         if (PROG[Jb]*nb < mn)
            ATL_laswp(tn, A0, lda, j+nb, mn, ipiv, 1);
      }
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][LSWP], ttsc);
      #endif
      if (mn < MN)
      {
         while (PROG[aNb] < oNb);
         #ifdef ATL_LAPROF
            gettime(ltsc, ctsc, Timers[id][CSync2S], ttsc);
         #endif
         j = cid*nb;
         A0 = A + ((j*((size_t)(lda))) SHIFT);
         A0 += (tslda SHIFT);
         for (Jb=cid; Jb<aNb; Jb+=qt, A0+=(nbqtlda SHIFT))
         {
            ATL_laswp(tn, A0, lda, mn, MN, ipiv, 1);
         }
      }
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][LSWP], ttsc);
      #endif
   #else /* new kind of lswp, doesn't work well with zero drain */
      /* first pass by first last pcol threads */
      int flpcol = Mmod(oNb, qt);
      cJB = oNb - qt + 1;
      mn = cJB * nb;
      mn = Mmin(mn, MN);
      if (flpcol == cid)
      {
         j = rid*nb;
         if (j+nb > mn) mn = j+nb;
         A0 = A + ((j*((size_t)(lda))) SHIFT);
         for (Jb=rid; Jb<cJB; Jb+=pt, j+=nbpt, A0+=(nbptlda SHIFT))
         {
            if (j+nb < mn) ATL_laswp(nb, A0, lda, j+nb, mn, ipiv, 1);
         }
      }
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][LSWP], ttsc);
      #endif
      /* wait for last panel */
      ATL_CBC2D_barrier(ATL_SYNC_GRID, pd->cbc, id);
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][CSync2S], ttsc);
      #endif
      /* second pass by all threads */
      j = id * nb;
      A0 = A + ((j*((size_t)(lda))) SHIFT);
      for (Jb=id; Jb<aNb; Jb+=nt, j+=nbnt, A0+=(nbntlda SHIFT))
      {
         if (Mmax(mn, j+nb) < MN)
            ATL_laswp(nb, A0, lda, Mmax(mn, j+nb), MN, ipiv, 1);
      }
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][LSWP], ttsc);
      #endif
   #endif
}

/*
 * Copy one full panel to local workspace.
 * m: rounded up
 * n: actual n
 * Aj: pointer to original storage for gpanel Jb
 * Wj: pointer to local workspace for gpanel Jb
 */
static void CopyInPanel(int Jb, int m, int rm, int n,
                  TYPE *Aj, const int lda, TYPE *Wj, ATL_GETRF_2D_t *pd)
{
   const int id = pd->id;
   const int rid = pd->rid, cid = pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   cmat2ablk_t cc2ca = pd->cc2ca;
   const int nbnb = nb*nb;
   const int nbpt = nb*pt;
   const int rn = nb - n; /* rn cols to zero */
   const int rnnb = rn*nb;
   int i, tm;
   TYPE *Ac, *Wc;
   #ifdef TREAL
      TYPE one = ATL_rone, zero = ATL_rzero;
   #elif defined(TCPLX)
      const TYPE one[2]={ATL_rone,ATL_rzero}, zero[2]={ATL_rzero,ATL_rzero};
   #endif

   if(!m) return; /* if nothing to copy, return immediately */

   if (rm) m -= nb;
   if (Jb)  /* later panels - copy to C-major (access) */
   {
      Wc = Wj;
      Ac = Aj + ((rid*nb) SHIFT);
      tm = nb;
      if (m==0) /* we have only one block */
      {
         tm = rm;
         rm = 0;
      }
      #ifndef FWD_COPY
         /* copy first block differently */
         if (rid > 0)
         {
            Wc += (nbnb SHIFT);
            if (tm<nb || n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(tm, n, one, Ac, lda, zero, Wc);
         }
         else
         {
            if (tm<nb || n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            Mjoin(PATL, gecopy)(tm, n, Ac, lda, Wc, nb);
            Wc += (nbnb SHIFT);
         }
         Wc += (nbnb SHIFT);
         Ac += (nbpt SHIFT);
         for (i=nb; i<m; i+=nb, Ac+=(nbpt SHIFT), Wc+=(nbnb SHIFT))
         {
            if (n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(nb, n, one, Ac, lda, zero, Wc);
         }
         if (rm)
         {
            Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(rm, n, one, Ac, lda, zero, Wc);
         }
      #else
         Ac += (((m/nb)*nbpt) SHIFT);
         Wc += (((m/nb)*nbnb) SHIFT);
         if (rm)
         {
            Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(rm, n, one, Ac, lda, zero, Wc);
         }
         Ac -= (nbpt SHIFT);
         Wc -= (nbnb SHIFT);
         for (i=m; i!=nb; i-=nb, Ac-=nbpt, Wc-=nbnb)
         {
            if (n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(nb, n, one, Ac, lda, zero, Wc);
         }
         /* now copy the first block */
         if (rid > 0)
         {
            if (tm<nb || n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            cc2ca(tm, n, one, Ac, lda, zero, Wc);
            Wc -= (nbnb SHIFT);
         }
         else
         {
            Wc -= (nbnb SHIFT);
            if (tm<nb || n<nb) Mjoin(PATL, gezero)(nb, nb, Wc, nb);
            Mjoin(PATL, gecopy)(tm, n, Ac, lda, Wc, nb);
         }
      #endif
   }
   #if !defined(USE_SGETF2)
   else     /* first panel - copy to block major */
   {
      Wc = Wj;
      Ac = Aj + ((rid*nb) SHIFT);
      #ifndef FWD_COPY
         for (i=0; i<m; i+=nb, Ac+=(nbpt SHIFT), Wc+=(nbnb SHIFT))
         {
            Mjoin(PATL, gecopy)(nb, n, Ac, lda, Wc, nb);
            Mjoin(PATL, gezero)(nb, rn, Wc+(rnnb SHIFT), nb);
         }
         if (rm)
         {
            Mjoin(PATL, gecopy)(rm, n, Ac, lda, Wc, nb);
            Mjoin(PATL, gezero)(nb-rm, n, Wc+(rm SHIFT), nb);
            Mjoin(PATL, gezero)(nb, rn, Wc+(rnnb SHIFT), nb);
         }
      #else
         Ac += (((m/nb)*nbpt) SHIFT);
         Wc += (((m/nb)*nbnb) SHIFT);
         if (rm)
         {
            Mjoin(PATL, gecopy)(rm, n, Ac, lda, Wc, nb);
            Mjoin(PATL, gezero)(nb-rm, n, Wc+(rm SHIFT), nb);
            Mjoin(PATL, gezero)(nb, rn, Wc+(rnnb SHIFT), nb);
         }
         Ac -= (nbpt SHIFT);
         Wc -= (nbnb SHIFT);
         for (i=m; i; i-=nb, Ac-=(nbpt SHIFT), Wc-=(nbnb SHIFT))
         {
            Mjoin(PATL, gecopy)(nb, n, Ac, lda, Wc, nb);
            Mjoin(PATL, gezero)(nb, rn, Wc+(rnnb SHIFT), nb);
         }
      #endif
   }
   #endif
}

static void ATL_getf2s(int id, int Jb, int m, int n, TYPE *W,
                        TYPE *Aj, int lda, int *ipiv, ATL_GETRF_2D_t *pd)
{
   const int rid = pd->rid, cid = pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   const int nbnb = nb * nb;
   const int active = Mmod(Jb,pt);
   ATL_GETRF_2D_t *pda, *pdb;
   int terr = 0;
   ATL_CBC2D_t *cbc = pd->cbc;
   TYPE *Ac = Aj + Jb*(nb SHIFT); /* Ac points to the diagonal block now */
   if (rid < active)
      Ac += (pt+rid-active)*(nb SHIFT);
   else
      Ac += (rid-active)*(nb SHIFT);
   if (Jb)
      bcAblk2cmat(m, n, W, myRBS, Ac, lda, nb, pt, pd->ca2cc);
   #if !defined(USE_SGETF2)
   else
      Mjoin(PATL, bcL2G_blkcpy)(m, n, W, nb, Ac, lda, pt);
   #endif

   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);

   if (rid == active)
   {
      terr = Mjoin(PATL, getrf)(CblasColMajor, pd->M-Jb*nb, n,
                                       Ac, lda, ipiv);
      if (terr && !pd->info2) pd->info2 = terr;
   }

   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);

   /* everyone read the error code from active guy */
   pdb = pd - id;
   pda = pdb + GID(active, cid, pt, qt);
   pd->info2 = pda->info2;

   if (Jb) W-=(nbnb SHIFT); /* for later panels, blocks are shifted */
   if (rid == active)
   {
      Mjoin(PATL, gecopy)(Mmin(m, nb), n, Ac, lda, W-(nbnb SHIFT), nb);
      bcRm2am(m-Mmin(m,nb), n, Ac+nb*(pt SHIFT), lda,
              W+myROFF, myRBS, nb, pt, pd->r2a);
   }
   else
      bcRm2am(m, n, Ac, lda, W-(nbnb SHIFT)+myROFF, myRBS, nb, pt, pd->r2a);

   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
}


/*
 * W: points to j00-th column and 0 or j00-th row.
 */
static void ATL_subtgetf2(int rank, int m00, int n, int j00, int p,
                               TYPE *W, int ldw, int *ipiv,
                               ATL_GETRF_2D_t *pd)
{
   const int nb = pd->nb;
   const int pt = pd->pt;
   TYPE *w = W, *v;
   int m0, mr, pivrank, m=m00;
   TYPE pivval, apv, apv2;
   volatile ATL_INT *maxindx = pd->maxindx, *stage=pd->stage;
   TYPE **WRKS = pd->wrks2;
   ATL_INT *ldws2 = pd->ldws2;
   ATL_INT locpiv, globpiv, i, j, k;
   #ifdef TCPLX
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #endif
   for (j=0; j < n; j++, w += (ldw SHIFT))
   {
      locpiv = cblas_iamax(m, w, 1);
      /*
       * Combine local pivot into global
       */
      if (!rank)
      {
         globpiv = IDXL2G(j+locpiv+j00, rank, nb, p); /* adjust for subpanel */
         pivrank = 0;
         apv = Mabs(w[locpiv SHIFT]);
         #ifdef TCPLX
            apv += Mabs(w[locpiv+locpiv+1]);
         #endif
         for (i=1; i < p; i++)
         {
            while(stage[i] < j);
            k = maxindx[i];
            apv2 = WRKS[i][((j+j00)*ldws2[i]+k)SHIFT];
            apv2 = Mabs(apv2);
            #ifdef TCPLX
               apv2 += Mabs(WRKS[i][(((j+j00)*ldws2[i]+k)<<1)+1]);
            #endif
            if (apv < apv2)
            {
               apv = apv2;
               globpiv = k;
               pivrank = i;
            }
            maxindx[i] = -1;
         }
         if (pivrank)
         {
            ipiv[j] = IDXL2G(globpiv, pivrank, nb, p);
            my_swap(n, W+(j SHIFT), ldw,
                       WRKS[pivrank]+((j00*ldws2[pivrank]+globpiv) SHIFT),
                       ldws2[pivrank]);
         }
         else
         {
            ipiv[j] = globpiv;
            if (globpiv != j)
               my_swap(n, W+(j SHIFT), ldw, W+((j+locpiv) SHIFT), ldw);
         }
         stage[0] = j;
         m--;                                           /* just finished */
         #ifdef TCPLX
            w += 2;                                     /* one row */
         #else
            w++;                                        /* one row */
         #endif
      }
      else /* all threads except 0 write their results, and await 0 */
      {
         maxindx[rank] = locpiv;
         stage[rank] = j;
         while (stage[0] < j);
      }
      #ifdef TCPLX
         v = &WRKS[0][((j+j00)*ldws2[0]+j+j00)SHIFT];
         if (*v != ATL_rzero || v[1] != ATL_rzero)
         {
            if (ATL_lapy2(v[0], v[1]) >= ATL_laSAFMIN)
            {
               TYPE inv[2];
               Mjoin(PATL,cplxinvert)(1, v, 1, inv, 1);
               cblas_scal(m, inv, w, 1);
            }
            else
            {
               Mjoin(PATL,cplxdivide)(m, v, w, 1, w, 1);
            }
         }
      #else
         pivval = WRKS[0][(j+j00)*ldws2[0]+j+j00];
         if (pivval != ATL_rzero)
         {
            if (Mabs(pivval) >= ATL_laSAFMIN)
               cblas_scal(m, ATL_rone/pivval, w, 1);
            else
            {
               for (k=0; k < m; k++)
                  w[k] /= pivval;
            }
         }
      #endif
      else /* pivot is zero, we have a singular matrix! */
         if (!pd->info2) pd->info2 = j+j00+1;   /* save if it is first */

      #ifdef TCPLX
         Mjoin(PATL,geru_L2)(m, n-j-1, none, w, 1,
                             WRKS[0]+(((j+j00)*(ldws2[0]+1)+ldws2[0])SHIFT),
                             ldws2[0], w+ldw+ldw, ldw);
      #else
         Mjoin(PATL,ger_L2)(m, n-j-1, ATL_rnone, w, 1,
                             WRKS[0]+(j+j00)*(ldws2[0]+1)+ldws2[0],
                             ldws2[0], w+ldw, ldw);
      #endif
   }
   stage[rank] = n;  /* let core 0 know we are done */
}

static void ATL_laswp_blk(int m, int n, int rs, int cs, TYPE *W, int ldw,
                          TYPE **WRKS, int *ldws, int *ipiv,
                          ATL_GETRF_2D_t *pd)
{
   const int nb = pd->nb;
   const int pt = pd->pt;
   int i, gpiv, lpiv, owner;
   TYPE *Ws, *Wc =  W+((cs*ldw)SHIFT);
   int r;
   if (n <= 0) return;
   for (i=0, r=rs; i<m; i++, r++)
   {
      gpiv = ipiv[i];
      lpiv = IDXG2L(gpiv, nb, pt);
      owner = IDXG2P(gpiv, nb, pt);
      if (owner == 0)
      {
         my_swap(n, Wc+(r SHIFT), ldw, Wc+(lpiv SHIFT), ldw);
      }
      else
      {
         Ws = WRKS[owner]+((cs*ldws[owner]+lpiv) SHIFT);
         my_swap(n, Wc+(r SHIFT), ldw, Ws, ldws[owner]);
      }
   }
}

/*
 * W: points to actual start of the block.
 * m: won't change during calls, need to call tgetf2 and gemm accordingly.
 */
static void ATL_rtgetf2(int rank, int Jb, int m, int n, int p, int j,
                        TYPE *W, int ldw, int *ipiv, ATL_GETRF_2D_t *pd)
{
   TYPE *A0, *B0, *C0;
   #ifdef TCPLX
      const TYPE one[2] = {ATL_rone, ATL_rzero};
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #else
      const TYPE one = ATL_rone;
      const TYPE none = ATL_rnone;
   #endif
   if (n > 4)
   {
      int nl, nr, nm, off;
      /* nl = (((n >> 1)+11)/12)*12; */
      nl = (((n >> 1)+3)/4)*4;
      nr = n - nl;
      ATL_rtgetf2(rank, Jb, m, nl, p, j, W, ldw, ipiv, pd);
      if (rank == 0)
      {
         /* right swap */
         ATL_laswp_blk(nl, nr,j, j+nl, W, ldw, pd->wrks2, pd->ldws2,
                       ipiv+j, pd);
         /* trsm */
         A0 = W + ((j*ldw + j) SHIFT);
         B0 = A0 + ((nl*ldw) SHIFT);
         my_trsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                 nl, nr, one, A0, ldw, B0, ldw);
      }
      ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, pd->id);
      /* gemm */
      nm = m;
      off = 0;
      if (rank == 0)
      {
         off = j+nl;
         nm -= (j+nl);
      }
      A0 = W + ((j*ldw + off) SHIFT);
      B0 = pd->wrks2[0] + (((j+nl)*pd->ldws2[0] + j) SHIFT);
      C0 = A0 + ((nl*ldw) SHIFT);
      if (rank < p)
         Mjoin(PATL, ammm)(CblasNoTrans, CblasNoTrans, nm, nr, nl,
                           none, A0, ldw, B0, pd->ldws2[0], one, C0, ldw);
      ATL_rtgetf2(rank, Jb, m, nr, p, j+nl, W, ldw, ipiv, pd);
      if (rank == 0)
      {
         /* left swap */
         ATL_laswp_blk(nr, nl, j+nl, j, W, ldw, pd->wrks2, pd->ldws2,
                       ipiv+j+nl, pd);
      }
   }
   else
   {
      int off = rank == 0 ? j : 0;
      pd->maxindx[rank] = -1;
      pd->stage[rank] = -1;
      A0 = W + ((j*ldw + off) SHIFT);
      ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, pd->id);
      if (rank < p)
         ATL_subtgetf2(rank, m-off, n, j, p, A0, ldw, ipiv+j, pd);
      ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, pd->id);
   }
}

static void ATL_btgetf2(int id, int Jb, int m00, int n, TYPE *W0,
                        TYPE *Aj, int lda, int *ipiv, ATL_GETRF_2D_t *pd)
{
   int rank = pd->rid;
   const int nb = pd->nb;
   const int pt = pd->pt;
   const int nbnb = nb*nb;
   const int active = Mmod(Jb, pt);
   const int nbpt = nb * pt;
   const int nbi = pd->nbi;
   ATL_CINT M = pd->M - Jb*nb, N = n, MN = Mmin(M, N);
   TYPE **WRKS = pd->wrks2;
   volatile ATL_INT *maxindx = pd->maxindx, *stage=pd->stage;
   int p = Mmin(pt, (M+nb-1)/nb);
   TYPE *a, *W, *Wc, *w, *v, *Ac;
   TYPE pivval, apv, apv2, pv2;
   TYPE *A0, *B0, *C0;
   int pivrank, m=m00;
   int m0, mr, nr;
   ATL_INT *ldws2 = pd->ldws2;
   ATL_INT locpiv, globpiv, k, j, i, ldw, nbildw, pj;
   #ifdef TCPLX
      const TYPE one[2] = {ATL_rone, ATL_rzero};
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #else
      const TYPE one = ATL_rone;
      const TYPE none = ATL_rnone;
   #endif
   if (p < 2)
   {
      ATL_getf2s(id, Jb, m00, n, W0, Aj, lda, ipiv, pd);
      return;
   }
   rank -= active;
   if (rank < 0) rank += pt; /* rank acts as vrank, with active as 0 */
   maxindx[rank] = -1;
   stage[rank] = -1;

/*
 * Make ldw's a multiple of 16 bytes that is not a power of 2; 0's ldw
 * is larger by mr than all other ldws (ldw1)
 */
   #if defined(DREAL) || defined(SCPLX)
      ldw = ((m+1)>>1)<<1;
      if (!(ldw & (ldw-1)))
         ldw += 2;
   #elif defined(SREAL)
      ldw = ((m+3)>>2)<<2;
      if (!(ldw & (ldw-1)))
         ldw += 4;
   #else
      ldw = m;
      if (!(ldw & (ldw-1)))
         ldw++;
   #endif
   w = W = WRKS[rank];
   ldws2[rank] = ldw;
   nbildw = nbi * ldw;
   if (pt > 1) ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, id);
   /* if (rank >= p) return; */

   /* everything is setup, now copy, then factorize */
   if (Jb) /* copy from C-major */
   {
      bcAblk2cmat(m, n, W0, myRBS, W, ldw, nb, 1, pd->ca2cc);
   }
   else /* copy from original column-major */
   {
      Mjoin(PATL, bcG2L_cpy)(m, n, Aj+((rank*nb) SHIFT), lda, W, ldw, nb, pt);
   }

   #ifndef NO_RECURSIVE_PANEL
      ATL_rtgetf2(rank, Jb, m, n, p, 0, W, ldw, ipiv, pd);
   #else
      pj = 0;
      nr = Mmod(n, nbi);
      for (j=0; j<(n-nr); j+=nbi, w+=(nbildw SHIFT))
      {
         /* getrfi */
         maxindx[rank] = -1;
         stage[rank] = -1;
         ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, id);
         if (rank < p)
            ATL_subtgetf2(rank, m-pj, nbi,j,p, w+(pj SHIFT), ldw, ipiv+j, pd);
         /* lswp */
         if (rank == 0)
         {
            /* ipiv is already adjusted for sub-panel */
            /* for (i=0; i<nbi; i++) ipiv[j+i] += j; */
            ATL_laswp_blk(nbi, j, j, 0, W, ldw, WRKS, ldws2, ipiv+j, pd);
         }
         if (j+nbi < N)
         {
            /* rswp n trsm */
            if (rank == 0)
            {
               ATL_laswp_blk(nbi, N-j-nbi, j, j+nbi, W, ldw,
                             WRKS, ldws2, ipiv+j, pd);
               B0 = w + ((nbildw+j) SHIFT);
               A0 = w + (j SHIFT);
               my_trsm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                  nbi, N-j-nbi, one, A0, ldw, B0, ldw);
               pj += nbi;
            }
            ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, id);
            /* gemm */
            A0 = w + (pj SHIFT);
            B0 = WRKS[0] + (((j+nbi)*ldws2[0] + j) SHIFT);
            C0 = w + ((nbildw + pj) SHIFT);
            if (rank < p)
               Mjoin(PATL, ammm)(CblasNoTrans,CblasNoTrans, m-pj,N-j-nbi,nbi,
                                 none, A0, ldw, B0, ldws2[0], one, C0, ldw);
         }
      }
      if (nr)
      {
         /* getrfi */
         maxindx[rank] = -1;
         stage[rank] = -1;
         ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, id);
         if (rank < p)
            ATL_subtgetf2(rank, m-pj, nr, j, p, w+(pj SHIFT), ldw, ipiv+j, pd);
         /* lswp */
         if (rank == 0)
         {
            /* ipiv is already adjusted for sub-panel */
            /* for (i=0; i<nbi; i++) ipiv[j+i] += j; */
            ATL_laswp_blk(nr, j, j, 0, W, ldw, WRKS, ldws2, ipiv+j, pd);
         }
      }
   #endif
   ATL_CBC2D_barrier(ATL_SYNC_COL, pd->cbc, id);

   m = m00;
   /* copy to A-major */
   if (Jb) W0 -= (nbnb SHIFT);
   if (!rank)
   {
      /* blk-maj for TRSM */
      Mjoin(PATL, gecopy)(nb, n, W, ldw, W0-(nbnb SHIFT), nb);
      bcRm2am(m-nb, n, W+(nb SHIFT), ldw, W0+myROFF, myRBS, nb, 1, pd->r2a);
   }
   else
   {
      bcRm2am(m, n, W, ldw, W0-(nbnb SHIFT)+myROFF, myRBS, nb, 1, pd->r2a);
   }
   #ifdef EARLY_A_COPY
      Ac = Aj + ((Jb*nb + rank*nb) SHIFT);
      Mjoin(PATL, bcL2G_cpy)(m00, n, W, ldw, Ac, lda, nb, pt);
   #endif
}

/*
 * Factors one panel.
 * Jb: gpanel no.
 * m: actual local m.
 * n: actual n for current panel.
 * W: must point to actual data.
 * Aj: points the beginning of the gpanel.
 * ipiv: already adjusted for current panel.
 */
static void ATL_getf2_CBC(int id, int Jb, int m, int n, TYPE *W,
                          TYPE *Aj, int lda, int *ipiv, ATL_GETRF_2D_t *pd)
{
   #if defined(USE_SGETF2)
      ATL_getf2s(id, Jb, m, n, W, Aj, lda, ipiv, pd);
   #else
      if (pd->pt==1)
         ATL_getf2s(id, Jb, m, n, W, Aj, lda, ipiv, pd);
      else
         ATL_btgetf2(id, Jb, m, n, W, Aj, lda, ipiv, pd);
   #endif
}


static void ATL_getrf_bcAmm_bk_2D(int id, ATL_GETRF_2D_t *pd)
{
   const int rid = pd->rid, cid = pd->cid;
   const int pt = pd->pt, qt = pd->qt, nt = pd->nt;
   const int nb = pd->nb;
   const int mu = pd->mu, nu = pd->nu;
   const int NMu = pd->NMu, NNu = pd->NNu;
   const int RMu = pd->RMu, RNu = pd->RNu;
   const int major = pd->major;
   const int aM = pd->M, aN = pd->N;
   const int MN = pd->MN;
   const int Rm = pd->Rm, Rn = pd->Rn;
   const int oM = pd->oM, oN = pd->oN;
   const int oMb = pd->oMb, oNb = pd->oNb;
   const int aMb = oMb-1, aNb = oNb-1;
   int *JB = pd->JB;
   TYPE *A = pd->A;
   const int lda = pd->lda;
   int *ipiv = pd->ipiv;
   TYPE **WRKS = pd->WRKS;
   int *ldws = pd->ldws;
   volatile int *PROG = (volatile int *)pd->PROG;
   ATL_CBC2D_t *cbc = pd->cbc;
   cm2am_t r2an = pd->r2a;
   cm2am_t c2a1 = pd->c2a;
   am2cm_t a2rn = pd->a2r;
   am2cm_t a2c1 = pd->a2c;
   ablk2cmat_t ca2cc = pd->ca2cc;
   cmat2ablk_t cc2ca = pd->cc2ca;
   ammkern_t amm1 = pd->amm_b1;
   const int nbnb = nb * nb;
   const int nbnb2 = nbnb*2;
   const int nbpt = nb * pt;
   const size_t nblda = nb*lda;
   const size_t nbqtlda = nblda *qt;
   int m, n, m_b, n_b, m0;
   void *vp;
   int nmu, nnu, rmu, rm, cn;
   int i, j, k;
   int ib, jb, kb, cjb, jbn, pcjb;
   int cJB, Jb, cJb, Jbn, Jbp;
   int pnl_selected;

   TYPE *W0, *Wj0, *Wj, *Wi, *Wk, *Wc;
   int ldw;
   size_t nbldw;
   TYPE *Ac, *Aj, *Aj0, *Ak;
   TYPE *A0, *B0, *C0, *nU;
   #ifdef TREAL
      TYPE one = ATL_rone, zero = ATL_rzero;
   #elif defined(TCPLX)
      const TYPE one[2]={ATL_rone,ATL_rzero}, zero[2]={ATL_rzero,ATL_rzero};
   #endif
   #ifdef ATL_LAPROF
      long long ltsc=0, ctsc=0, ttsc=0;
      rdtsc(ltsc);
   #endif

   /* Distribute rows and columns */
   NUMROC(m_b, rid, oMb, pt);
   NUMROC(n_b, cid, oNb, qt);
   m0 = m = m_b*nb;
   n = n_b*nb;
   ldws[id] = ldw = m + (3*nb);
   vp = malloc(ATL_MulBySize(((size_t)ldw)*n)+ATL_Cachelen);
   if (!vp) pd->info = NOT_ENOUGH_MEMORY;
   ATL_CBC2D_barrier(ATL_SYNC_GRID, cbc, id);
   if (id == 0)
   {
      for (i=0; i<nt; i++) if (pd[i].info) break;
      if (i < nt) for (j=0; j<nt; j++) pd[j].info = NOT_ENOUGH_MEMORY;
   }
   ATL_CBC2D_barrier(ATL_SYNC_GRID, cbc, id);
   if (pd->info) return;
   WRKS[id] = W0 = ATL_AlignPtr(vp);
   pd->m_bs[id] = m_b;
   pd->n_b = n_b;
   nbldw = nb * ldw;
   #ifdef USE_BCSWP
      pd->bp->larrs[rid] = W0 + ((nbnb + nbnb) SHIFT);
      pd->bp->lldps[rid] = (m_b+3)*nbnb;
   #endif

   nmu = NMu;
   pd->rmu = rmu = 0;
   rm = 0;
   if (Rm && (rid == Mmod(aMb,pt)))
   {
      rm = Rm;
      pd->rmu = rmu = RMu;
      m_b--;
      pd->m_bs[id]--;
      m0 = m0 + rm - nb;
   }

   /* Wj0 points to left-most panel's 1st block */
   for (jb=0, Jb=cid, Wj0=W0, Aj0=A+((cid*nblda) SHIFT); jb<n_b;
         jb++, Jb+=qt, Wj0+=(nbldw SHIFT), Aj0+=(nbqtlda SHIFT))
   {
      int active;
      while (1)
      {
         cjb = jb;
         cJb = Jb;
         cJB = *JB;
         /* find work */
         while ((cjb<n_b) && (PROG[cJb] >= 0) && (PROG[cJb] >= PROG[PROG[cJb]]))
         {
         #ifndef DISABLE_LOOKAHEAD
            cjb++;
            cJb += qt;
         #endif
         } /* found work */
         if (cjb >= n_b)
         {
            cjb = jb;
            cJb = Jb;
            while (PROG[cJb] >= PROG[PROG[cJb]]);
         }
         ATL_CBC2D_min(ATL_SYNC_COL, cbc, id, &cjb);
         cJb = BLKL2G(cjb, cid, qt);
         Wj = Wj0 + (((cjb-jb)*nbldw) SHIFT); /* Wj: points to current lpanel */
         Aj = Aj0 + (((cJb-Jb)*nblda) SHIFT); /* Aj: points to current gpanel */

         nnu = NNu;
         cn = nb;
         if(Rn && (cJb == aNb))
         {
            nnu = RNu;
            cn = Rn;
         }
         Wj += (nbnb SHIFT);    /* Move the pointer to leave the first block */

         if (PROG[cJb] < 0) /* need to copy this panel */
         {
            if (cJb != cJB) while (PROG[cJB] < 0); /* wait for active copy */
            #ifdef ATL_LAPROF
               gettime(ltsc, ctsc, Timers[id][CopyWait], ttsc);
            #endif
            CopyInPanel(cJb, m, rm, cn, Aj, lda, Wj, pd);
            #ifdef ATL_LAPROF
               gettime(ltsc, ctsc, Timers[id][Copy2L], ttsc);
            #endif
            if (pt > 1)
               ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
            if (PROG[cJb] < 0) PROG[cJb] = 0; /* race but same value */
         }
         else if (PROG[cJb] < PROG[PROG[cJb]]) /* we have more updates */
         {
            ApplyOneUpdate(cjb, cJb, PROG[cJb], nnu, cn, Wj, Aj, pd);
         }
         if (PROG[cJb] == cJb) /* time to do getf2 */
         {
            break; /* exit from this loop and apply getf2 */
         }
      } /* end of while(1) */
      /* now do the panel factorization */
      if (cJb)
      {
         ib = NEXT_GEMM_BLK(cJb-1, rid, pt);
         Wk = Wj + ((ib*nbnb + nbnb) SHIFT);
      }
      else
      {
         ib = 0;
         Wk = Wj;
      }
      i = ib*nb;
      active = Mmod(cJb, pt);
      #ifndef NO_GETF2
         ATL_getf2_CBC(id, cJb, m_b*nb+rm-i, cn, Wk, Aj, lda, ipiv+cJb*nb, pd);
         if (pd->info2) pd->info2 += cJb*nb;
      #else
         if (rid == active)
         {
            int I = cJb*nb;
            int tn = MN - I;
            tn = Mmin(tn, nb);
            for (i=0; i<tn; i++) ipiv[I+i] = i;
         }
      #endif
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][F2], ttsc);
      #endif
      /* sync to confirm panel is factorized */
      if (pt > 1)
         ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][CSync2], ttsc);
      #endif
      if (rid == active)
      {
         /* adjust ipiv if needed */
         int I = cJb*nb;
         int tn = MN - I;
         tn = Mmin(tn, nb);
         #ifdef USE_BCSWP
            ATL_bcIpivEncode(pd->bp, tn, I, I);
         #endif
         for (i=0; i<tn; i++) ipiv[I+i] += I;
         *JB = cJb+1;
         PROG[cJb] = *JB;
      }
      #if !defined(EARLY_A_COPY) \
         && !defined(NO_GETF2) && !defined(USE_SGETF2)
         if (pt > 1 && (aM-cJb*nb) > nb) /* copy wasn't done */
         {
            TYPE *tA = Aj + ((cJb*nb) SHIFT);
            int rank = rid - active;
            if (rank < 0) rank += pt;
            tA += ((rank*nb) SHIFT);
            if (cJb) Wk -= (nbnb SHIFT);
            i = ib*nb;
            if (rid == active)
            {
               int tm = Mmin(m_b*nb+rm-i, nb);
               Mjoin(PATL, gecopy)(tm, cn, Wk-(nbnb SHIFT), nb,
                                   tA, lda);
               bcAm2rm(m_b*nb+rm-i-tm, cn, Wk+myROFF, myRBS,
                       tA+(nbpt SHIFT), lda, nb, pt, pd->a2r);
            }
            else
            {
               bcAm2rm(m_b*nb+rm-i, cn, Wk-(nbnb SHIFT)+myROFF, myRBS,
                       tA, lda, nb, pt, pd->a2r);
            }
         }
         #ifdef ATL_LAPROF
            gettime(ltsc, ctsc, Timers[id][Copy2G], ttsc);
         #endif
      #endif
      #ifndef EARLY_B_COPY
         if (cJb) /* for later panels */
         {
            int last_active = active - 1;
            if (last_active < 0) last_active += pt;
            if (last_active == rid)
            {
               A0 = Aj + (((cJb-1)*nb) SHIFT);
               B0 = Wj + (((((cJb-1)/pt)-1)*nbnb) SHIFT);
               a2c1(nb, cn, one, A0, lda, B0);
            }
         }
         #ifdef ATL_LAPROF
            gettime(ltsc, ctsc, Timers[id][Copy2GU], ttsc);
         #endif
      #endif
      if (pt > 1)
         ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
      #ifdef ATL_LAPROF
         gettime(ltsc, ctsc, Timers[id][CSync2], ttsc);
      #endif
   } /* end of jb loop */
   ApplyLeftPivots(pd);
   #ifdef ATL_LAPROF
      gettime(ltsc, ctsc, Timers[id][LSWP], ttsc);
   #endif
   /* wait for everyone to finish before freeing memory */
   if (pt > 1)
      ATL_CBC2D_barrier(ATL_SYNC_COL, cbc, id);
   free(vp);
}


static void ATL_getrf_bcAmm_bk_tp(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *tpp = vpp;
   ATL_GETRF_2D_t *pd = tpp->PD;
   const int pt=pd->pt, qt=pd->qt, nt=pd->nt;
   int r = nt < ATL_NTHREADS ? vrank :
   #ifdef ATL_COLAFF
      /* col-major aff */
      GID(CID(rank, pt, qt), RID(rank, pt, qt), qt, pt);
   #else
      rank;
   #endif
   ATL_getrf_bcAmm_bk_2D(r, pd+r);
}

#ifdef ATL_CBC_NAME
   int Mjoin(PATL, tgetrf_bcAmm)(ATL_CINT major, ATL_CINT M, ATL_CINT N,
                                 TYPE *A, ATL_CINT lda, int *ipiv)
#else
   int Mjoin(PATL, tgetrf)(ATL_CINT major, ATL_CINT M, ATL_CINT N,
                           TYPE *A, ATL_CINT lda, int *ipiv)
#endif
{
   int ierr = 0, terr = 0;
   int i, j, k;
   int pt, qt, nt;
   int nb, nbi, oMb, oNb;
   int mu, nu, nmu, nnu;
   int rmu, rnu, rm, rn;
   int oM, oN, MN;
   cm2am_t r2a, c2a;
   am2cm_t a2r, a2c;
   ablk2cmat_t ca2cc;
   cmat2ablk_t cc2ca;
   ammkern_t amm_b1;
   ammkern_t amm_bn;

   void *vp;
   ATL_GETRF_2D_t* pd;
   TYPE **WRKS;
   int *ldws;
   int *m_bs;
   int *PROG;
   int JB = 0;
   int IDX = -1;
   ATL_CBC2D_t cbc;

   ATL_TGETF2_M_t *lups;
   ATL_INT *maxindx, *stage;
   TYPE **wrks2;
   int *ldws2;
   size_t ldw2, pan2sz;
   TYPE *T;
   void *vp2;
   int *ipiv2;
   #ifdef USE_BCSWP
      ATL_bcpiv_t *bp;
   #endif
   if (M < N) return CASE_NOT_HANDLED;

   #if 0
      SelectBlockAndGrid(M, N, &nb, &nbi, &pt, &qt);
   #else
      IDX = Mjoin(PATL,getrf_bcAmm_info)(M, N, &nb, &pt, &qt);
   #endif
   SelectGenKernels(IDX, &nb, &mu, &nu, &r2a, &c2a, &a2r, &a2c,
                    &ca2cc, &cc2ca, &amm_b1, &amm_bn);


   oMb = (M+nb-1)/nb;
   oNb = (N+nb-1)/nb;
   pt = Mmin(pt, oMb);
   qt = Mmin(qt, oNb);
   nt = pt * qt;

   if (nt < 2) return Mjoin(PATL, getrf)(major, M, N, A, lda, ipiv);

   oM = oMb * nb;
   oN = oNb * nb;
   MN = Mmin(M, N);
   nmu = nb / mu;
   nnu = nb / nu;
   rm = Mmod(M, nb);
   rn = Mmod(N, nb);
   rmu = (rm+mu-1)/mu;
   rnu = (rn+nu-1)/nu;

   vp = malloc( ((sizeof(ATL_GETRF_2D_t) + sizeof(TYPE*)
               + (sizeof(int)*2)) * nt)
               + (oNb * sizeof(int)) + 128);
   if (!vp) return NOT_ENOUGH_MEMORY;
   pd = ATL_AlignPtr(vp);
   WRKS = (TYPE**)(pd + nt);
   ldws = (int*)(WRKS + nt);
   m_bs = (int*)(ldws + nt);
   PROG = ATL_AlignPtr((int*)(m_bs + nt));
   for(i=0; i<oNb; i++) PROG[i] = -1;

   ATL_CBC2D_barrier_init(&cbc, pt, qt);

   ldw2 = ((oMb + pt + pt - 1) / pt) * pt * nb; /* 1 xtra blk for all */
   pan2sz = ( sizeof(ATL_TGETF2_M_t)       /* for the struct */
            + sizeof(ATL_INT)*2           /* for 2 state vars */
            + sizeof(TYPE*) ) * pt        /* for wrks2 ptrs */
            + (ATL_MulBySize(((size_t)(ldw2))*nb)   /* for wrks2 data */
            + sizeof(ATL_INT)             /* for ldws2 */
            + ATL_Cachelen ) * pt;
   vp2 = malloc(pan2sz*qt + sizeof(int)*nb);
   ipiv2 = (int*)(vp2+qt*pan2sz);
   if (!vp2) return NOT_ENOUGH_MEMORY;
   #ifdef USE_BCSWP
      bp = ATL_bcIpivInit(MN, ipiv, 1, pt, qt, oM, nb, mu, nu);
   #endif

   for (j=0; j<qt; j++)
   {
      wrks2 = (TYPE**)(vp2 + j*pan2sz);
      ldws2 = (ATL_INT*)(wrks2 + pt);
      lups = (ATL_TGETF2_M_t*)(ldws2 + pt);
      maxindx = (ATL_INT*)(lups + pt);
      stage = (ATL_INT*)(maxindx + pt);
      T = (TYPE*)(stage + pt);

      for (i=0, k=j; i<pt; i++, k+=qt)
      {
         pd[k].rid = i;
         pd[k].cid = j;
         pd[k].id = k;
         pd[k].pt = pt;
         pd[k].qt = qt;
         pd[k].nt = nt;

         pd[k].nb = nb;
         pd[k].nbi = nbi;
         pd[k].mu = mu;
         pd[k].nu = nu;
         pd[k].NMu = nmu;
         pd[k].NNu = nnu;
         pd[k].RMu = rmu;
         pd[k].RNu = rnu;

         pd[k].r2a = r2a;
         pd[k].c2a = c2a;
         pd[k].a2r = a2r;
         pd[k].a2c = a2c;
         pd[k].ca2cc = ca2cc;
         pd[k].cc2ca = cc2ca;
         pd[k].amm_b1 = amm_b1;
         pd[k].amm_bn = amm_bn; /* needed for complex */

         pd[k].M = M;
         pd[k].N = N;
         pd[k].MN = MN;
         pd[k].A = A;
         pd[k].lda = lda;
         pd[k].ipiv = ipiv;
         pd[k].WRKS = WRKS;
         pd[k].ldws = ldws;
         pd[k].m_bs = m_bs;
         pd[k].info = 0;

         pd[k].Rm = rm;
         pd[k].Rn = rn;
         pd[k].oM = oM;
         pd[k].oM = oM;
         pd[k].oN = oN;
         pd[k].oMb = oMb;
         pd[k].oNb = oNb;

         pd[k].JB = &JB;
         pd[k].PROG = PROG;
         pd[k].cbc = &cbc;

         pd[k].lups = lups;
         pd[k].maxindx = maxindx;
         pd[k].stage = stage;
         pd[k].wrks2 = wrks2;
         pd[k].ldw2 = ldw2;
         pd[k].ldws2 = ldws2;
         pd[k].info2 = 0;
         wrks2[i] = ATL_AlignPtr(T +
               i*(((size_t)(ldw2))*(nb SHIFT) + (ATL_Cachelen/sizeof(TYPE))));
         #ifdef USE_BCSWP
            pd[k].bp = bp+j;
         #endif
      }
   }

   ATL_goParallel(nt, ATL_getrf_bcAmm_bk_tp, NULL, pd, NULL);

   ATL_CBC2D_barrier_destroy(&cbc);
   /* first check for malloc error */
   ierr = pd[0].info;
   if (!ierr) /* no malloc error */
   {
      /* check for singular matrix error */
      for (i=0; i<qt; i++)
      {
         if (pd[i].info2 && terr)
            terr = Mmin(terr, pd[i].info2);
         else if (pd[i].info2)
            terr = pd[i].info2;
      }
      ierr = terr;
   }
   #ifdef USE_BCSWP
      free(bp);
   #endif

   free(vp);
   free(vp2);

   #if defined(ATL_LAPROF) && defined(SHOW_PROF)
      PrintProfilingInfo_pt(nt, NTimers, Timers, TimerNames, -1);
   #endif
   return (ierr);
}
