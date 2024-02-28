#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_tbitvec.h"
#include "atlas_gatmctr.h"
#include "atlas_tperf.h"
#include Mstr(Mjoin(ATLAS_PRE,amm_syrk.h))
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
#ifdef Conj_
   #define syrkBlk Mjoin(PATL,herkBlk)
   #define syrk_K Mjoin(PATL,herk_KLoop)
#else
   #define syrkBlk Mjoin(PATL,syrkBlk)
   #define syrk_K  Mjoin(PATL,syrk_KLoop)
#endif
#define DEBUG 0
#ifndef DEBUG
   #define DEBUG 0
#endif
void syrkBlk(ipinfo_t*,int,size_t,size_t,const TYPE *,cm2am_t,ablk2cmat_t,
             const SCALAR,TYPE*,TYPE*,TYPE*,TYPE*,TYPE*,TYPE*,TYPE*,TYPE*,
             TYPE*,TYPE*);
void syrk_K(ipinfo_t*, int flg, size_t d, const TYPE*A, cm2am_t,
            ablk2cmat_t, const SCALAR, TYPE *C, TYPE *rS, TYPE *wS, TYPE *wA,
            TYPE *wB, TYPE *rC, TYPE *wC, TYPE *wU);
/*
 * Translates C coordinate in lower triangular matrix to Cblk #.
 * nm is the number of Mblocks,
 * (i,j) is the coordinate of the block
 */
#define Mcoord2cblk(i_, j_, nm_) ((((nm_)+(nm_)-j-1)*j)>>1 - (j_) + (i_) - 1)
/*
 * Translates block number b_ to (i,j) coordinates assuming first panel
 * has nm_ blocks in it (including diagonal block)
 */
#define Mcblk2coord(NM_, B_, I_, J_) \
{ \
   register int n_ = (NM_)-1, b_=(B_), j_; \
   for (j_=0; b_ >= n_; j_++) \
   { \
      b_ -= n_; \
      n_--; \
   } \
   (J_) = j_; \
   (I_) = j_ + b_ + 1; \
}

static void DoBlkWtCpy
(
   unsigned const int rank,   /* my thread rank */
   ATL_tsyrk_ammN_t *pd,      /* problem def structure */
   unsigned const long dblk,  /* diagonal block of C to compute */
   unsigned long kblk,        /* kblk to compute */
   unsigned const int flg,    /* 0: amm_b0: else: amm_b1 */
   TYPE *rS,                  /* real wrk for SYRK's A (unused for real) */
   TYPE *iS,                  /* imag wrk for SYRK's A (real for real */
   TYPE *rC,                  /* real wrk for C (unused for real) */
   TYPE *iC                   /* imag wrk for C (real for real) */
)
{
   ipinfo_t *ip=pd->ip;
   const TYPE *A;             /* ptr to original matrix */
   TYPE *wA, *wB;             /* ptr to copy space */
   #ifdef TCPLX
      TYPE *rA, *rB;
      ammkern_t amm_bn, amm_b1, ammR, ammC;
   #else
      ammkern_t amm;          /* amm kern to use  */
   #endif
   void *KdonCtr = ATL_AddBytesPtr(pd->KdonCtr, dblk*pd->KDCgap);
   const size_t skp=pd->pansz0, ndiag=pd->ndiag, lda=ip->lda;
   size_t i, k;
   const int UPPER=flg&1;
   const int kb0=ip->kb;
   int kb = kb0, KB=kb0, b=ip->nb;
   int kk, nmu=ip->nmu, nnu=ip->nnu;
/*
 * K cleanup needs to use special kern
 */
   if (!kblk)
   {
      #ifdef TCPLX
         if (flg&4)
            ammR = ammC = ip->ammK1_b0;
         else
         {
            ammR = ip->ammK1_bn;
            ammC = ip->ammK1_b1;
         }
         amm_b1 = ip->ammK1_b1;
         amm_bn = ip->ammK1_bn;
      #else
         amm = (flg&4) ? ip->ammK1_b0 : ip->ammK1_b1;
      #endif
      KB = ip->KB0;
      kb = ip->kb0;
   }
   else
   {
      #ifdef TCPLX
         if (flg&4)
            ammR = ammC = ip->amm_b0;
         else
         {
            ammR = ip->amm_bn;
            ammC = ip->amm_b1;
         }
         amm_b1 = ip->amm_b1;
         amm_bn = ip->amm_bn;
      #else
         amm = (flg&4) ? ip->amm_b0 : ip->amm_b1;
      #endif
   }
   if (dblk == ndiag-1)
   {
      nmu = ip->nmuF;
      nnu = ip->nnuF;
      b = ip->nF;
   }
/*
 * If last K block of different size, take it from end of vector to avoid
 * screwing up alignment
 */
   if (kb0 != ip->kb0)
      kblk = (kblk) ? kblk-1 : ip->nfkblks;
/* fprintf(stderr, "%d: dblk=%d, kblk=%d\n", rank, dblk, kblk);*/
   A = IdxA_ip(ip, pd->A, dblk, kblk);
   if (UPPER)
   {
      if (dblk != ndiag-1)
         wA = IdxAw_ip(ip, pd->wA, dblk, kblk);
      else
         wA = iS;
      if (dblk)
         wB = IdxBw_ip(ip, pd->wAt, kblk, dblk) - skp;
      else
         wB = iS;
   }
   else         /* Lower */
   {
      if (dblk != ndiag-1)
         wB = IdxBw_ip(ip, pd->wAt, kblk, dblk);
      else
         wB = iS;
      if (dblk)
          wA = IdxAw_ip(ip, pd->wA, dblk, kblk) - skp;
      else
         wA = iS;
   }
   #if DEBUG > 1
      fprintf(stderr, "%d: DIAMM(%ld,%ld), wA=%p, wB=%p\n", rank, dblk, kblk,
              wA, wB);
   #endif

   #ifdef TCPLX
      if (dblk < ip->nfnblks)
      {
         rA = wA + ip->szA;
         rB = wB + ip->szB;
      }
      else
      {
         rA = wA + ip->pszA;
         rB = wB + ip->pszB;
      }
      ip->a2blk(kb, b, ip->alpA, A, lda, rA, wA);
      ip->b2blk(kb, b, ip->alpB, A, lda, rA, wB);
   #else
      ip->a2blk(kb, b, ip->alpA, A, lda, wA);
      ip->b2blk(kb, b, ip->alpB, A, lda, wB);
   #endif
/*
 * Signal copy of this K-panel of A/A' is complete iff I just copied the
 * last block
 */
   if (ATL_gatmctr_dec(KdonCtr, rank) == 1)
      ATL_tSetBitBV(pd->cpydonBV, dblk);
   #ifdef TCPLX
      ammR(nmu, nnu, KB, wA, wB, rC, rA, wB, iC);
      ammC(nmu, nnu, KB, rA, wB, iC, rA, rB, iC);
      amm_bn(nmu, nnu, KB, rA, rB, rC, wA, rB, iC);
      amm_b1(nmu, nnu, KB, wA, rB, iC, wA, rB, iC);
   #else
      amm(nmu, nnu, KB, wA, wB, iC, wA, wB, iC);
   #endif
}

static void WriteAmmDiag
   (ATL_tsyrk_ammN_t *pd, const long d, const int B, const SCALAR beta,
    const TYPE *W)
/*
 * Copies from BxB block-major W to upper/lower portion of pd->C
 */
{
   ipinfo_t *ip=pd->ip;
   TYPE *c;
   #ifdef TCPLX
      const SCALAR ONE=ip->ONE;
   #else
      #define ONE ATL_rone
   #endif
   c = IdxC_ip(ip, pd->C, d, d);
   if (pd->flg & 1)  /* Upper */
   {
      int j;
      const size_t ldc=(ip->ldc)SHIFT;

      for (j=0; j < B; j++, c += ldc, W += (B SHIFT))
      {
         Mjoin(PATL,axpby)(j+1, ONE, W, 1, beta, c, 1);
         #ifdef Conj_
            c[j+j+1] = ATL_rzero;
         #endif
      }
   }
   else              /* Lower */
   {
      const size_t ldp1=(ip->ldc+1)SHIFT;
      int j;

      for (j=0; j < B; j++, c += ldp1, W += (B+1)SHIFT)
      {
         Mjoin(PATL,axpby)(B-j, ONE, W, 1, beta, c, 1);
         #ifdef Conj_
            c[1] = ATL_rzero;
         #endif
      }
   }
}
#ifndef TCPLX
   #undef ONE
#endif

static void DoBlksWtCpy
(
   unsigned const int rank,   /* my thread rank */
   ATL_tsyrk_ammN_t *pd,      /* problem def structure */
   unsigned const long dblk,  /* diagonal block of C to compute */
   unsigned long kctr,        /* non-zero Kbeg ctr */
   TYPE *rS,                  /* real wrk for SYRK's A (unused for real) */
   TYPE *iS,                  /* imag wrk for SYRK's A (real for real */
   TYPE *rC,                  /* real wrk for C (unused for real) */
   TYPE *iC,                  /* imag wrk for C (real for real) */
   TYPE *U                    /* Upper reflection space */
)
{
   ipinfo_t *ip=pd->ip;
   ablk2cmat_t syblk2c=pd->syblk2c_b1;
   cm2am_t sya2blk = pd->sya2blk;
   const size_t ldc = ip->ldc;
   size_t skpA, skpB;
   #ifdef TCPLX
      const TYPE *beta=ip->ONE;
      #define ONE ip->ONE
      const TYPE ZERO[2] = {ATL_rzero, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
      TYPE beta = ONE;
   #endif
   ATL_UINT szC=ip->szC;
   const unsigned int flg = pd->flg;
   unsigned int nb;
   const size_t nkblks=pd->nkblks, ndiag=pd->ndiag;
   long kblk = nkblks - kctr, b;
   void *KbegCtr = ATL_AddBytesPtr(pd->KbegCtr, dblk*pd->KBCgap);
   void *KdonCtr = ATL_AddBytesPtr(pd->KdonCtr, dblk*pd->KDCgap);
   void *cpydonBV = pd->cpydonBV;
   TYPE *wA0=pd->wA, *wB0=pd->wAt, *wA, *wB;
   const TYPE *A = pd->A;
   void *Cdmut;

   if (sya2blk)
   {
      if (flg & 1) /* Upper */
      {
         skpA = 0;
         skpB = pd->pansz0;
         if (dblk == ndiag-1)
            wA0 = NULL;
         if (!dblk)
            wB0 = NULL;
      }
      else         /* Lower */
      {
         skpA = pd->pansz0;
         skpB = 0;
         if (dblk == ndiag-1)
            wB0 = NULL;
         if (!dblk)
            wA0 = NULL;
      }
      wA = (wA0) ? (IdxAw_ip(ip, wA0, dblk, kblk) - skpA) : NULL;
      wB = (wB0) ? (IdxBw_ip(ip, wB0, kblk, dblk) - skpB) : NULL;
      syrkBlk(ip, flg|4, dblk, kblk, A, sya2blk, NULL, ZERO, NULL,
              rS, iS, wA, wA, wB, wB, rC, iC, U);
      if (ATL_gatmctr_dec(KdonCtr, rank) == 1)
      {
         ATL_tSetBitBV(cpydonBV, dblk);
         #if DEBUG > 1
            ATL_tPrintBV(stderr, "00", cpydonBV);
         #endif
      }
   }
   else
      DoBlkWtCpy(rank, pd, dblk, kblk, flg|4, rS, iS, rC, iC);
   while ((kctr = ATL_gatmctr_dec(KbegCtr, rank)))
   {
      long kk;
      kblk = nkblks - kctr;

      if (sya2blk)
      {
         wA = (wA0) ? (IdxAw_ip(ip, wA0, dblk, kblk) - skpA) : NULL;
         wB = (wB0) ? (IdxBw_ip(ip, wB0, kblk, dblk) - skpB) : NULL;

         syrkBlk(ip, flg, dblk, kblk, A, sya2blk, NULL, ZERO, NULL,
                 rS, iS, wA, wA, wB, wB, rC, iC, U);
         kk = ATL_gatmctr_dec(KdonCtr, rank);
         if (kk == 1)
            ATL_tSetBitBV(cpydonBV, dblk);
      }
      else
         DoBlkWtCpy(rank, pd, dblk, kblk, flg, rS, iS, rC, iC);
   }
/*
 * For upper or amm, write to U ptr before seizing mutex.  After getting mutex,
 * we will reflect/copy this matrix to Upper one
 */
   if (!sya2blk)  /* copy from access-major iC to block-major U */
   {
      nb = (dblk < ip->nfnblks) ? ip->nb : ip->pnb;
      #ifdef TCPLX
         pd->syblk2c(nb, nb, ip->alpC, rC, iC, ZERO, U, nb);
      #else
         pd->syblk2c(nb, nb, ip->alpC, iC, ATL_rzero, U, nb);
      #endif
   }
   else if (flg&1)  /* Upper must use workspace */
      syrkBlk(ip, flg, dblk, 0, NULL, NULL, pd->syblk2c, ZERO, NULL,
              NULL, NULL, NULL, NULL, NULL, NULL, rC, iC, U);
/*
 * Now, seize mutex for diagonal block of original C, and copy back out
 * only above/below diagonal
 * NOTE: because I'm using bitvecs, have to use tSetBitBV, which contains
 *       another mutex lock.  If I kept betaBV in a char array instead,
 *       could reduce to only needing Cdmut here.
 */
   Cdmut = ATL_AddBytesPtr(pd->Cdmuts, dblk*pd->dCinc);
   ATL_lock(Cdmut);
   if (!ATL_SetBitBV(pd->dbetaBV, dblk))   /* I'm first, so tell rest of  */
   {                                       /* thrs don't apply beta */
      beta = pd->beta;                     /* because I'm going to do it */
      syblk2c = pd->syblk2c;
   }
   if (sya2blk)
   {
      if (flg&1)  /* upper just needs to reflect from U to C */
         syrkBlk(ip, flg, dblk, 0, NULL, NULL, NULL, beta, pd->C,
                 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, U);
      else /* Lower writes directly to C from rC,iC */
         syrkBlk(ip, flg, dblk, 0, NULL, NULL, syblk2c, beta, pd->C,
                 NULL, NULL, NULL, NULL, NULL, NULL, rC, iC, NULL);
   }
   else
      WriteAmmDiag(pd, dblk, nb, beta, U);
   ATL_unlock(Cdmut);
}
#undef ONE

/*
 * In this phase, we work only on nshar diagonal blocks, while copying both
 * A & A'.  For init diag work, we parallelize both N & K dims so that the copy
 * is done as quickly as possible.  Threads coming in first choose differing
 * diag blks; diagonal blocks are dealt out using the diCtr local
 * counter (which starts at nshar).
 * Once all shared diagonal blocks are dealt out, new threads will start using
 * the atomic ctr array KbegCtr array to share K work for each diagonal.
 * both KbegCtr & KdonCtr are nshar-len arrays of atomic counters.  Each
 * counter starts at nkblks.  Once the block pointed to by KbegCtr is
 * completely copied, the copying array increments the KdonCtr.  Only one
 * core per diag will get KdonCtr == 1 after doing his copy, and this
 * core will to set the appropriate bit in cpydonBV, which is a nnblks-length
 * global threaded bit vector.  If the kth bit is set, that means the
 * kth row of A & kth col of A' has been copied.
 * Once a thread gets KbegCtr for a particular diag of 0, it means there's
 * no more work for this block of C, and so it will seize the appropriate
 * Cdmuts mutex which protects each diagonal block of C, and write its
 * finished contribution out to C.  The first such thread to ever seize
 * the mutex will scope dbetaBV to find if this diagonal block needs beta
 * applied, while later threads will use beta=1.
 * Eventually, all init diagonal work is finished, and the first processor to
 * get 0 for all dCtr & KbegCtr requests will set NOINCPY=1, so later
 * threads don't have to query all the counters to know they should proceed
 * to working on non-init diagonals or non-diagonal block.
 */
static void DoSharDiag(const int P, const int rank, ATL_tsyrk_ammN_t *pd,
                       TYPE *rS, TYPE *iS, TYPE *rC, TYPE *iC, TYPE *U)
{
   const size_t KBinc = pd->KBCgap;
   void *KbegCtr;
   const unsigned long nshar=pd->nshar;
   unsigned long jstart, k;
   int DIAG=1;

   if (nshar >= P)
      jstart = rank * (nshar/P);
   else
      jstart = rank - (rank/nshar)*nshar;

   #if DEBUG > 1
      fprintf(stderr, "%d: start DoSharDiag\n", rank);
   #endif
   while (!(pd->NOINICPY))
   {
      long d=0;
/*
 *    Find which diagonal block to work on, and then which k blk to use
 */
      if (DIAG)
      {
         d = ATL_atmctr_dec(pd->diCtr);
         if (d)
         {
            k = nshar - d;
            KbegCtr = ATL_AddBytesPtr(pd->KbegCtr, k*KBinc);
            k = ATL_gatmctr_dec(KbegCtr, rank);
            if (!k)     /* if no more K work to do */
               d = 0;   /* can't work on this diag after all */
         }
      }
/*
 *    If all diagonal blocks currently being worked on by threads, find
 *    one that I can help with.
 */
      if (!d)
      {
         unsigned long i;
         DIAG = 0;
         for (i=0; i < nshar; i++)
         {
            unsigned long j = i+jstart;
            j = (j < nshar) ? j : j-nshar;
            KbegCtr = ATL_AddBytesPtr(pd->KbegCtr, j*KBinc);
            k = ATL_gatmctr_dec(KbegCtr, rank);
            d = nshar-j;
            if (k)
              goto FOUNDDK;
         }
         pd->NOINICPY = 1;    /* no initial work left to assign */
         #if DEBUG > 1
            fprintf(stderr, "%d: DONE  DoSharDiag\n", rank);
         #endif
         return;             /* so done with copy of A/A' & this func */
      }
/*
 *    If I reach here, I've got a valid d & k;  and I'll call a routine
 *    that continues to grab blocks from this diag along K until all K
 *    is done; it will then write the answer back to the original C, and
 *    return to this loop to see if it can help with another diag.
 */
      FOUNDDK:
         DoBlksWtCpy(rank, pd, nshar-d, k, rS, iS, rC, iC, U);
   }
   #if DEBUG > 1
      fprintf(stderr, "%d: DONE  DoSharDiag\n", rank);
   #endif
}

static void Do1Diag(const int rank, ATL_tsyrk_ammN_t *pd, long d,
                    TYPE *rS, TYPE *iS, TYPE *rC, TYPE *iC, TYPE *U)
/*
 * This routine computes 1 diagonal block using SYRK
 */
{
   ipinfo_t *ip=pd->ip;
   TYPE *wA=pd->wA, *wB=pd->wAt;
   const size_t skp=pd->pansz0;
   const unsigned long ndiag=pd->ndiag;
   const unsigned int flg=pd->flg;
   #if DEBUG > 1
      fprintf(stderr, "%d: start Do1Diag %ld\n", rank, d);
   #endif
   if (flg&1)
   {
      if (d != ndiag-1)
         wA = IdxAw_ip(ip, wA, d, 0);
      else
         wA = NULL;
      if (d)
         wB = IdxBw_ip(ip, wB, 0, d) - skp;
      else
         wB = NULL;
   }
   else         /* Lower */
   {
      if (d != ndiag-1)
         wB = IdxBw_ip(ip, wB, 0, d);
      else
         wB = NULL;
      if (d)
          wA = IdxAw_ip(ip, wA, d, 0) - skp;
      else
         wA = NULL;
   }
   #if DEBUG > 1
      fprintf(stderr, "%d:C(%ld,%ld), wA=%p, wB=%p\n", rank, d, d, wA, wB);
   #endif
   if (pd->sya2blk)
      syrk_K(ip, flg, d, pd->A, pd->sya2blk, pd->syblk2c, pd->beta, pd->C,
             rS, iS, wA, wB, rC, iC, U);
   else  /* using gemm instead of SYRK! */
   {
      const int B0=ip->nb, B = (d < ip->nfnblks) ? B0 : ip->nF;
      const SCALAR beta=pd->beta;
      const TYPE *lA;
      #ifdef TCPLX
         const TYPE ONE[2]={ATL_rone,ATL_rzero}, ZERO[2]={ATL_rzero,ATL_rzero};
      #else
         #define ONE ATL_rone
         #define ZERO ATL_rzero
      #endif
      TYPE *c;
      int MV=3;

      if (!wA)
      {
         wA = iS;
         MV = 2;
      }
      else if (!wB)
      {
         MV = 1;
         wB = iS;
      }
      lA = IdxA_ip(ip, pd->A, d, 0);
      Mjoin(PATL,iploopsK)(ip, d, d, lA, lA, NULL, MV, wA, wB, rC, iC,
                           pd->beta, NULL);
/*
 *    Copy from access-major rC to block-major U, then update Upper/lower Cblk
 */
      #ifdef TCPLX
         pd->syblk2c(B, B, ip->alpC, rC, iC, ZERO, U, B);
      #else
         pd->syblk2c(B, B, ip->alpC, iC, ZERO, U, B);
      #endif
      c = IdxC_ip(ip, pd->C, d, d);
      if (pd->flg & 1)  /* Upper */
      {
         int j;
         const size_t ldc=(ip->ldc)SHIFT;

         for (j=0; j < B; j++, c += ldc, U += (B SHIFT))
         {
            Mjoin(PATL,axpby)(j+1, ONE, U, 1, beta, c, 1);
            #ifdef Conj_
               c[j+j+1] = ATL_rzero;
            #endif
         }
      }
      else              /* Lower */
      {
         const size_t ldp1=(ip->ldc+1)SHIFT;
         int j;

         for (j=0; j < B; j++, c += ldp1, U += (B+1)SHIFT)
         {
            Mjoin(PATL,axpby)(B-j, ONE, U, 1, beta, c, 1);
            #ifdef Conj_
               c[1] = ATL_rzero;
            #endif
         }
      }
   }
   ATL_tSetBitBV(pd->cpydonBV, d);
   if (!ATL_tInfoBV(pd->cpydonBV, ATL_TBV_NUNSET))
      pd->cpydone = 1;
   #if DEBUG > 1
      fprintf(stderr, "%d: DONE  Do1Diag %ld\n", rank, d);
   #endif
}
#ifndef TCPLX
   #undef ONE
   #undef ZERO
#endif

static int compPossCols
(
   void *gcpanDonBV,
   ATL_BV_t *lcpyDonBV, /* set bits are ready to to be computed */
   ATL_BV_t *lColBV     /* OUTPUT: ready cols still to be computed */
)
/*
 * Computes all column coordinates that can legally be computed given
 * cpyDonBV ready input blks, subtracting out any already-finished colpans
 * indicated by gcpanDonBV
 * RETURNS: nunset in gcpanDonBV
 *
 */
{
   long nunset, N, Ne, Nr, i;

   nunset = ATL_tGlb2locBV(lColBV, gcpanDonBV, 0);
   N = lColBV[0];
   Ne = N >> shBV;
/*
 * We can only compute those columns that are: (1) not already computed &
 * (2) have the Kpanel copied.
 */
   for (i=0; i < Ne; i++)
      lColBV[i+1] = (~(lColBV[i+1])) & lcpyDonBV[i+1];
   Nr = N & modmskBV;
   if (Nr)
      lColBV[i+1] = ((~(lColBV[i+1])) & lcpyDonBV[i+1]) & ((1L<<Nr)-1);

   return(nunset);
}

static int compPossRows
(
   const int UPPER,
   long ndiag,
   long J,              /* colpan to compute rows of */
   void *gcblkBV,
   ATL_BV_t *lcpyDonBV, /* set bits are ready to to be computed */
   ATL_BV_t *lrowBV     /* OUTPUT: rows in col J that can be computed */
)
/*
 * RETURNS: 0 if no rows possible, else 1.
 */
{
   long ret=0;
   if (UPPER)  /* lower blks don't exist, j'th blk diag */
   {
      ATL_BV_t N=(ndiag+bpiBV-1)>>shBV, k;
      if (!J)        /* no non-diag blks in 1st Upper colpan */
         return(0);
      k = J>>shBV;
      lrowBV[k+1] = allsetBV;
      if (!ATL_tGlb2locBV(lrowBV, gcblkBV, 0))  /* if no unset bits in cblkBV */
         return(0);                   /* col J has no computable blocks left! */
      for (k++; k < N; k++)
         lrowBV[k+1] = allsetBV;
   }
   else /* Lower, upper j blks don't exist, j'th blk diag */
   {
      ATL_BV_t *lp = lrowBV+1, N=(J+bpiBV)>>shBV, k;
      if (J == ndiag-1)  /* no non-diag blks in last colpan */
         return(0);
      for (k=0; k < N; k++)
         lp[k] = allsetBV;
      k = ATL_tGlb2locBV(lrowBV, gcblkBV, J+1);
      if (!k)            /* if no unset bits in cblkBV */
         return(0);      /* col J has no computable blocks left! */
   }
/*
 * To be a computable row block, cblkBV must not be set (not already done)
 * and lcpyDone must be set (or input K-panel not yet copied)
 */
   {
      const ATL_BV_t n=ndiag>>shBV, nr=ndiag&modmskBV;
      ATL_BV_t i;

      for (i=1; i <= n; i++)
      {
         ATL_BV_t tmp;
         tmp = (~lrowBV[i]) & lcpyDonBV[i];
         ret |= tmp;
         lrowBV[i] = tmp;
      }
      if (nr)
      {
         ATL_BV_t tmp;
         tmp = ((~lrowBV[i]) & lcpyDonBV[i])&((1L<<nr)-1);
         ret |= tmp;
         lrowBV[i] = tmp;
      }
   }
   return(ret ? 1:0);
}

static long FindColWithCount
   (const long ndiag, const int UPPER, long p0, const long maxpos,
    void *gCpanDonBV, void *gcblkBV0, size_t strd, ATL_BV_t *lcpyDonBV,
    ATL_BV_t *lcolsBV, ATL_BV_t *lrowsBV)
{
   long J, p=p0-1, KEEPON = 0;

   SEARCH_JS:
      if (++p == maxpos)
         return(-1);
      J = ATL_FindFirstSetBitBV(lcolsBV, p);
      J = (J < maxpos) ? J : -1;
      if (J != -1)   /* possible candidate */
      {
         if (ATL_tIsBitSetBV(gCpanDonBV, J))  /* if this J already done */
         {
            p = J;                            /* goto next */
            goto SEARCH_JS;                   /* and search again */
         }
      }
/*
 *    If we've got a candidate col, compute lrowsBV and see if we can do work!
 */
      if (J != -1)
      {
         void *gcblkBV;
         gcblkBV = ATL_AddBytesPtr(gcblkBV0,J*strd);
         if (!compPossRows(UPPER, ndiag, J, gcblkBV, lcpyDonBV, lrowsBV))
         {
            if (J != maxpos-1)  /* if more blks to search */
            {
               p = J;           /* start search from following J */
               goto SEARCH_JS;  /* and keep looking */
            }
            else                /* no blocks left */
               J = -1;          /* stop search with failure */
         }
      }
/*
 * If there are unsearched bits prior to p0, and we have so far not found
 * a candidate column, try again starting from 0
 */
   if (p0 && J == -1) /* if unsearched bits prior to p0, and we failed */
      return(FindColWithCount(ndiag, UPPER, 0, p0, gCpanDonBV, gcblkBV0,
                              strd, lcpyDonBV, lcolsBV, lrowsBV));
   return(J);
}
/*
 * This function is does non-diagonal work when A/A^T have been copied,
 * and switches to diagonal work/copying when it can't find ready work
 * to do.  The idea is to reduce max bus load by intermixing copy and
 * no-copy code.  It returns once the copy is complete, so we can go
 * to lower-overhead code that doesn't need to check copy dependencies.
 */
static int DoDepBlks(const int rank, ATL_tsyrk_ammN_t *pd,
                     TYPE *rS, TYPE *iS, TYPE *rC, TYPE *iC, TYPE *U)
{
   ipinfo_t *ip=pd->ip;
   const size_t stride = pd->LOCgap, Cstride=pd->Cgap;
   size_t skpA, skpB;
   ATL_BV_t *lcolsBV=ATL_AddBytesPtr(pd->locBVs,(rank+(rank<<1))*stride);
   ATL_BV_t *lrowsBV=ATL_AddBytesPtr(lcolsBV, stride);
   ATL_BV_t *lcpyDonBV=ATL_AddBytesPtr(lrowsBV, stride);
   void *gcpyDonBV = pd->cpydonBV, *gCpanDonBV=pd->cpanDonBV;
   void *gcblkBV0=pd->cblkBV;
   TYPE *wA = pd->wA, *wB = pd->wAt, *C = pd->C;
   const unsigned long ndiag=pd->ndiag;
   long nun, p0, J, nunCpy=(-2), nunCpan=(-2);
   #ifdef TCPLX
      const TYPE *beta=pd->beta;
   #else
      const TYPE beta=pd->beta;
   #endif
   unsigned long mycnt=0;
   const ablk2cmat_t blk2c=ip->blk2c;
   const int UPPER = pd->flg & 1;
   int RECALC;

   lcolsBV[0] = lrowsBV[0] = lcpyDonBV[0] = ndiag;
   p0 = ATL_tGetLocalBoundsBV(gcpyDonBV, rank, NULL);
   if (UPPER)  /* Upper skips 1st colpan B */
   {
      skpA = 0;
      skpB = pd->pansz0;
   }
   else        /* Lower skips 1st rowpan A */
   {
      skpA = pd->pansz0;
      skpB = 0;
   }
   if (!pd->NODWORK && rank < pd->ncpDiag)  /* If I'm a thread asked to */
      goto DODIAG;                          /* copy before computing dep */
   DEPLOOP:
      if (pd->cpydone || pd->DONE)  /* if copying or full computation done */
         return(mycnt);             /* this routine is complete */
/*
 *    We must recalculate computable cols if the number of copied blocks
 *    or the number of finished columns changes
 */
      nun = ATL_tInfoBV(gcpyDonBV, ATL_TBV_NUNSET);
      if (ndiag-nun < 2)   /* Must have 2 kpans copied to do anything! */
         goto DODIAG;
      if (nun == 0)  /* if all blocks copied, go to faster code */
      {
         pd->cpydone = 1;
         return(mycnt);   /* that skip dependency checks! */
      }
      if (nun == nunCpy)  /* number of copied blocks hasn't changed */
      {
         nun = ATL_tInfoBV(gCpanDonBV, ATL_TBV_NUNSET);
         if (!nun)          /* if no Cpans left undone */
         {
            pd->DONE = 1;
            return(mycnt);  /* all non-diag work complete, so leave */
         }
         RECALC = (nun != nunCpan);
         nunCpan = nun;
      }
      else
      {
         RECALC = 1;      /* # of copied blks changed, so recalc! */
         nunCpy = ATL_tGlb2locBV(lcpyDonBV, gcpyDonBV, 0);
         if (nunCpy == 0)
         {
            pd->cpydone = 1;
            return(mycnt);
         }
      }
      if (RECALC)
         nunCpan = compPossCols(gCpanDonBV, lcpyDonBV, lcolsBV);
      J = FindColWithCount(ndiag, UPPER, p0, ndiag, gCpanDonBV, gcblkBV0,
                           Cstride, lcpyDonBV, lcolsBV, lrowsBV);
      if (J == -1)    /* if no cols are computable */
         goto DODIAG; /* go copy more A/A^T by computing diag blks */
/*
 *    If we've got a column, work on it until out of good row blocks
 */
      {
         void *gcblkBV;
         long I=0, DIDSOME=0;

         gcblkBV = ATL_AddBytesPtr(gcblkBV0, J*Cstride);
         while ((I = ATL_FindFirstSetBitBV(lrowsBV, I)) != -1)
         {
            const long II = UPPER ? I : I-J-1;
            if (!ATL_tSetBitBV(gcblkBV, II)) /* 0=reserved it */
            {
               TYPE *wa, *wb;
               wa = IdxAw_ip(ip, wA, I, 0) - skpA;
               wb = IdxBw_ip(ip, wB, 0, J) - skpB;
               Mjoin(PATL,iploopsK)(ip, I, J, NULL, NULL, IdxC_ip(ip, C, I, J),
                                    3, wa, wb, rC, iC, beta, blk2c);
               #if DEBUG > 1
                  fprintf(stderr, "%d: C(%ld,%ld), wA=%p, wB=%p\n",
                          rank, I, J, wa, wb);
               #endif
               #if DEBUG
                  mycnt++;
               #endif
               DIDSOME = 1;
            }
            if (++I == ndiag)
               break;
            if (pd->cpydone)
               return(mycnt);
         }
/*
 *       If I computed some blocks, see if this column is finished
 *       If column is finished, see if whole problem is finished.
 *
 *       This is unsafe, and can result in not marking completed column.
 *       This is OK, because eventually we'll exit this function due to
 *       the copying finishing, and DoNonDiag will mark column safely.
 */
         if (DIDSOME)
         {
            if (!ATL_tInfoBV(gcblkBV, ATL_TBV_NUNSET))
            {
               ATL_tSetBitBV(gCpanDonBV, J);
               if (!ATL_tInfoBV(gCpanDonBV, ATL_TBV_NUNSET))
               {
                  pd->DONE = 1;
                  return(mycnt);
               }
            }
         }
         else             /* someone bogarted all my rows, try reducing */
            goto DODIAG;  /* contention by producing more possible cols */
      }
   goto DEPLOOP;
DODIAG:
   if (!pd->NODWORK)
   {
      J = ATL_gatmctr_dec(pd->dCtr, rank);
      if (J)
         Do1Diag(rank, pd, ndiag-J, rS, iS, rC, iC, U);
      else
         pd->NODWORK = 0;
   }
   goto DEPLOOP;
}

static long DoNonDiag
(
   const int rank,       /* worker # 0 <= rank < P */
   ATL_tsyrk_ammN_t *pd, /* problem definition */
   const unsigned int P, /* # of workers; rank < P */
   TYPE *rC,             /* real C, same as iC for real */
   TYPE *iC              /* imag C, same as rC for real */
)
{
   ipinfo_t *ip=pd->ip;
   const unsigned long ndiag=pd->ndiag;
   const int UPPER = (pd->flg & 1);
   TYPE *wA = pd->wA, *wB = pd->wAt;
   TYPE *C = pd->C;
   void *cpanDonBV = pd->cpanDonBV;
   long mycnt=0;
   #ifdef TCPLX
      const TYPE *beta=pd->beta;
   #else
      const TYPE beta=pd->beta;
   #endif
   const ablk2cmat_t blk2c=ip->blk2c;
#if 0   /* ignore this until the high-overhead cleanup debugged */
/*
 * First make one pass through cols, using only multiples of my rank
 */
   if (ndiag >= P)
   {
      unsigned int ngrab=1, ngot;
      unsigned long msk=1, old, i;
      for (j=0; j < ndiag; j += P)
      {
         const long c=j+rank;
         long r;
         void *cblkBV;
         if (c >= ndiag)
            break;
         cblkBV = ATL_AddBytesPtr(pd->cblkBV, pd->Cgap*c);

         if (ngrab != pd->ngrab)
         {
            ngrab = pd->ngrab;
            msk = (ngrab < bpiBV) ? ((1L<<ngrab)-1):allsetBV;
         }
/*
 *       Attempt to get ngrab rowblks of C from this c's cblkBV starting @ r
 */
         r = ATL_tFindUnsetBitBV(cblkBV, rank);
         if (r == -1)
         {
            if (!ATL_tIsBitSetBV(cpanDonBV, c))
               ATL_tSetBitBV(cpanDonBV, c);
            continue;
         }
         ngot = ngrab;
         old = ATL_tSetRangeBV(cblkBV, &ngot, r, msk);
         if (ngot)
         {
            unsigned int k;
            r = (UPPER) ? r : (1+c+r);
            for (k=0; k < ngot; k++)
            {
               if ((msk>>k)&(1L))
               {
                  Mjoin(PATL,iploopsK)(ip, k+r, c, NULL, NULL,
                     IdxC_ip(ip, C, r+k, c), 3, IdxAw_ip(ip, wA, r+k, 0),
                     IdxBw_ip(ip, wB, 0, c), rC, iC, beta, blk2c);
               }
            }
         }
      }
   }
#endif
/*
 * Now, all cores to into mode where they steal work from each other,
 * taking only 1 Cblk at a time.  Threads start with differing cblkBV regions
 * to minimize contention.
 */
   pd->ngrab = 1;  /* tell other cores to stop bogarting the Cblks */
   if (!pd->DONE)
   {
      void *cblkBV;
      const long b = ndiag / P, bs=b*rank, ndi=(UPPER) ? ndiag : ndiag-1;
      long c;
      size_t skpA, skpB;

      if (UPPER)  /* Upper skips 1st colpan B */
      {
         skpA = 0;
         skpB = pd->pansz0;
      }
      else        /* Lower skips 1st rowpan A */
      {
         skpA = pd->pansz0;
         skpB = 0;
      }
      while ( (c=ATL_tFindUnsetBitBV(cpanDonBV, rank)) != -1)
      {
         long r;
         TYPE *wa, *wb;
         int DIDSOME=0;

         if (pd->DONE)
            return(mycnt);
         cblkBV = ATL_AddBytesPtr(pd->cblkBV, pd->Cgap*c);
         DOROWS:
            r = ATL_tSetUnsetBitBV(cblkBV, rank);
            if (r < 0)
            {
               if (!ATL_tIsBitSetBV(cpanDonBV, c))
                  ATL_tSetBitBV(cpanDonBV, c);
               continue;
            }
            r = (UPPER) ? r : 1+c+r;
            wa = IdxAw_ip(ip, wA, r, 0) - skpA;
            wb = IdxBw_ip(ip, wB, 0, c) - skpB;
            #if DEBUG > 1
               fprintf(stderr, "%u:Nondiag(%ld,%ld), wA=%p, wB=%p\n",
                       rank, r, c, wa, wb);
            #endif
            Mjoin(PATL,iploopsK)(ip, r, c, NULL, NULL, IdxC_ip(ip, C, r, c),
                                 3, wa, wb, rC, iC, beta, blk2c);
            DIDSOME = 1;
            mycnt++;
         goto DOROWS;
      }
   }
   return(mycnt);
}

static void DoWorkN(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tsyrk_ammN_t *pd = pp->PD;
   long ndiag=0, ndep=0, nindep=0;
   const unsigned int P = pp->nworkers, szC = Mmax(pd->ip->szC,pd->szCs);
   TYPE *iS = pd->wC + vrank*pd->wrksz, *rS, *rC, *iC, *U;

   iS = ATL_AlignPtr(iS);
   #ifdef TCPLX
      rS = iS + pd->szS;
      iC = rS + pd->szS;
      iC = ATL_AlignPtr(iC);
      rC = iC + szC;
      U  = rC + szC;
   #else
      rS = iS;
      iC = rS + pd->szS;
      rC = iC = ATL_AlignPtr(iC);
      U  = iC + szC;
   #endif
   if (pd->sya2blk && !(pd->flg & 1))
      U = NULL;

   if (!(pd->NOINICPY))
      DoSharDiag(P, vrank, pd, rS, iS, rC, iC, U);
   if (!pd->cpydone && !pd->DONE)
      ndep = DoDepBlks(vrank, pd, rS, iS, rC, iC, U);
   if (!pd->DONE)
      nindep = DoNonDiag(vrank, pd, P, rC, iC);
   #if DEBUG
      fprintf(stdout, "%d: Cblks=%ld+%ld, ndiag=%ld,%d Dblks=%ld(%d)\n",
              vrank, ndep, nindep, pd->ndiag, pd->nshar, ndiag, pd->ncpDiag);
   #endif
}
/*
 * SYRK where all parallelism comes from blocks of N, built atop amm directly
 *    if (TA == AtlasNoTrans)
 *       C = alpha * A*A' + beta*C
 *    else
 *       C = alpha * A'*A + beta*C
 *    C is an upper or lower symmetric NxN matrix,
 *    A is a dense rectangular NxK (NoTrans) or KxN (Trans) matrix
 * RETURNS: 0 if operation performed, non-zero otherwise.
 *   This routine assumes it can copy all of A up-front to simplify parallelism.
 *   Will return non-zero if memory cannot be allocated, on the assumption
 *   it is called from recursive implementation that can recur until malloc
 *   succeeds.
 */
#ifdef Conj_
   #define tsyrk_amm_N Mjoin(PATL,therk_amm_N)
#else
   #define tsyrk_amm_N Mjoin(PATL,tsyrk_amm_N)
#endif
int tsyrk_amm_N
(
   const enum ATLAS_UPLO Uplo,
   const enum ATLAS_TRANS TA,
   ATL_CINT N,
   ATL_CINT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CINT lda,
   const SCALAR beta,
   TYPE *C,
   ATL_CINT ldc
)
{
   void *vp=NULL, *vp0;
   size_t szA, szC, szKpan, sz, extra, szThr, szAll;
   size_t ndiag, ncblks;
   double stim;
   unsigned int i, P, idx, nb, nbS, kbS, szS, szPO, nshar=0, ncpDiag, n;
   unsigned int szBETABV=0, DISYRK=0;
   ATL_tsyrk_ammN_t pd;
   ipinfo_t ip;
   #ifdef ATL_AVX
      #define symul 1.34
   #else
      #define symul 1.2
   #endif
   #define parpen (2.0*ATL_tstart_sec)
   #ifdef Conj_
      const enum ATLAS_TRANS TB = (TA == AtlasNoTrans) ?
                                  AtlasConjTrans : AtlasNoTrans;
   #else
      const enum ATLAS_TRANS TB = (TA == AtlasNoTrans) ?
                                  AtlasTrans : AtlasNoTrans;
   #endif
/*
 * Demand at least 2 large blocks along some dimension
 */
   n = ATL_sqAMM_LASTNB<<1;
   n = Mmin(200, n);
   if (N < n && K < n)
      goto DO_SERIAL;
   pd.flg = (Uplo == AtlasUpper);
   pd.flg |= (TA == AtlasNoTrans) ? 2 : 0;
   #if 1
      if (K >= (ATL_sqAMM_LASTKB<<3))
         idx = ATL_sqAMM_NCASES-1;
      else
         idx = Mjoin(PATL,GetSyrkIdx)(pd.flg, N, K, symul);
   #else
      idx = 0;
      nb = 4;
   #endif
   stim = Mjoin(PATL,sSyrkTimeEst)(idx, pd.flg, N, K, symul);
   if (stim <= parpen+ATL_tstart_sec+4.0*ATL_tstartgap_sec)
      goto DO_SERIAL;
   nb = Mjoin(PATL,sqGetAmmInfoInt)('K', idx);
   P = ATL_NTHREADS;
   if (N <= nb)  /* only possible parallelism is along K */
   {
      nshar = 1;                      /* share all K blocks */
      szKpan = (K+nb-1) / nb;
      P = Mmin(szKpan, ATL_NTHREADS); /* limit parallism to k blocks */
      goto GOT_NB;
   }
   ndiag = nb*ATL_NTHREADS;
   if (ndiag < N)            /* is there enough paralleism along N? */
      goto GOT_NB;
   if (ndiag < K)  /* is there enough parallism from K? */
   {
      nshar = 1;   /* means make all blks shared */
      goto GOT_NB;
   }
/*
 * See if we can create enough parallelism by sharing nshar initial
 * K-panels, which will unleash at least ncblks parallelism afterwords
 */
   szKpan = (K+nb-1) / nb;
   ndiag = (N+nb-1) / nb;
   if (nshar == 1)
      nshar = ndiag;
   ncblks = ((ndiag-1)*ndiag)>>1;
   if (szKpan*ndiag >= P)
   {
      nshar = P / szKpan;
      if (nshar < 2)      /* need at least 2 panels complete */
         nshar = 2;       /* before non-diag comp begins */
      goto GOT_NB;
   }
/*
 * If we reach here, need to try reducing nb to increase parallelism;
 * If this doesn't work, reduce parallelism.
 */
   nshar = 1;         /* get max parallism from K */
   nbS = nb;
   n = ATL_lcm(Mjoin(PATL,sqGetAmmInfoInt)('m', idx),
               Mjoin(PATL,sqGetAmmInfoInt)('n', idx));
   if (n <= nb)
   {
      P = szKpan;
      goto GOT_NB;
   }
   nb = (nb/n)*n;
   TRY_ANOTHER_NB:
      ndiag = N / nb;
      ncblks = ((ndiag-1)*ndiag)>>1;
      if (ndiag*szKpan+ncblks >= ATL_NTHREADS)
         goto GOT_NB;
      nbS = nb;
      nb -= n;
   if (nb > 24)
      goto TRY_ANOTHER_NB;

   nb = nbS;
   ndiag = N / nb;
   ncblks = ((ndiag-1)*ndiag)>>1;
   P = ndiag*szKpan+ncblks;
GOT_NB:
   if (P > 1)
   {
      double trem, ptim;
      trem = stim;
      ptim = parpen + ATL_tstart_sec;
      for (i=0; i < P; i++)
      {
         double t1;
         t1 = (i+1)*ATL_tstartgap_sec;
         if (trem > t1)
         {
            trem -= t1;
            ptim += ATL_tstartgap_sec;
         }
         else
         {
            ptim += trem/(i+1);
            break;
         }
      }
      ptim += trem / P;
      if (stim <= 1.2*ptim) /* don't accept less than 20% win */
         goto DO_SERIAL;    /* since parallel times hugely variable */
      else
         P = (ptim < stim) ? i : 1;
      #if DEBUG > 1
         fprintf(stdout, "   predicted speedup = %.2f (%d)\n", stim/ptim, i);
      #endif
      if (P < 1)
         goto DO_SERIAL;
      Mjoin(PATL,sqFillInIPInfo)(&ip, idx, TA, TB, N, N, K, lda, lda, ldc,
                                 alpha, beta, nb);
      ndiag = ip.nfnblks + ip.npnblks;
      ncblks = ((ndiag-1)*ndiag)>>1;
      ncpDiag = 0;
/*
 * If N <= NB, always use SYRK, in this case all computation done by SYRK.
 * -> This means we can only copy A once, while doing half the computation.
 * If we do SYRK for diagonals, but GEMM for off-diagonals, we must make
 * 3 copies of of A: 2 for GEMM, one for SYRK (if SYRK's copy is the same
 * as GEMM's A or B copy, this isn't true, but this is only rarely the case,
 * since SYRK always uses K-major wt MU=NU, which is a loser on modern x86).
 * -> Cost of extra copy is N*K*(Mt) (Mt = time to copy from mem)
 * -> Diagonal computation savings is nb*N*K (4*nb*N*K for cplx)
 * ==> Use syrk when NB > Mt; we estimate Mt as 16*P
 *     -> compute speed scales perfectly with P, but Mt doesn't
 */
      #ifdef TCPLX
         DISYRK = (N <= nb || nb > (ATL_NTHREADS<<2));
      #else
         DISYRK = (N <= nb || nb > (ATL_NTHREADS<<4));
      #endif
   }
   if (P < 2)
   {
DO_SERIAL:
      #ifdef Conj_
         Mjoin(PATL,herk)(Uplo, TA, N, K, *alpha, A, lda, *beta, C, ldc);
      #else
         Mjoin(PATL,syrk)(Uplo, TA, N, K, alpha, A, lda, beta, C, ldc);
      #endif
      return(0);
   }
   #if DEBUG
   printf("%u: N=%d,%d (%d*%d+%d) F=(%d,%d;%d,%d), kb=%d SY=%d, nshar=%ld\n",
          P, (int)N, (int)K, (int)ip.nfnblks, ip.nb, ip.pnb, ip.mF, ip.nF,
          ip.nmuF, ip.nnuF, ip.kb, DISYRK, nshar);
   #endif
/*
 * Set up parallel data structure
 */
   szKpan = ip.nfnblks ? ip.szA : ip.pszA;
   szKpan *= ip.nfkblks+1;
   if (ndiag > 1)
   {
      szA = szKpan * (ndiag-1);
      szAll = ATL_MulBySize(szA + szA) + 2*ATL_Cachelen;
   }
   else
      szA = szAll = 0;
   extra = (ip.mu<<1)*ip.nu;
/*
 * Each thread needs: spc C (max(SYRK,GEMM)) + spc for U (if Upper) + 1 blk syrk
 */
   nb = (ip.nfnblks) ? ip.nb : ip.pnb;
   szThr = ip.szC;
   if (DISYRK)   /* need extra workspace if doing SYRK */
   {
      extra = Mmax(extra, (ATL_SYRKK_NU+ATL_SYRKK_NU)*ATL_SYRKK_NU);
      nbS = ((nb+ATL_SYRKK_NU-1)/ATL_SYRKK_NU)*ATL_SYRKK_NU;
      kbS = ((ip.kb+ATL_SYRKK_KU-1)/ATL_SYRKK_KU)*ATL_SYRKK_KU;
      szC = ((nbS+1)*nbS)>>1; /* only need lower tri blks, not full nnu*nnu */
      szC *= ((ATL_SYRKK_NU*ATL_SYRKK_NU+ATL_SYRKK_VLEN-1)/ATL_SYRKK_VLEN)
             * ATL_SYRKK_VLEN;
      if (Uplo == AtlasUpper)
         extra = Mmax(extra, ip.szC);
      szS = nbS * kbS;
      szS = ((szS+ATL_SYRKK_VLEN-1)/ATL_SYRKK_VLEN)*ATL_SYRKK_VLEN;
      szThr = Mmax(szThr, szC);         /* spc for amm/syrk C */
      szThr = ATL_MulBySize(szS + szThr + extra) + ATL_Cachelen + ATL_sizeof-1;
      szThr = ATL_DivBySize(szThr);
   }
   else
   {
      szS = (ip.nfnblks) ? ip.szA : ip.pszA;
      szC = (ip.nfnblks) ? ip.nb : ip.pnb;
      szC *= szC;
      szC = Mmax(szC, extra);
      szThr += szS + szC;
      szC = 0;
      szThr += ATL_DivBySize(ATL_Cachelen+ATL_Cachelen);
   }
   pd.dCgap = pd.Cgap = ATL_tSizeofBV(ndiag, P);
   pd.KBCgap = ATL_gatmctr_sizeof(P, ATL_GAC_PUB);
   pd.KDCgap = pd.KBCgap; /* ATL_gatmctr_sizeof(P, ATL_GAC_PUB); */
   pd.dCinc = sizeof(ATL_lock_t);
   pd.dCinc = (size_t) ATL_AlignSafeLS(pd.dCinc);
   pd.nkblks = ip.nfkblks + 1;
   pd.LOCgap = (ndiag+bpiBV-1) >> shBV;
   pd.LOCgap = (pd.LOCgap+1)*sizeof(ATL_BV_t);
   pd.LOCgap = (size_t) ATL_AlignSafeLS(pd.LOCgap);
   pd.diCgap = ATL_atmctr_sizeof;
   if (nshar)  /* do I need space for nshar-based work? */
   {
      szBETABV = (((nshar+bpiBV-1)>>shBV)+1)*sizeof(ATL_BV_t);
      szBETABV = (size_t) ATL_AlignSafeLS(szBETABV);
      szPO = nshar*(pd.KBCgap+pd.KDCgap+pd.dCinc+pd.diCgap) + szBETABV;
   }
   else
      szPO = 0;
   szPO += ndiag*pd.Cgap;                    /* sz of cblkBV */
   szPO += pd.KBCgap;                        /* dCtr */
   szPO += pd.dCgap;                         /* cpydonBV */
   szPO += (P+(P<<1))*pd.LOCgap;             /* each core has 3 locBVs */
   szPO += szBETABV;                         /* dbetaBV */
   sz = ATL_MulBySize(szAll + P*szThr) + szPO + ATL_SAFELS;

   if (sz <= ATL_PTMAXMALLOC)
      vp = malloc(sz);
   if (!vp)                              /* if I'm over malloc limit */
      return(2);                            /* return and recur on K */

   pd.wrksz = szThr SHIFT;
   pd.KbegCtr  = ATL_AlignSafeLS(vp);
   pd.KdonCtr  = ATL_AddBytesPtr(pd.KbegCtr, pd.KBCgap*nshar);
   pd.cpydonBV = ATL_AddBytesPtr(pd.KdonCtr, pd.KDCgap*nshar);
   pd.locBVs   = ATL_AddBytesPtr(pd.cpydonBV, pd.dCgap);
   pd.cblkBV   = ATL_AddBytesPtr(pd.locBVs, (P+(P<<1))*pd.LOCgap);
   pd.Cdmuts   = ATL_AddBytesPtr(pd.cblkBV, pd.dCgap*ndiag);
   pd.diCtr    = ATL_AddBytesPtr(pd.Cdmuts, pd.dCinc*nshar);
   pd.dCtr     = ATL_AddBytesPtr(pd.diCtr, pd.diCgap);
   pd.dbetaBV  = ATL_AddBytesPtr(pd.dCtr, pd.KBCgap);
   if (ndiag > 1)
   {
      pd.wA    = ATL_AddBytesPtr(pd.dbetaBV, szBETABV);
      pd.wA    = ATL_AlignPtr(pd.wA);
      pd.wAt   = pd.wA + (szA SHIFT);
      pd.wAt   = ATL_AlignPtr(pd.wAt);
      pd.wC    = pd.wAt + (szA SHIFT);
   }
   else
   {
      pd.wA    = pd.wAt = NULL;
      pd.wC    = ATL_AddBytesPtr(pd.dbetaBV, szBETABV);
   }
   pd.wC       = ATL_AlignPtr(pd.wC);

   pd.ip = &ip;
   pd.ndiag = ndiag;
   pd.nshar = nshar;
   pd.ncpDiag = ncpDiag;
   pd.ncblks = ncblks;
   pd.beta = beta;
   pd.szCs = szC;
   pd.szS = szS;
   pd.NODWORK = 0;
   pd.DONE = (ndiag > 1) ? 0 : 1;
   pd.ngrab = 4;  /* adjust this later */

   if (DISYRK)
   {
      pd.sya2blk=IS_COLMAJ(TA)?Mjoin(PATL,a2blk_syrkT):Mjoin(PATL,a2blk_syrkN);
      if (Uplo == AtlasLower)
      {
         pd.pansz0 = (size_t)IdxAw_ip(&ip, NULL, 1, 0);
         if (SCALAR_IS_NONE(alpha))
         {
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_aNb1);
            if (SCALAR_IS_ONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_aNb1);
            else if (SCALAR_IS_NONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_aNbN);
            else
               pd.syblk2c = SCALAR_IS_ZERO(beta) ?
                  Mjoin(PATL,SyrkIntoC_aNb0) : Mjoin(PATL,SyrkIntoC_aNbX);
         }
         else if (SCALAR_IS_ONE(alpha))
         {
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_a1b1);
            if (SCALAR_IS_ONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_a1b1);
            else if (SCALAR_IS_NONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_a1bN);
            else
               pd.syblk2c = SCALAR_IS_ZERO(beta) ?
                  Mjoin(PATL,SyrkIntoC_a1b0) : Mjoin(PATL,SyrkIntoC_a1bX);
         }
         else
         {
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_aXb1);
            if (SCALAR_IS_ONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_aXb1);
            else if (SCALAR_IS_NONE(beta))
               pd.syblk2c = Mjoin(PATL,SyrkIntoC_aXbN);
            else
               pd.syblk2c = SCALAR_IS_ZERO(beta) ?
                  Mjoin(PATL,SyrkIntoC_aXb0) : Mjoin(PATL,SyrkIntoC_aXbX);
         }
      }
      else
      {
         pd.pansz0 = (size_t)IdxBw_ip(&ip, NULL, 0, 1);
         if (SCALAR_IS_NONE(alpha))
         {
            pd.syblk2c = Mjoin(PATL,SyrkIntoC_aNb0);
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_aNb1);
         }
         else if (SCALAR_IS_ONE(alpha))
         {
            pd.syblk2c = Mjoin(PATL,SyrkIntoC_a1b0);
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_a1b1);
         }
         else
         {
            pd.syblk2c = Mjoin(PATL,SyrkIntoC_aXb0);
            pd.syblk2c_b1 = Mjoin(PATL,SyrkIntoC_aXb1);
         }
      }
   }
   else
   {
      int ialp;

      if (Uplo == AtlasLower)
         pd.pansz0 = (size_t)IdxAw_ip(&ip, NULL, 1, 0);
      else
         pd.pansz0 = (size_t)IdxBw_ip(&ip, NULL, 0, 1);
      if (SCALAR_IS_NONE(ip.alpC))
         ialp = -1;
      else
         ialp = SCALAR_IS_ONE(ip.alpC) ? 1 : 2;
      pd.syblk2c = Mjoin(PATL,sqGetAmmInfoPtr)(idx, 1, ialp, 0);
      pd.sya2blk = NULL;
   }
   pd.pansz0 /= sizeof(TYPE);
   pd.NODWORK = (ndiag > nshar) ? 0 : 1;
   pd.NOINICPY = !nshar;
   pd.cpydone = 0;
   pd.A = A;
   pd.lda = lda;
   pd.C = C;
   pd.ldc = ldc;
/*
 * Initialize all parallel overhead structures
 */
   {
      void *vp1;

      for (vp0=pd.Cdmuts,i=0; i < nshar; i++,vp0=ATL_AddBytesPtr(vp0,pd.dCinc))
         ATL_lock_init(vp0);

      vp0 = pd.KbegCtr;
      vp1 = pd.KdonCtr;
      for (i=0; i < nshar; i++)
      {
         ATL_gatmctr_init(vp0, P, pd.nkblks, ATL_GAC_PUB);
         ATL_gatmctr_init(vp1, P, pd.nkblks, ATL_GAC_PUB);
         vp0 = ATL_AddBytesPtr(vp0, pd.KBCgap);
         vp1 = ATL_AddBytesPtr(vp1, pd.KDCgap);
      }
      if (ndiag > nshar)
         ATL_gatmctr_init(pd.dCtr, P, ndiag-nshar, ATL_GAC_PUB);
      if (nshar)
         ATL_atmctr_init(pd.diCtr, nshar);

      vp0 = pd.cblkBV;
      for (i=0; i < ndiag; i++)
      {
         int n = (Uplo == AtlasUpper) ? i : ndiag-i-1;
         ATL_tInitBV(vp0, n, P);
         vp0 = ATL_AddBytesPtr(vp0, pd.Cgap);
      }
      ATL_tInitBV(pd.cpydonBV, ndiag, P);
      ATL_InitBV(nshar, pd.dbetaBV, 0);
/*
 *    cpanDonBV is set to the colpan with no non-diag work: 1st for Upper,
 *    last for Lower
 */
      if (Uplo == AtlasUpper)
      {
         pd.cpanDonBV = pd.cblkBV;
         ATL_tInitBV(pd.cpanDonBV, ndiag, P);
         ATL_tScopeBitBV(pd.cpanDonBV, 0, 514);
      }
      else  /* Uplo == AtlasLower */
      {
         pd.cpanDonBV = ATL_AddBytesPtr(pd.cblkBV, (ndiag-1)*pd.Cgap);
         ATL_tInitBV(pd.cpanDonBV, ndiag, P);
         ATL_tScopeBitBV(pd.cpanDonBV, ndiag-1, 514);
      }
   }

   ATL_goParallel(P, DoWorkN, NULL, &pd, NULL);
   #ifdef DEBUG2
      ATL_assert(ATL_FindFirstUnsetBitBV(pd.cblkBV, 0) == -1);
   #endif

/*
 * Free allocated structures and return;
 */
   for (vp0=pd.Cdmuts,i=0; i < nshar; i++,vp0=ATL_AddBytesPtr(vp0,pd.dCinc))
      ATL_lock_destroy(vp0);
   free(vp);
   return(0);
}

#ifdef Conj_
   #define tsyrk_amm Mjoin(PATL,therk_amm)
#else
   #define tsyrk_amm Mjoin(PATL,tsyrk_amm)
#endif
void tsyrk_amm
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CSZT N,
    ATL_CSZT K, const SCALAR alpha, const TYPE *A, ATL_CSZT lda,
    const SCALAR beta, TYPE *C, ATL_CSZT ldc)

{
   #ifdef TCPLX
      const TYPE ONE[2]={ATL_rone, ATL_rzero};
      size_t kmul = (Trans==AtlasNoTrans||Trans==AtlasConj) ? lda+lda : 2;
   #else
      #define ONE ATL_rone
      size_t kmul = (Trans==AtlasNoTrans) ? lda : 1;
   #endif
/*
 *  Recur on K until tsyrk_amm can allocate enough space
 */
   if (tsyrk_amm_N(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc))
   {
      unsigned int kL = K>>1, kR=K-kL;
      const TYPE *a = A + kL*kmul;
/*
 *    This stopping criteria should never happen, but it's here in case
 *    we have a system where you can't malloc much of anything, where we'll
 *    just try serial
 */
      if (kL <= ATL_sqAMM_LASTKB)
      {
         Mjoin(PATL,syrk)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
         return;
      }
      tsyrk_amm(Uplo, Trans, N, kL, alpha, A, lda, beta, C, ldc);
      tsyrk_amm(Uplo, Trans, N, kR, alpha, a, lda, ONE,  C, ldc);
   }
}
#ifndef TCPLX
   #undef ONE
#endif
