/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 Siju Samuel
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 */

#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_ptalias_lapack.h"
#include "atlas_lapack.h"
#include "atlas_lvl3.h"
#include "atlas_qrrmeth.h"
/* #include Mstr(Mjoin(ATLAS_PRE,oprk_perf.h)) */
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#if ATL_VWopgen_MAX_KB <= 88
   #define ATL_QRKB ATL_VWopgen_100KB
   #define ATL_QRLCM ATL_VWopgen_100LCMMN
#elif ATL_VWopgen_99KB <= 88
   #define ATL_QRKB ATL_VWopgen_99KB
   #define ATL_QRLCM ATL_VWopgen_99LCMMN
#else
   #define ATL_QRKB ATL_VWopgen_98KB
   #define ATL_QRLCM ATL_VWopgen_98LCMMN
#endif

#ifdef ATL_USEPTHREADS
   #include "atlas_threads.h"
   #include "atlas_taffinity.h"
   #include "atlas_tcacheedge.h"
#else
   #include "atlas_cacheedge.h"
#endif


int ATL_gerqr(ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT LDA, TYPE  *TAU,
               TYPE *ws_RQ2, TYPE *ws_T, ATL_CINT LDT,
               TYPE *WORKM, const int buildT)
{
   int top, bottom, buildT_temp;
   int topMN;
   int I, INFO, IINFO, lbuilt, rbuilt, method;
   int LDA2 = LDA SHIFT;                    /* for complex LDA *2             */
   int LDT2 = LDT SHIFT;                    /* for complex LDT *2             */
   ATL_CINT minMN = Mmin(M, N);

   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #else
      #define ONE ATL_rone
   #endif

   if (M < 1 || N < 1) return(0);           /* Nothing to do.                 */
   METHOD(method, N, M, LDA);                   /* Find the method.           */
   #if !defined(ATL_USEPTHREADS)
   if (method == 2 || method == 3) method=1;    /* Don't PCA if no affinity.  */
   #endif

   switch(method)                               /* Based on method;           */
   {
      case 0:  /* RECURSION. */

      /*
       * Choose a smart recursive column partitioning based on M:
       */
         if (minMN >= 2*ATL_QRKB) /* big prob, put remainder on right */
         {
            topMN = ((minMN>>1)/ATL_QRKB)*ATL_QRKB;
            bottom = minMN - topMN;
            top  = M - bottom;
         }
         else /* small prob, keep M mult of MU (MU more critical than NU)     */
         {
            bottom = ((minMN>>1)/ATL_QRLCM)*ATL_QRLCM;
            topMN = minMN - bottom;
            top = M - bottom;
         }

         if (top==0 || bottom==0)      /* If too small for that,              */
         {
            bottom = (minMN>>1);       /* Stop trying to be clever.           */
            topMN = minMN - bottom;
            top = M - bottom;
         }

      /*----------------------------------------------------------------------*/
      /* On the bottom half, we use the same workspaces.                      */
      /* Because we know we have a top hand side we must always               */
      /* build T, so we can multiply by Q before doing the top side.          */
      /*----------------------------------------------------------------------*/
         ATL_gerqr(bottom,N, (A+(top SHIFT)), LDA, (TAU+(topMN SHIFT)), ws_RQ2,
                   ( ws_T+(topMN SHIFT)+topMN*LDT2), LDT, WORKM, 1);

      /*----------------------------------------------------------------------*/
      /* Now we must adjust the top hand side according to our T.             */
      /* We must apply H'                                                     */
      /*----------------------------------------------------------------------*/

         ATL_larfb(CblasRight, CblasNoTrans, LABackward, LARowStore,
                   top, N, bottom, (A +(top SHIFT))  , LDA,
                   (ws_T+(topMN SHIFT)+topMN*LDT2), LDT, A, LDA, WORKM, M);

      /*----------------------------------------------------------------------*/
      /* On the top  half,                                                    */
      /*----------------------------------------------------------------------*/
         ATL_gerqr(top, N-bottom,(A), LDA, (TAU),
                   ws_RQ2,
                   ( ws_T), LDT, WORKM, buildT);

      /*----------------------------------------------------------------------*/
      /* If we build T, the bottom side must be completely built, and         */
      /* the top side should be partially built. We need to fill in           */
      /* the lower  left  hand block, 'bottom' rows by 'top' columns.         */
      /* The formula is -T2 * (Y2 * Y1^T) * T1.                               */
      /* The routine is in ATL_larft.c.                                       */
      /*----------------------------------------------------------------------*/

         if (buildT )
         {
            ATL_larft_block(LABackward, LARowStore,
                            N, minMN, minMN-bottom, bottom,
                            A+((M -minMN) SHIFT), LDA, ws_T, LDT);
         }

         return(0);

      case 1: /* SERIAL (single core mode) */
         if (minMN >= 4)
         {
            Mjoin(PATL,gemoveT)(N, minMN, ONE, A+((M-minMN) SHIFT),LDA,WORKM,N);
            ATL_geql2(N, minMN, WORKM, N, TAU, ws_RQ2);
            Mjoin(PATL,gemoveT)(minMN,N,ONE,WORKM, N, A+((M-minMN) SHIFT), LDA);
            /* make conjugate  of TAU */
            #ifdef TCPLX
               Mjoin(PATLU,scal)(minMN, ATL_rnone, TAU+1, 2);
            #endif
         }
         else
         {
            ATL_gerq2(minMN, N, A+((M -minMN) SHIFT) , LDA, TAU, ws_RQ2);
         }

         if (buildT  || M > minMN)
         {
            ATL_larft(LABackward, LARowStore, N, minMN,
                      A+((M -minMN) SHIFT), LDA, TAU, ws_T, LDT);
         }
         break; /* END CASE */

      #if defined(ATL_USEPTHREADS)  /* Cases 2 & 3 only for parallel. */
      case 2: /* PCA COPY */
      case 3: /* PCA NOCOPY (but does not exist for RQ) */
         if (buildT || (M > minMN) )
         {
            buildT_temp = 1;
         }
         else
         {
            buildT_temp = buildT;
         }

         /* Here minMN, N are reversed */
         ATL_tgerq2(N, minMN,  A+((M-minMN) SHIFT), LDA, TAU, ws_RQ2,
                    ws_T, LDT, WORKM, buildT_temp, 1);
         break; /* END CASE */
      #endif /* defined(ATL_USEPTHREADS) */
   } /* END SWITCH on method */

   /* Common code for cases Serial, PCA COPY, PCA NOCOPY */
   if (M > minMN)
   {
      ATL_larfb(CblasRight, CblasNoTrans,
                LABackward, LARowStore, M-minMN, N, minMN, A+((M -minMN) SHIFT),
                LDA, (ws_T), LDT, A, LDA, WORKM, M);
   }

   return(0);
} /* END ATL_gerqr */


