/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 Rakib Hasan
 * Code contributers : Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"

#define ATL_ntrcopy ATL_ntrcopyLTUN
#define ATL_ntrsm ATL_ntrsmLT
#define ATL_ntrsm_RR ATL_ntrsmLT_RR
#ifndef Upper_
   void ATL_utrsmR_LT
#else
   void ATL_utrsmR_UN
#endif
   (sminfo_t *ip, const enum ATLAS_DIAG Diag,
    ATL_CINT R, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    TYPE *X, ATL_CINT ldx, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w)
{
   cm2am_t l2a = ip->a2blk, r2a = ip->b2blk;
   ammkern_t amm_b0 = ip->amm_b0, amm_b1 = ip->amm_b1;
   ATL_CSZT ainc = ip->incA;
   const int IsTrans = (ainc != (1 SHIFT));
   const int IsConj = (ip->bv & 8) != 0;
   #ifdef TCPLX
      ammkern_t amm_bn = ip->amm_bn;
      TYPE ONE[3] = {ATL_rone, ATL_rzero, ATL_rzero}, *ZERO=ONE+1;
      TYPE NONE[2] = {ATL_rnone, ATL_rzero};
   #else
      #define ONE ATL_rone
      #define ZERO ATL_rzero
      #define NONE ATL_rnone
   #endif
   const int MU = ip->mu, NU = ip->nu;
   const int NUNU = NU*NU, MUNU = MU*NU;
   const int UnitAlpha = SCALAR_IS_ONE(alpha);
   const int NeedCopyA = (A != NULL);
   TYPE *l, *x=X;
   int mb = (N + NU - 1) / NU;
   const int iRWoff = mb*MUNU;
   ATL_SZT r;
   int i, ix=-1, mbi;

   for (r=0; r < R; r += MU, x += (MU SHIFT))
   {
      int nu, mu = R - r;
      const int DoCopy = (!r && NeedCopyA);
      TYPE *d, *Ac = ((TYPE*)A), *xc = x;
      mu = Mmin(mu, MU);
      nu = Mmin(NU, N);
      /* do the first triangle */
      if (Diag == AtlasUnit) d = NULL;
      else d = diag;
      if (DoCopy)
         ATL_ntrcopy(IsConj, IsTrans, nu, Ac, lda, L, NU, d);
      if (nu == NU)
      {
         if (!UnitAlpha) { ATL_nscal(mu, alpha, xc, ldx); }
         ATL_ntrsm(mu, d, L, NU, xc, ldx);
      }
      else
      {
         if (!UnitAlpha) { ATL_nscal_RR(mu, nu, alpha, xc, ldx); }
         ATL_ntrsm_RR(mu, nu, d, L, NU, xc, ldx);
      }

      for (i=nu, mbi=1, l=L+(NUNU SHIFT);
            i < N; i += NU, mbi++, l+=(NUNU SHIFT))
      {
         nu = Mmin(N-i, NU);
         Ac += (NU SHIFT) * (lda+1);
         xc += (NU SHIFT) * ldx;
         #ifdef TCPLX
         {
            TYPE *rL, *iL, *rR, *iR;
            iR = RW + (mbi-1)*MUNU;
            rR = iR + iRWoff;
            rL = l + mbi*NUNU;
            iL = l;
            r2a(NU, mu, ONE, xc-(NU SHIFT)*ldx, ldx, rR, iR);
            if (DoCopy) /* do the copy as needed for in-cache comp. */
               l2a(i, nu, NONE, Ac-i*ainc, lda, rL, iL);
            iR = RW;
            rR = iR + iRWoff;
            #ifdef USE_TRANS
               amm_b0(1, 1, i, iL, iR, w, rL, iR, w+MUNU);
               amm_b0(1, 1, i, rL, iR, w+MUNU, rL, rR, w);
               amm_bn(1, 1, i, rL, rR, w, iL, rR, w+MUNU);
               amm_b1(1, 1, i, iL, rR, w+MUNU, l+(mbi SHIFT)*NUNU, RW, w);
            #else
               amm_b0(1, 1, i, iR, iL, w, rR, iL, w+MUNU);
               amm_b0(1, 1, i, rR, iL, w+MUNU, rR, rL, w);
               amm_bn(1, 1, i, rR, rL, w, iR, rL, w+MUNU);
               amm_b1(1, 1, i, iR, rL, w+MUNU, RW, l+(mbi SHIFT)*NUNU, w);
            #endif
            ATL_ntrsm_cpC(mu, nu, ONE, w, w+MUNU, alpha, xc, ldx);
         }
         #else
            r2a(NU, mu, ONE, xc-(NU*ldx), ldx, RW+(mbi-1)*MUNU);
            if (DoCopy) /* do the copy as needed for in-cache comp. */
            {
               l2a(i, nu, NONE, Ac-i*ainc, lda, l);
            }
            #ifdef USE_TRANS
               amm_b0(1, 1, i, l, RW, w, l+mbi*NUNU, RW, w);
            #else
               amm_b0(1, 1, i, RW, l, w, RW, l+mbi*NUNU, w);
            #endif
            ATL_ntrsm_cpC(mu, nu, ONE, w, alpha, xc, ldx);
         #endif
         /* now do the diagonal trsm */
         l += (mbi SHIFT) * NUNU;
         if (Diag == AtlasUnit) d = NULL;
         else d = diag + (i SHIFT);
         if (DoCopy)
            ATL_ntrcopy(IsConj, IsTrans, nu, Ac, lda, l, NU, d);
         if (nu == NU)
            ATL_ntrsm(mu, d, l, NU, xc, ldx);
         else
            ATL_ntrsm_RR(mu, nu, d, l, NU, xc, ldx);
      }
   }
}
