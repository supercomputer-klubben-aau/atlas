/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 Rakib Hasan
 * Code contributers : Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"

#define ATL_ntrcopy ATL_ntrcopyLNUT
#define ATL_ntrsm ATL_ntrsmLN
#define ATL_ntrsm_RR ATL_ntrsmLN_RR
#ifdef Upper_
   void ATL_utrsmL_UT
#else
   void ATL_utrsmL_LN
#endif
   (sminfo_t *ip, const enum ATLAS_DIAG Diag,
    ATL_CINT N, ATL_CINT R, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    TYPE *X, ATL_CINT ldx, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w)
{
   cm2am_t l2a = ip->a2blk;
   cm2am_t r2a = ip->b2blk;
   ammkern_t amm_b0 = ip->amm_b0, amm_b1 = ip->amm_b1;
   ATL_CSZT ainc = ip->incA;
   const int IsTrans = (ainc == (1 SHIFT));
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
   const int MUMU = MU*MU, MUNU = MU*NU;
   const int UnitAlpha = SCALAR_IS_ONE(alpha);
   const int NeedCopyA = (A != NULL);
   TYPE *l, *x=X;
   int mb = (N + MU - 1) / MU;
   const int iRWoff = mb*MUNU;
   ATL_SZT r;
   int i, ix=-1, mbi;

   for (r=0; r < R; r += NU, x += (NU SHIFT)*ldx)
   {
      int mu, nu = R - r;
      const int DoCopy = (!r && NeedCopyA);
      TYPE *d, *Ac = ((TYPE*)A), *xc = x;
      nu = Mmin(nu,NU);
      mu = Mmin(MU, N);

      /* do the first triangle */
      if (Diag == AtlasUnit) d = NULL;
      else d = diag;
      if (DoCopy)
         ATL_ntrcopy(IsConj, IsTrans, mu, Ac, lda, L, MU, d);
      if (mu == MU)
      {
         if (!UnitAlpha) { ATL_nscal(nu, alpha, xc, ldx); }
         ATL_ntrsm(nu, d, L, MU, xc, ldx);
      }
      else
      {
         if (!UnitAlpha) { ATL_nscal_RR(mu, nu, alpha, xc, ldx); }
         ATL_ntrsm_RR(mu, nu, d, L, MU, xc, ldx);
      }

      for (i=mu, mbi=1, l=L+(MUMU SHIFT);
            i < N; i += MU, mbi++, l+=(MUMU SHIFT))
      {
         mu = Mmin(N-i, MU);
         Ac += (MU SHIFT) * (lda+1);
         xc += (MU SHIFT);
         #ifdef TCPLX
         {
            TYPE *rL, *iL, *rR, *iR;
            iR = RW + (mbi-1)*MUNU;
            rR = iR + iRWoff;
            rL = l + mbi*MUMU;
            iL = l;
            r2a(MU, nu, ONE, xc-(MU SHIFT), ldx, rR, iR);
            if (DoCopy) /* do the copy as needed for in-cache comp. */
               l2a(i, mu, NONE, Ac-i*ainc, lda, rL, iL);
            iR = RW;
            rR = iR + iRWoff;
            #ifndef USE_TRANS
               amm_b0(1, 1, i, iL, iR, w, rL, iR, w+MUNU);
               amm_b0(1, 1, i, rL, iR, w+MUNU, rL, rR, w);
               amm_bn(1, 1, i, rL, rR, w, iL, rR, w+MUNU);
               amm_b1(1, 1, i, iL, rR, w+MUNU, l+(mbi SHIFT)*MUMU, RW, w);
            #else
               amm_b0(1, 1, i, iR, iL, w, rR, iL, w+MUNU);
               amm_b0(1, 1, i, rR, iL, w+MUNU, rR, rL, w);
               amm_bn(1, 1, i, rR, rL, w, iR, rL, w+MUNU);
               amm_b1(1, 1, i, iR, rL, w+MUNU, RW, l+(mbi SHIFT)*MUMU, w);
            #endif
            ATL_ntrsm_cpC(mu, nu, ONE, w, w+MUNU, alpha, xc, ldx);
         }
         #else
            r2a(MU, nu, ONE, xc-MU, ldx, RW+(mbi-1)*MUNU);
            if (DoCopy) /* do the copy as needed for in-cache comp. */
               l2a(i, mu, NONE, Ac-i*ainc, lda, l);
            #ifndef USE_TRANS
               amm_b0(1, 1, i, l, RW, w, l+mbi*MUMU, RW, w);
            #else
               amm_b0(1, 1, i, RW, l, w, RW, l+mbi*MUMU, w);
            #endif
            ATL_ntrsm_cpC(mu, nu, ONE, w, alpha, xc, ldx);
         #endif
         /* now do the diagonal trsm */
         l += (mbi SHIFT) * MUMU;
         if (Diag == AtlasUnit) d = NULL;
         else d = diag + (i SHIFT);
         if (DoCopy)
            ATL_ntrcopy(IsConj, IsTrans, mu, Ac, lda, l, MU, d);
         if (mu == MU)
            ATL_ntrsm(nu, d, l, MU, xc, ldx);
         else
            ATL_ntrsm_RR(mu, nu, d, l, MU, xc, ldx);
      }
   }
}
