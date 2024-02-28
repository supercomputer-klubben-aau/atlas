/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from cyclic-column major to contiguous access-major A */
void Mjoin(PATL, bcRm2am)
(
   int m,
   int n,
   TYPE *A,
   int lda,
   TYPE *rW,
   ATL_CSZT rbs,
   #ifdef TCPLX
      TYPE *iW,
      ATL_CSZT ibs,
   #endif
   int nb,
   int nt,
   cm2am_t r2an
)
{
   const int nbnt = nb * nt;
   const int m_ = (m/nb)*nb;
   const int mr = m - m_;
   TYPE *Ac, *rWc;
   int i;
   #ifdef TREAL
      TYPE none=ATL_rnone;
      for (i=0, rWc=rW, Ac=A; i<m_; i+=nb, rWc+=rbs, Ac+=nbnt)
      {
         r2an(n, nb, none, Ac, lda, rWc);
      }
      if (mr)
      {
         r2an(n, mr, none, Ac, lda, rWc);
      }
   #elif defined(TCPLX)
      TYPE *iWc;
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
      for (i=0, rWc=rW, iWc=iW, Ac=A; i<m_;
            i+=nb, rWc+=rbs, iWc+=ibs, Ac+=(nbnt SHIFT))
      {
         r2an(n, nb, none, Ac, lda, rWc, iWc);
      }
      if (mr)
      {
         r2an(n, mr, none, Ac, lda, rWc, iWc);
      }
   #endif
}
