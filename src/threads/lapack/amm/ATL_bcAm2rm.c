/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from access-major A to cyclic column major */
void Mjoin(PATL, bcAm2rm)
(
   int m,
   int n,
   TYPE *rW,
   ATL_CSZT rbs,
   #ifdef TCPLX
      TYPE *iW,
      ATL_CSZT ibs,
   #endif
   TYPE *A,
   int lda,
   int nb,
   int nt,
   am2cm_t a2rn
)
{
   /* copy back data from A-major to Col-major original storage */
   const int nbnt = nb * nt;
   int m_ = (m/nb)*nb;
   int mr = m - m_;
   TYPE *rWc, *Ac;
   int i;
   #ifdef TREAL
      TYPE none=ATL_rnone;
      if (m <= 0) return;
      for (i=0, rWc=rW, Ac=A; i<m_; i+=nb, rWc+=rbs, Ac+=nbnt)
      {
         a2rn(n, nb, none, Ac, lda, rWc);
      }
      if (mr)
      {
         a2rn(n, mr, none, Ac, lda, rWc);
      }
   #elif defined(TCPLX)
      TYPE *iWc;
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
      if (m <= 0) return;
      for (i=0, rWc=rW, iWc=iW, Ac=A; i<m_;
            i+=nb, rWc+=rbs, iWc+=ibs, Ac+=(nbnt SHIFT))
      {
         a2rn(n, nb, none, Ac, lda, rWc, iWc);
      }
      if (mr)
      {
         a2rn(n, mr, none, Ac, lda, rWc, iWc);
      }
   #endif
}
