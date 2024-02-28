/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from C-access-major to cyclic-column major */
void Mjoin(PATL, bcAblk2cmat)
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
   ablk2cmat_t ca2cc
)
{
   const int nbnt = nb * nt;
   const int m_ = (m/nb)*nb;
   const int mr = m - m_;
   TYPE *Ac, *rWc;
   int i;
   #ifdef TREAL
      TYPE one=ATL_rone, zero=ATL_rzero;
      for (i=0, rWc=rW, Ac=A; i<m_; i+=nb, rWc+=rbs, Ac+=nbnt)
      {
         ca2cc(nb, n, one, rWc, zero, Ac, lda);
      }
      if (mr)
      {
         ca2cc(mr, n, one, rWc, zero, Ac, lda);
      }
   #elif defined(TCPLX)
      TYPE *iWc;
      const TYPE one[2]={ATL_rone,ATL_rzero}, zero[2]={ATL_rzero,ATL_rzero};
      for (i=0, rWc=rW, iWc=iW, Ac=A; i<m_;
            i+=nb, rWc+=rbs, iWc+=ibs, Ac+=(nbnt SHIFT))
      {
         ca2cc(nb, n, one, rWc, iWc, zero, Ac, lda);
      }
      if (mr)
      {
         ca2cc(mr, n, one, rWc, iWc, zero, Ac, lda);
      }
   #endif
}
