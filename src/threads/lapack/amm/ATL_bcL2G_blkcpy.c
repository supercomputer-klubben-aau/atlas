/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from block-major to cyclic-column major */
void Mjoin(PATL, bcL2G_blkcpy)
(
   int m,
   int n,
   TYPE *W,
   int nb,
   TYPE *A,
   int lda,
   int nt
)
{
   const int nbnb = nb * nb;
   const int nbnt = nb * nt;
   const int m_ = (m/nb)*nb;
   const int mr = m - m_;
   TYPE *Ac, *Wc;
   int i;
   for (i=0, Wc=W, Ac=A; i<m_; i+=nb, Wc+=(nbnb SHIFT), Ac+=(nbnt SHIFT))
   {
      Mjoin(PATL, gecopy)(nb, n, Wc, nb, Ac, lda);
   }
   if (mr)
   {
      Mjoin(PATL, gecopy)(mr, n, Wc, nb, Ac, lda);
   }
}
