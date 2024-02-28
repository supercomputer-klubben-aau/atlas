/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from cyclic column major to column major */
void Mjoin(PATL, bcL2G_cpy)
(
   int m,
   int n,
   TYPE* W,
   int ldw,
   TYPE* A,
   int lda,
   int nb,
   int nt
)
{
   const int nbnt = nb * nt;
   int m0 = (m/nb)*nb;
   int mr = m - m0;
   int i;
   TYPE *Wc, *Ac=A;
   if (m <= 0) return;
   Ac += ((m0*nt) SHIFT);
   Wc = W + (m0 SHIFT);
   if (mr)
   {
      Mjoin(PATL, gecopy)(mr, n, Wc, ldw, Ac, lda);
   }
   Ac -= (nbnt SHIFT);
   Wc -= (nb SHIFT);
   for (i=m0; i; i-=nb, Ac-=(nbnt SHIFT), Wc-=(nb SHIFT))
   {
      Mjoin(PATL, gecopy)(nb, n, Wc, ldw, Ac, lda);
   }
}
