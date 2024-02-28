/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 Md Rakib Hasan
 * Code contributers : Md Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bcamm.h"

/* copy a panel from column major to cyclic column major */
void Mjoin(PATL, bcG2L_cpy)
(
   int m,
   int n,
   TYPE* A,
   int lda,
   TYPE* W,
   int ldw,
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
      Mjoin(PATL, gecopy)(mr, n, Ac, lda, Wc, ldw);
   }
   Ac -= (nbnt SHIFT);
   Wc -= (nb SHIFT);
   for (i=m0; i; i-=nb, Ac-=(nbnt SHIFT), Wc-=(nb SHIFT))
   {
      Mjoin(PATL, gecopy)(nb, n, Ac, lda, Wc, ldw);
   }
}
