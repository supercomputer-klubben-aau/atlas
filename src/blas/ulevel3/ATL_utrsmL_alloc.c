/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2017 Rakib Hasan
 * Code contributers : Rakib Hasan, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_amm.h"
void* Mjoin(PATL,utrsmL_alloc)
   (sminfo_t*ip, int N, TYPE **Diag, TYPE **L, TYPE **R, TYPE **w)
{
   const int MU = ip->mu, NU=ip->nu;
   int mb = (N + MU - 1) / MU;
   const int MUMU = MU*MU;
   const int MUNU = MU*NU;
   void *vp;
   vp = malloc( ATL_MulBySize( N + (MUMU*mb*(mb+1)/2) + mb*MUNU + 2*MUNU )
                           + 4*ATL_Cachelen);
   if (vp)
   {
      *Diag = ATL_AlignPtr(vp);
      *L = (*Diag) + (N SHIFT);
      *L = ATL_AlignPtr(*L);
      *R = (*L) + ((MUMU*mb*(mb+1)/2) SHIFT);
      *R = ATL_AlignPtr(*R);
      *w = (*R) + ((mb*MUNU) SHIFT);
      *w = ATL_AlignPtr(*w);
   }
   return(vp);
}
