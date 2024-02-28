#include "atlas_tbitvec.h"
unsigned long ATL_tGetLocalBoundsBV(void *vp, unsigned int rank,
                                    unsigned long *LB)
{
   unsigned long *lp = ATL_AlignSafeLS(vp);
   const unsigned long P=lp[0], nlrg=lp[5], nsml=lp[6], B=lp[7];
   unsigned long off, nblk;

   rank = rank - (rank/P)*P;
   nblk = Mmin(rank, nlrg);
   off = nblk*(B+1);
   nblk -= rank;
   off += nblk*B;
   if (LB)
      *LB = off + ((rank < nlrg) ? B : B-1);
   return(off);
}
