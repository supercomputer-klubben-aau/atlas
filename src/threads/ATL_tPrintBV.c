/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
void ATL_tPrintBV(FILE *fpout, char *nm, void *bv)
{
   if (bv)
   {
      long *gp=ATL_AlignSafeLS(bv), *tp;
      const unsigned long sumSz=gp[2], lckSz=gp[3], bvSz=gp[4];
      const int P=gp[0];
      unsigned int i;

      fprintf(fpout,
         "%s: P=%u, N=%lu, [sum,lck,bv]Sz=[%lu,%lu,%lu], blk=[%d,%d], b=%lu\n",
         nm?nm:"tbv", P, gp[1], sumSz, lckSz, bvSz, gp[5], gp[6], gp[7]);

      tp = gp + sumSz + P*lckSz;
      for (i=0; i < P; i++)
      {
         const unsigned long N=tp[1], nelt=(N+bpiBV-1)>>shBV;
         unsigned long k;
         fprintf(fpout, "   %u: [N,UN]=[%lu,%lu] bv=%lx", i, N, tp[0],
                 N?tp[2]:0);
         for (k=1; k < nelt; k++)
            fprintf(fpout, ",%lx", tp[2+k]);
         fprintf(fpout, "\n");
      }
   }
   else
      fprintf(fpout, "%s: NULL\n", nm ? nm : "tbv");
}
