/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bitvec.h"
#define ATL_GLOBELT 1
#include "atlas_amm.h"
#undef ATL_GLOBELT
/*
 * Print a matrix in access-major storage
 * flg is bitvec:
 *  0 : if set, w is matrix A, else B or C
 *  1 : if set, w is matrix B, else A or C
 */
void Mjoin(PATL,ipprint)(FILE *fp, ipinfo_t *ip, int flg, char *nm,
                         ATL_iptr_t M, ATL_iptr_t N, const TYPE *w)
{
   ATL_iptr_t i, j;
   #ifdef TCPLX
      void (*getElt)(ipinfo_t*, const TYPE*, ATL_iptr_t,ATL_iptr_t, TYPE*);
   #else
      TYPE (*getElt)(ipinfo_t*, const TYPE*, ATL_iptr_t,ATL_iptr_t);
   #endif
   ATL_assert(flg&3);  /* only printing A & B right now */
   getElt = (flg&1) ? IdxAwElt_ip : IdxBwElt_ip;
   fprintf(fp, "\n%s :\n", nm);
   if (!N)
      return;
   for (i=0; i < M; i++)
   {
      #ifdef TCPLX
         TYPE da[2];
         getElt(ip, w, i, 0, da);
         fprintf(fp, "%+.2e,%+.2e", *da, da[1]);
      #else
         TYPE d;
         d = getElt(ip, w, i, 0);
         fprintf(fp, "%+.2e", d);
      #endif
      for (j=1; j < N; j++)
      {
         #ifdef TCPLX
            getElt(ip, w, i, j, da);
            fprintf(fp, " %+.2e,%+.2e", *da, da[1]);
         #else
            d = getElt(ip, w, i, j);
            fprintf(fp, " %+.2e", d);
         #endif
      }
      fprintf(fp, "\n");
   }
   fprintf(fp, "\n");
}

