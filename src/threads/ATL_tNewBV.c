/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_tbitvec.h"
#define UL unsigned long
static UL compSz(UL nbits, unsigned int *P, UL *SUMSZ, UL *LCKSZ, UL *BVSZ,
                 UL *LNBITS, UL *NLRG, UL *NSML)
#undef UL
{
   unsigned long sumSz, lckSz, bvSz, bvElts, lnbits, nlrg, nsml;
   unsigned int p=*P;
   bvSz = nbits >> 2;
   p = *P;
   p = Mmin(p,bvSz);   /* want at least 4 bits per locBV! */
   p = Mmax(p,1);
   lnbits = nbits / p;
   nlrg = nbits - lnbits * p;
   nsml = p - nlrg;
   sumSz = 8*sizeof(long) ;
   sumSz = (unsigned long) ATL_AlignSafeLS(sumSz);
   bvElts = (nlrg) ? lnbits+1 : lnbits;
   bvElts = (bvElts+bpiBV-1)>>shBV;
   bvSz = (bvElts + 2) * sizeof(long);
   bvSz = (unsigned long) ATL_AlignSafeLS(bvSz);
   lckSz = sizeof(ATL_lock_t);
   lckSz = (unsigned long) ATL_AlignSafeLS(lckSz);
   *P = p;
   *SUMSZ = sumSz;
   *LCKSZ = lckSz;
   *BVSZ  = bvSz;
   *LNBITS = lnbits;
   *NLRG = nlrg;
   *NSML = nsml;
   return(bvElts);
}

size_t ATL_tSizeofBV(unsigned long nbits, unsigned int P)
{
   unsigned long sumSz, lckSz, bvSz, lnbits, nlrg, nsml;

   compSz(nbits, &P, &sumSz, &lckSz, &bvSz, &lnbits, &nlrg, &nsml);
   return(sumSz + P*(lckSz + bvSz) + ATL_SAFELS);
}

void ATL_tInitBV(void *vp, unsigned long nbits, unsigned int P)
{
   if (vp)
   {
      unsigned long sumSz, lckSz, bvSz, bvElts, lnbits, nlrg, nsml, i;
      long *lp=ATL_AlignSafeLS(vp);

      bvElts = compSz(nbits, &P, &sumSz, &lckSz, &bvSz, &lnbits, &nlrg, &nsml);
      lp = ATL_AlignSafeLS(vp);
      sumSz = ATL_lDivBySize(sumSz);
      lckSz = ATL_lDivBySize(lckSz);
      bvSz  = ATL_lDivBySize(bvSz);
      lp[ATL_TBV_P] = P;
      lp[ATL_TBV_GNBITS] = nbits;
      lp[ATL_TBV_SUMSZ] = sumSz;
      lp[ATL_TBV_LCKSZ] = lckSz;
      lp[ATL_TBV_BVSZ] = bvSz;
      lp[ATL_TBV_NLGB] = nlrg;
      lp[ATL_TBV_NSMB] = nsml;
      lp[ATL_TBV_B] = lnbits;
/*
 *    Initialize all locks
 */
      lp += sumSz;
      for (i=0; i < P; i++, lp += lckSz)
         ATL_lock_init(lp);
/*
 *    Initialize all locBVs & return
 */
      for (i=0; i < P; i++, lp += bvSz)
      {
         unsigned long k;
         lp[0] = lp[1] = (i < nlrg) ? lnbits+1 : lnbits;
         for (k=0; k < bvElts; k++)
            lp[k+2] = 0;
      }
   }
}

/*
 * 3 areas, each rounded up to 128 bytes:
 * 1. summary area:
 *    <P> <gnbits> <sumSz> <lckSz> <bvSz> <nLrgBlks> <nSmlBlks> <b>
 * 2. local lock storage area : *each lock* rounded up to 128 bytes
 * 3. locBV area: <nunset> [serial BV] -> <nleft> <nbits> <BV> => 128 rounded
 */
void *ATL_tNewBV(unsigned long nbits, unsigned int P)
{
   unsigned long sumSz, lckSz, bvSz, lnbits, nlrg, nsml;
   void *vp;
   long *lp;

   compSz(nbits, &P, &sumSz, &lckSz, &bvSz, &lnbits, &nlrg, &nsml);
/*
 * Allocate memory, and fill in summary area
 */
   vp = malloc(sumSz + P*(lckSz + bvSz) + ATL_SAFELS);
   ATL_assert(vp);
   ATL_tInitBV(vp, nbits, P);
   return(vp);
}
