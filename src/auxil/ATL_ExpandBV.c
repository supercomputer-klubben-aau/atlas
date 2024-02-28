/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
ATL_BV_t *ATL_ExpandBV(ATL_BV_t *bv, unsigned long newbits)
{
   ATL_BV_t *v;
   unsigned int i, nbit;
   unsigned int neltO, neltN;  /* old and new # ints to store BV */

   if (bv)
      nbit = *bv;
   else
      nbit = 0;
   if (newbits <= nbit)           /* if it already has the required length */
      return(bv);                 /* return same BV */
   neltO = (nbit+bpiBV-1)>>shBV;  /* present # of ints in BV part of array */
   neltN = (newbits+bpiBV)>>shBV; /* how many ints are needed for new nbits? */
/*
 * If increased length doesn't need more ints to store, just change maxbits
 * and return present BV.  New bits were zeroed in original malloc.
 */
   if (neltO == neltN)
   {
      *bv = newbits;
      return(bv);
   }
/*
 * If I reach here, I have to get a longer bitvec, use realloc.
 */
   bv = realloc(bv, sizeof(ATL_BV_t)*(neltN+1));
   *bv = newbits;
/*
 * Any remainder bits in original BV are already 0 by original calloc.
 * realloc does not zero memory, however, so any new elements must now be
 * manually zeroed so that all bits start unset.
 */
   v = bv + 1 + neltO;
   for (i=neltN-neltO; i; i--)
      *v++ = 0;
   return(bv);
}
