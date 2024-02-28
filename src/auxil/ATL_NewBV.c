/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
/*
 * A bitvector consists of an integer array.  The first element of the
 * array indicates the number of bitvecs stored in the array, while any
 * following array entries store 32-bit ints that comprise the full bitvec.
 * Assume 32 bits per int for portability.
 */
ATL_BV_t *ATL_NewBV
(
   unsigned long nmax               /* max # of bits BV needs to store */
)
{
   ATL_BV_t *bv=NULL, nelt;

   if (nmax)
   {
      nelt = (nmax+bpiBV-1) >> shBV;       /* nelt = CEIL(nmax/nbits(ATL_BV_t)) */
      bv = calloc(nelt+1, sizeof(ATL_BV_t)); /* vector begins all unset */
      ATL_assert(bv);
      *bv = nmax;                       /* 1st array elt max # of bits for BV */
   }
   return(bv);
}
