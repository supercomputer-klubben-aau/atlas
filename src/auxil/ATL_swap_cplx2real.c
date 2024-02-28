/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
/*
 * Does a swap of complex vector C with two real vectors.  R gets the
 * real components of C, while I gets the imaginary.
 */
void Mjoin(PATL,swap_cplx2real)(ATL_CUINT N, TYPE *C, ATL_CINT incc,
                                TYPE *R, ATL_CINT incR, TYPE *I, ATL_CINT incI)
{
   ATL_CINT incC=incc+incc;
   ATL_UINT i;
   for (i=N; i; i--, C += incC, R += incR, I += incI)
   {
      const register TYPE rv=(*C), iv=C[1];
      *C = *R; C[1] = *I;
      *R = rv; *I = iv;
   }
}
