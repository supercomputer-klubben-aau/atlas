/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,gescal)
   (const int M, const int N, const SCALAR beta, TYPE *C, const int ldc)
/*
 * C <- beta*C
 */
{
#ifdef TREAL
   if (beta == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
   else if (beta == ATL_rone) return;
   else Mjoin(PATL,gescal_bX)(M, N, beta, C, ldc);
#else
   TYPE rbeta = *beta;
   if (beta[1] == ATL_rzero)
   {
      if (rbeta == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (rbeta == ATL_rone) return;
      else Mjoin(PATL,gescal_bXi0)(M, N, beta, C, ldc);
   }
   else Mjoin(PATL,gescal_bX)(M, N, beta, C, ldc);
#endif
}
