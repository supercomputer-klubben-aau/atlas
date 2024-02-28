/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 * Code contributers : R. Clint Whaley, Antoine P. Petitet
 */

#include "f77wrap_lapack.h"
#include "atlas_lapack.h"
void F77WRAP_ILAENV(const F77_INTEGER *ispec, const F77_INTEGER *irout,
                    const F77_INTEGER *iopts,
                    const F77_INTEGER *N1, const F77_INTEGER *N2,
                    const F77_INTEGER *N3, const F77_INTEGER *N4,
                    F77_INTEGER *iret)
{
   #ifdef ATL_USEPTHREADS
      *iret = ATL_itlaenv(*ispec, *irout, *iopts, *N1, *N2, *N3, *N4);
   #else
      *iret = ATL_ilaenv(*ispec, *irout, *iopts, *N1, *N2, *N3, *N4);
   #endif
}
