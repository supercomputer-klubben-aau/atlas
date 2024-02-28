/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"
#ifdef ATL_USEPTHREADS
   #include "atlas_ptalias3.h"
#endif
#include "atlas_level3.h"
#include "cblas.h"

#include <stdarg.h>
int cblas_errprn(int ierr, int info, char *form, ...)
{
   va_list argptr;

   va_start(argptr, form);
#ifdef GCCWIN
   vprintf(form, argptr);
#else
   vfprintf(stderr, form, argptr);
#endif
   va_end(argptr);
   return(Mmin(ierr,info));
}
