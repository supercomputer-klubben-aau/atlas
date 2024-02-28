/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2012, 2009 Siju Samuel
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 */

/*
 * This is the C translation of the standard LAPACK Fortran routine:
 *      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
 *
 *      It is important to use the appropriate sqrt function for the
 *      precision given in order to test against the LAPACK fortran
 *      reference library; in particular due to high serial data
 *      dependence the SGEHRD routine can drift far enough on larger
 *      problems to appear to be in error; about 0.0004 in a single
 *      result element. (Modification by Tony Castaldo, 12/29/2011).
 *
 *  Purpose
 *  =======
 *
 *  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
 *  overflow.
 *
 *  Arguments
 *  =========
 *
 *  X       (input) single/double precision
 *  Y       (input) single/double precision
 *          X and Y specify the values x and y.
 *
 */
#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_lapack.h"
#include "math.h"
/*
 * C90/C89 did not define sqrtf, and so compiling with -ansi will cause
 * it to be unprototyped.  All modern compilers provide it, so just
 * prototype it and use it for FORTRAN/LAPACK compatibility.
 */
#if !defined(__STDC_VERSION__) || __STDC_VERSION__ < 199900
   float sqrtf(float);
#endif

TYPE  ATL_lapy2(TYPE X, TYPE Y)
{
/*
 * Usual case where neither number is a NAN
 */
   if (X == X && Y == Y)
   {
      const TYPE XABS = Mabs(X), YABS = Mabs(Y);
      TYPE  ONE=1.0, ZERO=0.0, W, Z, TEMP;
      if (XABS < YABS)
      {
         W = YABS;
         Z = XABS;
      }
      else
      {
         W = XABS;
         Z = YABS;
      }
      if (Z == ZERO)  /* NOTE: if Z != 0, then W != 0 also, since W >= Z */
         return(W);
      TEMP = Z/W;

      TEMP = ONE + TEMP*TEMP;
      #if defined(SREAL) || defined(SCPLX)
         TEMP = sqrtf(TEMP);
      #else
         TEMP = sqrt(TEMP);
      #endif
      TEMP *= W;
      return(TEMP);
   }
   else if (X == X)  /* Y is a NaN */
      return(Y);     /* return it */
   return(X);        /* X is a NaN, so return it */
}                                               /* END ATL_?lapy2             */


