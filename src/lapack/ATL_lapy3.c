/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2012, 2009 Siju Samuel
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 */

/*-----------------------------------------------------------------------------
 *  This is the C translation of the standard LAPACK Fortran routine:
 *      DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
 *      NOTE : ATL_lapy3.c  will get compiled to
 *             single precision complex (ATL_clapy3.o)  and
 *              double precision complex (ATL_zlapt3.o)
 *
 *   Purpose
 *   =======
 *   ATL_lapy3 returns sqrt(x**2+y**2+z**2), taking care not to cause
 *         unnecessary overflow.
 *
 *   Arguments
 *   =========
 *
 *   X       (input) single/double precision
 *   Y       (input) single/double precision
 *   Z       (input) single/double precision
 *                 X, Y and Z specify the values x, y and z.
-----------------------------------------------------------------------------*/
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

TYPE  ATL_lapy3(const TYPE X, const TYPE Y, const TYPE Z)
{
   TYPE  ONE=1.0, ZERO=0.0, W, Wtemp,  XABS, YABS, ZABS, TEMP;

   XABS = Mabs(X);
   YABS = Mabs(Y);
   ZABS = Mabs(Z);

/* W : get the maximum absolute value from x,y,z                              */
   Wtemp = (XABS<YABS)?YABS:XABS;
   W = (Wtemp<ZABS)?ZABS:Wtemp;

   if (W == ZERO)
   {
/*    W can be zero for max(0,nan,0). Adding all three entries                */
/*    together will make sure  NaN will not disappear.                        */

      return( XABS + YABS + ZABS);
   }
   else
   {
      TEMP =( XABS / W )*( XABS / W ) +
            ( YABS / W )*( YABS / W ) +
            ( ZABS / W )*( ZABS / W ) ;

      #if defined(SREAL) || defined(SCPLX)
         return (W * sqrtf(TEMP));           /* Use single precision sqrt.    */
      #else
         return (W * sqrt(TEMP));            /* Use double precision sqrt.    */
      #endif
   }
}                                            /* END ATL_?lapy3                */
