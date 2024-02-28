/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 Siju Samuel
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 */

/*
 * This is the C translation of the standard LAPACK Fortran routine:
 *      SUBROUTINE C/Z LACGV( N, X, INCX )
 *
 *  Purpose
 *  =======
 *
 *  ATL_lacgv.c conjugates a complex vector of length N.
 *
 *  Arguments
 *  =========
 *
 *  N       (input) INTEGER
 *          The length of the vector X.  N >= 0.
 *
 *  X       (input/output) complex array, dimension
 *                         (1+(N-1)*abs(INCX))
 *          On entry, the vector of length N to be conjugated.
 *          On exit, X is overwritten with conjg(X).
 *
 *          NOTE : complex numbers are stored as,
 *          real(single/complex), imaginary(single/complex)
 *          in concequtive memory locations.
 *
 *  INCX    (input) INTEGER
 *          The spacing between successive elements of X.
 *
 *  NOTE : rewritten by RCW to just call SCAL of underlying type.
 *
 */
#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_lapack.h"

/* Compiled only to precisions single complex and double complex.             */
void ATL_lacgv(ATL_CINT N, TYPE *X, ATL_CINT INCX)
{
   ATL_CINT incX = (INCX >= 0) ? (INCX+INCX) : ((-INCX)<<1);
   Mjoin(PATLU,scal)(N, ATL_rnone, X+1, incX); /* conj imag elts */
}
