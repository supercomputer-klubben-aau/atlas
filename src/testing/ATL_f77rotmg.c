/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 Antoine P. Petitet
 */
#include "atlas_misc.h"
#include "atlas_tst.h"
#include "atlas_f77blas.h"

void Mjoin( PATL, f77rotmg )
(
   TYPE                      * D1,
   TYPE                      * D2,
   TYPE                      * X1,
   const SCALAR              Y1,
   TYPE                      * PARAM
)
{
   TYPE                      Y10 = Y1;

   F77rotmg( D1, D2, X1, &Y10, PARAM );
}
