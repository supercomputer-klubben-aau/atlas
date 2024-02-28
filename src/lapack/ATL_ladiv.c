/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 Siju Samuel
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 */

/*
 * This is the C translation of the standard LAPACK Fortran routine:
 *      void      FUNCTION ZLADIV( X, Y, Z)
 *
 * ATL_ladiv.c :
 *             void ATL_cmladiv( TYPE *X, TYPE *Y, TYPE  *Z)
 *     NOTE :  a) ATL_ladiv.c will get compiled to two  precisions
 *                single precision complex,   double precision complex
 *
 *             b) This should be called only for real numbers
 *
 * Purpose
 * =======
 *
 * Z := X / Y, where X and Y are complex.  The computation of X / Y
 * will not overflow on an intermediary step unless the results
 * overflows.
 *
 *        This performs complex division in  real arithmetic
 *
 *                               a + i*b
 *                    p + i*q = ---------
 *                               c + i*d
 *
 *        The algorithm is due to Robert L. Smith and can be found
 *         in D. Knuth, The art of Computer Programming, Vol.2, p.195
 *
 * Arguments
 * =========
 *
 *         X       (input)
 *         Y       (input)
 *                 The complex scalars X and Y ( pointer to X and Y).
 *         Z       (input/output) is the output  ( pointer  to Z)
 *
 */
#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_lapack.h"


void ATL_ladiv(const TYPE *X, const TYPE *Y, TYPE  *Z)
{
   TYPE   E, F;

/* If           X[0], X[1], Y[0], Y[1], &Z[0], &Z[1]  is mapped to            */
/* real numbers  A,    B,    C,    D,    *P,   *Q                             */
/* the computation is as below                                                */
/*                                                                            */
/*   if ( fabs(D) < fabs( C) ) {                                              */
/*         E = D / C ;                                                        */
/*         F = C + D*E ;                                                      */
/*         *P = ( A+B*E ) / F ;                                               */
/*         *Q = ( B-A*E ) / F ;                                               */
/*                                                                            */
/*   } else{                                                                  */
/*         E = C / D ;                                                        */
/*         F = D + C*E ;                                                      */
/*         *P = ( B+A*E ) / F ;                                               */
/*         *Q = ( -A+B*E ) / F ;                                              */
/*   }                                                                        */

   if ( Mabs(Y[1])  < Mabs(Y[0]) )
   {
      E = Y[1]/Y[0];
      F = Y[0] + Y[1]*E;
      *(Z)   = ( X[0] + X[1]*E ) / F ;
      *(Z+1) = ( X[1] - X[0]*E ) / F ;
   }
   else
   {
      E = Y[0]/Y[1];
      F = Y[1] + Y[0]*E;
      *(Z)   = (X[1]+ X[0]*E ) / F ;
      *(Z+1) = (-X[0] + X[1]*E ) / F ;
   }
}                                           /* END AL_ladiv                   */
