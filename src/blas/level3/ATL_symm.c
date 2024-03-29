/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.9.24 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Contributor(s) : R. Clint Whaley
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
/*
 * Include files
 */
#include "atlas_rblas3.h"
#include "atlas_kernel3.h"
#include "atlas_lvl3.h"

#ifndef SYMM_NB
   #define SYMM_NB 48
#endif

void Mjoin( PATL, symm_APR )
(
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   ATL_CSZT                   M,
   ATL_CSZT                   N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   ATL_CSZT                   LDA,
   const TYPE                 * B,
   ATL_CSZT                   LDB,
   const SCALAR               BETA,
   TYPE                       * C,
   ATL_CSZT                   LDC
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, symm_APR )  performs one of the matrix-matrix operations
 *
 *    C := alpha * A * B + beta * C,
 *
 * or
 *
 *    C := alpha * B * A + beta * C,
 *
 * where alpha and beta are scalars,  A  is a symmetric matrix and B and
 * C are m by n matrices.
 *
 * This is a  recursive  version of the  algorithm.  For a more detailed
 * description of  the arguments of this function, see the reference im-
 * plementation in the  ATLAS/src/blas/reference directory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
#ifdef TREAL
   TYPE                       alpha0 = (TYPE)(ALPHA), beta0 = (TYPE)(BETA);
   const TYPE                 one = ATL_rone;
#else
   TYPE                       one[2] = { ATL_rone, ATL_rzero };
#endif
   TYPE                       * alpha, * beta;
   RC3_FUN_SYMM_T             ATL_rsymm;
   RC3_SYMM_T                 type;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( SCALAR_IS_ZERO( ALPHA ) ) && ( SCALAR_IS_ONE( BETA ) ) ) ) return;

   if( SCALAR_IS_ZERO( ALPHA ) )
   { Mjoin( PATL, gescal )( M, N, BETA, C, LDC ); return; }
#ifdef TREAL
   type.size    = sizeof( TYPE );           type.one = (void *)(&one);
   type.TgemmNN = Mjoin( PATL, gemmNN_RB );
   alpha        = &alpha0;                  beta     = &beta0;
#else
   type.size    = sizeof( TYPE[2] );        type.one = (void *)(one);
   type.TgemmNN = Mjoin( PATL, gemmNN_RB );
   alpha = (TYPE *)(ALPHA);                 beta     = (TYPE *)(BETA);
#endif

   if( SIDE == AtlasLeft )
   {
      type.Tgemm = Mjoin( PATL, gemmTN_RB );
      if( UPLO == AtlasUpper )
      { type.Tsymm = Mjoin( PATL, symmLU ); ATL_rsymm = ATL_rsymmLU; }
      else
      { type.Tsymm = Mjoin( PATL, symmLL ); ATL_rsymm = ATL_rsymmLL; }
   }
   else
   {
      type.Tgemm = Mjoin( PATL, gemmNT_RB );
      if( UPLO == AtlasUpper )
      { type.Tsymm = Mjoin( PATL, symmRU ); ATL_rsymm = ATL_rsymmRU; }
      else
      { type.Tsymm = Mjoin( PATL, symmRL ); ATL_rsymm = ATL_rsymmRL; }
   }

   ATL_rsymm( &type, M, N, ((void *)alpha), ((void *)A), LDA, ((void *)B),
              LDB, ((void *)beta), ((void *)C), LDC, SYMM_NB );
/*
 * End of Mjoin( PATL, symm_APR )
 */
}
