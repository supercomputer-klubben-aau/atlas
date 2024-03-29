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
#include "atlas_misc.h"
#include "atlas_rblas3.h"
#include "atlas_kernel3.h"
#include "atlas_level3.h"

void ATL_rher2kUC
(
   RC3_HER2K_T                * RTYP,
   ATL_CSZT                   N,
   ATL_CSZT                   K,
   const void                 * ALPHA,
   const void                 * Calph,
   const void                 * A,
   ATL_CSZT                   LDA,
   const void                 * B,
   ATL_CSZT                   LDB,
   const void                 * BETA,
   void                       * C,
   ATL_CSZT                   LDC,
   ATL_CSZT                   RB
)
{
/*
 * Purpose
 * =======
 *
 * ATL_rher2kUC performs the following Hermitian rank-2k operation
 *
 *    C := alpha * conjg( A' ) * B + conjg( alpha ) * conjg( B' ) * A +
 *         beta * C,
 *
 * where  alpha  and  beta are scalars with  beta  real,  C is an n by n
 * (upper) Hermitian matrix and  A  and  B are k by n matrices.
 *
 * This is a type-less recursive version of the algorithm.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   size_t                     size;
   int                        n1, n2;
/* ..
 * .. Executable Statements ..
 *
 */
#ifndef SYR2K_REC
   if( !RTYP->Ther2k( N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) ) return;
#endif
   if( ( n1 = N - RB ) <= 0 )
   {
      ATL_assert( RTYP->Ther2k( N, K, ALPHA, A, LDA, B, LDB, BETA,
                                C, LDC ) == 0 );
      return;
   }

   n2 = N - ( n1 = RB + ( n1 / ( RB << 1 ) ) * RB ); size = RTYP->size;

   ATL_rher2kUC( RTYP, n1, K, ALPHA, Calph, A, LDA, B, LDB, BETA, C, LDC,
                 RB );

   RTYP->Tgemm( n1, n2, K, ALPHA, A, LDA, Mrc3( B, 0, n1, LDB, size ),
                LDB, BETA, Mrc3( C, 0, n1, LDC, size ), LDC );

   RTYP->Tgemm( n1, n2, K, Calph, B, LDB, Mrc3( A, 0, n1, LDA, size ),
                LDA, RTYP->one, Mrc3( C, 0, n1, LDC, size ), LDC );

   ATL_rher2kUC( RTYP, n2, K, ALPHA, Calph, Mrc3( A, 0, n1, LDA, size ),
                 LDA, Mrc3( B, 0, n1, LDB, size ), LDB, BETA, Mrc3( C, n1,
                 n1, LDC, size ), LDC, RB );
/*
 * End of ATL_rher2kUC
 */
}
