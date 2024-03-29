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

void ATL_rherkLC
(
   RC3_HERK_T                 * RTYP,
   ATL_CSZT                   N,
   ATL_CSZT                   K,
   const void                 * ALPHA,
   const void                 * A,
   ATL_CSZT                   LDA,
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
 * ATL_rherkLC performs the following Hermitian rank-k operation
 *
 *    C := alpha * conjg( A' ) * A + beta * C,
 *
 * where alpha and beta are real scalars,  C is an n by n (lower) Hermi-
 * tian matrix and  A is an k by n matrix.
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
   if( ( n1 = N - RB ) <= 0 )
   { if(!RTYP->Therk( N, K, ALPHA, A, LDA, BETA, C, LDC)) return; }

   n2 = N - ( n1 = RB + ( n1 / ( RB << 1 ) ) * RB ); size = RTYP->size;

   ATL_rherkLC( RTYP, n1, K, ALPHA, A, LDA, BETA, C, LDC, RB );

   RTYP->Tgemm( n2, n1, K, ALPHA, Mrc3( A, 0, n1, LDA, size ), LDA,
                A, LDA, BETA, Mrc3( C, n1, 0, LDC, size ), LDC );

   ATL_rherkLC( RTYP, n2, K, ALPHA, Mrc3( A, 0, n1, LDA, size ), LDA,
                BETA, Mrc3( C, n1, n1, LDC, size ), LDC, RB );
/*
 * End of ATL_rherkLC
 */
}
