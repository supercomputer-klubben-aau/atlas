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

void ATL_rtrsmRUN
(
   RC3_TRSM_T                 * RTYP,
   ATL_CSZT                   M,
   ATL_CSZT                   N,
   const void                 * ALPHA,
   const void                 * A,
   ATL_CSZT                   LDA,
   void                       * B,
   ATL_CSZT                   LDB,
   ATL_CSZT                   RB
)
{
/*
 * Purpose
 * =======
 *
 * ATL_rtrsmRUN solves the matrix equation
 *
 *    X * A = alpha * B,
 *
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit, upper triangular matrix. The matrix X is overwritten on B.
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
   { RTYP->Ttrsm( M, N, ALPHA, A, LDA, B, LDB ); return; }

   n2 = N - ( n1 = RB + ( n1 / ( RB << 1 ) ) * RB ); size = RTYP->size;

   ATL_rtrsmRUN( RTYP, M, n1, ALPHA, A, LDA, B, LDB, RB );

   RTYP->Tgemm( M, n2, n1, RTYP->negone, B, LDB, Mrc3( A, 0, n1, LDA,
                size ), LDA, ALPHA, Mrc3( B, 0, n1, LDB, size ), LDB );

   ATL_rtrsmRUN( RTYP, M, n2, RTYP->one, Mrc3( A, n1, n1, LDA, size ),
                 LDA, Mrc3( B, 0, n1, LDB, size ), LDB, RB );
/*
 * End of ATL_rtrsmRUN
 */
}
