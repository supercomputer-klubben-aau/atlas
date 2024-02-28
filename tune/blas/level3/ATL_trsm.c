/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2004 Antoine P. Petitet
 * Code contributers : Antoine P. Petitet, R. Clint Whaley
 */

/*
 * Include files
 */
#include "atlas_rblas3.h"
#include "atlas_kernel3.h"
#include "atlas_lvl3.h"

extern int ATL_trsmtstNB;
#ifdef TRSM_NB
   #undef TRSM_NB
#endif
#define TRSM_NB ATL_trsmtstNB

void Mjoin( PATL, trsm )
(
   const enum ATLAS_SIDE      SIDE,
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const int                  LDA,
   TYPE                       * B,
   const int                  LDB
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, trsm )  solves one of the matrix equations
 *
 *    op( A ) * X = alpha * B,   or  X * op( A ) = alpha * B,
 *
 * where alpha is a scalar, X and B are m by n matrices, A is a unit, or
 * non-unit, upper or lower triangular matrix and op( A ) is one of
 *
 *    op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
 *
 * The matrix X is overwritten on B.
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
   TYPE                       alpha0 = (TYPE)(ALPHA);
   const TYPE                 negone = ATL_rnone, one = ATL_rone;
#else
   TYPE                       negone[2] = { ATL_rnone, ATL_rzero },
                              one   [2] = { ATL_rone,  ATL_rzero };
#endif
   TYPE                       * alpha;
   RC3_FUN_TRSM_T             ATL_rtrsm;
   RC3_TRSM_T                 type;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) ) return;

   if( SCALAR_IS_ZERO( ALPHA ) )
   { Mjoin( PATL, gescal )( M, N, ALPHA, B, LDB ); return; }

#ifdef TREAL
   type.size   = sizeof( TYPE );    type.one = (void *)(&one);
   type.negone = (void *)(&negone); alpha    = &alpha0;
#else
   type.size   = sizeof( TYPE[2] ); type.one = (void *)one;
   type.negone = (void *)negone;    alpha    = (TYPE *)(ALPHA);
#endif

   if( SIDE == AtlasLeft )
   {
      if( TRANS == AtlasNoTrans )
      {
         type.Tgemm = Mjoin( PATL, gemmNN_RB );
         if( UPLO == AtlasUpper )
         {
            ATL_rtrsm = ATL_rtrsmLUN;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLUNN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLUNU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmLLN;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLLNN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLLNU );
         }
      }
#ifdef TREAL
      else
#else
      else if( TRANS == AtlasTrans )
#endif
      {
         type.Tgemm = Mjoin( PATL, gemmTN_RB );
         if( UPLO == AtlasUpper)
         {
            ATL_rtrsm = ATL_rtrsmLUT;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLUTN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLUTU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmLLT;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLLTN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLLTU );
         }
      }
#ifdef TCPLX
      else
      {
         type.Tgemm = Mjoin( PATL, gemmCN_RB );
         if( UPLO == AtlasUpper )
         {
            ATL_rtrsm = ATL_rtrsmLUC;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLUCN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLUCU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmLLC;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmLLCN );
            else                       type.Ttrsm = Mjoin( PATL, trsmLLCU );
         }
      }
#endif
   }
   else
   {
      if( TRANS == AtlasNoTrans )
      {
         type.Tgemm = Mjoin( PATL, gemmNN_RB );
         if( UPLO == AtlasUpper )
         {
            ATL_rtrsm = ATL_rtrsmRUN;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRUNN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRUNU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmRLN;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRLNN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRLNU );
         }
      }
#ifdef TREAL
      else
#else
      else if( TRANS == AtlasTrans )
#endif
      {
         type.Tgemm = Mjoin( PATL, gemmNT_RB );
         if( UPLO == AtlasUpper )
         {
            ATL_rtrsm = ATL_rtrsmRUT;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRUTN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRUTU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmRLT;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRLTN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRLTU );
         }
      }
#ifdef TCPLX
      else
      {
         type.Tgemm = Mjoin( PATL, gemmNC_RB );
         if( UPLO == AtlasUpper )
         {
            ATL_rtrsm = ATL_rtrsmRUC;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRUCN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRUCU );
         }
         else
         {
            ATL_rtrsm = ATL_rtrsmRLC;
            if( DIAG == AtlasNonUnit ) type.Ttrsm = Mjoin( PATL, trsmRLCN );
            else                       type.Ttrsm = Mjoin( PATL, trsmRLCU );
         }
      }
#endif
   }

   ATL_rtrsm( &type, M, N, (void *)(alpha), (void *)(A), LDA, (void *)(B),
              LDB, TRSM_NB );
/*
 * End of Mjoin( PATL, trsm )
 */
}
