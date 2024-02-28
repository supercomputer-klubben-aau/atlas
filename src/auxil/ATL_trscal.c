/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2010 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin( PATL, trscal )
(
   const enum ATLAS_UPLO UPLO,
   ATL_CINT M,
   ATL_CINT N,
   const SCALAR ALPHA,
   TYPE *A,
   ATL_CINT lda
)
{
/*
 * Scales a trapezoidal MxN matrix A by the scalar alpha
 */
#ifndef TCPLX
   register ATL_INT i, j;
   ATL_CINT mn = Mmin(M,N);
   const register TYPE alpha = ALPHA;

   if (alpha == ATL_rone || mn < 1)
      return;
   if( UPLO == AtlasLower )
   {
      if (alpha == ATL_rzero)
      {
         for (j=0; j < mn; j++, A += lda)
            for (i=j; i < M; i++)
               A[i] = alpha;
      }
      else
      {
         for (j=0; j < mn; j++, A += lda)
            for (i=j; i < M; i++)
               A[i] *= alpha;
      }
   }
   else  /* Upper matrix */
   {
      if (alpha == ATL_rzero)
      {
         for (j=0; j < mn; j++, A += lda)
            for (i=0; i <= j; i++)
               A[i] = alpha;
      }
      else
      {
         for (j=0; j < mn; j++, A += lda)
            for (i=0; i <= j; i++)
               A[i] *= alpha;
      }
      if (N > mn)  /* scale rectangular portion */
         Mjoin(PATL,gescal)(M, N-mn, alpha, A, lda);
   }
}
#else
   ATL_CINT M2 = M+M, incA = lda+lda, mn = Mmin(M,N);
   register ATL_INT i, j;
#ifdef TCPLX
   register TYPE ra, ia;
   register const TYPE rb = ALPHA[0], ib = ALPHA[1];
#endif
   if( UPLO == AtlasLower )
   {
      if (ib == ATL_rzero)              /* real scalar */
      {
         if (rb == ATL_rzero)           /* scale by zero */
         {
            for (j=0; j < mn; j++, A += incA)
            {
               for (i=j+j; i < M2; i++)
                  A[i] = rb;
            }
         }
         else if (rb == ATL_rone)       /* no scaling to be done */
            return;
         else                           /* scale by real scalar */
         {
            for (j=0; j < mn; j++, A += incA)
            {
               for (i=j+j; i < M2; i++)
                  A[i] *= rb;
            }
         }
      }
      else /* must apply complex scalar */
      {
         for (j=0; j < mn; j++, A += incA)
         {
            for (i=j+j; i < M2; i += 2)
            {
               ra = A[i];
               ia = A[i+1];
               A[i] = ra*rb - ia*ib;
               A[i+1] = ra*ib + ia*rb;
            }
         }
      }
   }
   else  /* Upper matrix */
   {
      if (ib == ATL_rzero)              /* real scalar */
      {
         if (rb == ATL_rzero)           /* scale by zero */
         {
            for (j=1; j <= mn; j++, A += incA)
            {
               for (i=0; i < j+j; i++)
                  A[i] = ATL_rzero;
            }
         }
         else if (rb == ATL_rone)       /* no scaling to be done */
            return;
         else                           /* scale by real scalar */
         {
            for (j=1; j <= mn; j++, A += incA)
            {
               for (i=0; i < j+j; i++)
                  A[i] *= rb;
            }
         }
      }
      else                              /* scale by complex scalar */
      {
         for (j=1; j <= mn; j++, A += incA)
         {
            for (i=0; i < j+j; i += 2)
            {
               ra = A[i];
               ia = A[i+1];
               A[i] = ra*rb - ia*ib;
               A[i+1] = ra*ib + ia*rb;
            }
         }
      }
/*
 *    Finish off any remaining rectangular portion
 */
      if (N > mn)
         Mjoin(PATL,gescal)(M, N-mn, ALPHA, A, lda);
   }  /* end if over matrix type */
}
#endif
