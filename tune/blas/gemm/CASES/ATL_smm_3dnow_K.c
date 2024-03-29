/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2000 Peter Soendergaard
 */

/*****************************************************************************/
/*                             ATL_mm_3dnow_2K.c                             */
/*****************************************************************************/
#if (KB == 2)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                             ATL_mm_3dnow_4K.c                             */
/*****************************************************************************/
#elif (KB == 4)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                             ATL_mm_3dnow_6K.c                             */
/*****************************************************************************/
#elif (KB == 6)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                             ATL_mm_3dnow_8K.c                             */
/*****************************************************************************/
#elif (KB == 8)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_10K.c                             */
/*****************************************************************************/
#elif (KB == 10)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_12K.c                             */
/*****************************************************************************/
#elif (KB == 12)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_14K.c                             */
/*****************************************************************************/
#elif (KB == 14)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_16K.c                             */
/*****************************************************************************/
#elif (KB == 16)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_18K.c                             */
/*****************************************************************************/
#elif (KB == 18)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_20K.c                             */
/*****************************************************************************/
#elif (KB == 20)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_22K.c                             */
/*****************************************************************************/
#elif (KB == 22)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_24K.c                             */
/*****************************************************************************/
#elif (KB == 24)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_26K.c                             */
/*****************************************************************************/
#elif (KB == 26)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_28K.c                             */
/*****************************************************************************/
#elif (KB == 28)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_30K.c                             */
/*****************************************************************************/
#elif (KB == 30)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_32K.c                             */
/*****************************************************************************/
#elif (KB == 32)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_34K.c                             */
/*****************************************************************************/
#elif (KB == 34)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_36K.c                             */
/*****************************************************************************/
#elif (KB == 36)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_38K.c                             */
/*****************************************************************************/
#elif (KB == 38)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_40K.c                             */
/*****************************************************************************/
#elif (KB == 40)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_42K.c                             */
/*****************************************************************************/
#elif (KB == 42)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_44K.c                             */
/*****************************************************************************/
#elif (KB == 44)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_46K.c                             */
/*****************************************************************************/
#elif (KB == 46)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_48K.c                             */
/*****************************************************************************/
#elif (KB == 48)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_50K.c                             */
/*****************************************************************************/
#elif (KB == 50)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_52K.c                             */
/*****************************************************************************/
#elif (KB == 52)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_54K.c                             */
/*****************************************************************************/
#elif (KB == 54)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_56K.c                             */
/*****************************************************************************/
#elif (KB == 56)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_58K.c                             */
/*****************************************************************************/
#elif (KB == 58)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_60K.c                             */
/*****************************************************************************/
#elif (KB == 60)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_62K.c                             */
/*****************************************************************************/
#elif (KB == 62)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_64K.c                             */
/*****************************************************************************/
#elif (KB == 64)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_66K.c                             */
/*****************************************************************************/
#elif (KB == 66)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_68K.c                             */
/*****************************************************************************/
#elif (KB == 68)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_70K.c                             */
/*****************************************************************************/
#elif (KB == 70)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_72K.c                             */
/*****************************************************************************/
#elif (KB == 72)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_74K.c                             */
/*****************************************************************************/
#elif (KB == 74)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_76K.c                             */
/*****************************************************************************/
#elif (KB == 76)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_78K.c                             */
/*****************************************************************************/
#elif (KB == 78)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_80K.c                             */
/*****************************************************************************/
#elif (KB == 80)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_82K.c                             */
/*****************************************************************************/
#elif (KB == 82)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_84K.c                             */
/*****************************************************************************/
#elif (KB == 84)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+4;
   const int incBm = -KB+4;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-4; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_86K.c                             */
/*****************************************************************************/
#elif (KB == 86)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+6;
   const int incBm = -KB+6;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-6; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_88K.c                             */
/*****************************************************************************/
#elif (KB == 88)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+8;
   const int incBm = -KB+8;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-8; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_90K.c                             */
/*****************************************************************************/
#elif (KB == 90)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+10;
   const int incBm = -KB+10;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-10; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_92K.c                             */
/*****************************************************************************/
#elif (KB == 92)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+12;
   const int incBm = -KB+12;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-12; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_94K.c                             */
/*****************************************************************************/
#elif (KB == 94)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+14;
   const int incBm = -KB+14;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-14; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_96K.c                             */
/*****************************************************************************/
#elif (KB == 96)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+16;
   const int incBm = -KB+16;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-16; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+2,reg4);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+2+KB,reg5);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+2,reg6);
         vec_mul_mr(pB0+2,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+2+KB,reg7);
         vec_mul_mr(pB0+2,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+4,reg4);
         vec_mul_mr(pB0+2+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+4+KB,reg5);
         vec_mul_mr(pB0+2+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+4,reg6);
         vec_mul_mr(pB0+4,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+4+KB,reg7);
         vec_mul_mr(pB0+4,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+6,reg4);
         vec_mul_mr(pB0+4+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+6+KB,reg5);
         vec_mul_mr(pB0+4+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+6,reg6);
         vec_mul_mr(pB0+6,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+6+KB,reg7);
         vec_mul_mr(pB0+6,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+8,reg4);
         vec_mul_mr(pB0+6+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+8+KB,reg5);
         vec_mul_mr(pB0+6+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+8,reg6);
         vec_mul_mr(pB0+8,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+8+KB,reg7);
         vec_mul_mr(pB0+8,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+10,reg4);
         vec_mul_mr(pB0+8+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+10+KB,reg5);
         vec_mul_mr(pB0+8+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+10,reg6);
         vec_mul_mr(pB0+10,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+10+KB,reg7);
         vec_mul_mr(pB0+10,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+12,reg4);
         vec_mul_mr(pB0+10+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+12+KB,reg5);
         vec_mul_mr(pB0+10+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+12,reg6);
         vec_mul_mr(pB0+12,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+12+KB,reg7);
         vec_mul_mr(pB0+12,reg5);
         vec_add_rr(reg4,reg0);
         vec_mov_mr(pA0+14,reg4);
         vec_mul_mr(pB0+12+KB,reg6);
         vec_add_rr(reg5,reg1);
         vec_mov_mr(pA0+14+KB,reg5);
         vec_mul_mr(pB0+12+KB,reg7);
         vec_add_rr(reg6,reg2);
         vec_mov_mr(pA0+14,reg6);
         vec_mul_mr(pB0+14,reg4);
         vec_add_rr(reg7,reg3);
         vec_mov_mr(pA0+14+KB,reg7);
         vec_mul_mr(pB0+14,reg5);
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+14+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+14+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

/*****************************************************************************/
/*                            ATL_mm_3dnow_98K.c                             */
/*****************************************************************************/
#elif (KB == 98)
#define THREEDNOW
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   float *pC1 = C+(ldc SHIFT);
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 2*KB-KB+2;
   const int incBm = -KB+2;
   const int incCm = (2 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 2*KB;
   const int incCn = ((ldc*2-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/
   vec_enter();

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_zero(reg0);
         vec_zero(reg1);
         vec_zero(reg2);
         vec_zero(reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC1,reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC1+(1 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr(pA0,reg4);
         vec_mul_mr(pB0,reg4);
         vec_mov_mr(pA0+KB,reg5);
         vec_mul_mr(pB0,reg5);
         vec_mov_mr(pA0,reg6);
         vec_mov_mr(pA0+KB,reg7);
         align();
         for (k=0; k<KB-2; k+=16)
         {
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+2,reg4);
            vec_mul_mr(pB0+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+2+KB,reg5);
            vec_mul_mr(pB0+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+2,reg6);
            vec_mul_mr(pB0+2,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+2+KB,reg7);
            vec_mul_mr(pB0+2,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+4,reg4);
            vec_mul_mr(pB0+2+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+4+KB,reg5);
            vec_mul_mr(pB0+2+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+4,reg6);
            vec_mul_mr(pB0+4,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+4+KB,reg7);
            vec_mul_mr(pB0+4,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+6,reg4);
            vec_mul_mr(pB0+4+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+6+KB,reg5);
            vec_mul_mr(pB0+4+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+6,reg6);
            vec_mul_mr(pB0+6,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+6+KB,reg7);
            vec_mul_mr(pB0+6,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+8,reg4);
            vec_mul_mr(pB0+6+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+8+KB,reg5);
            vec_mul_mr(pB0+6+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+8,reg6);
            vec_mul_mr(pB0+8,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+8+KB,reg7);
            vec_mul_mr(pB0+8,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+10,reg4);
            vec_mul_mr(pB0+8+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+10+KB,reg5);
            vec_mul_mr(pB0+8+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+10,reg6);
            vec_mul_mr(pB0+10,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+10+KB,reg7);
            vec_mul_mr(pB0+10,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+12,reg4);
            vec_mul_mr(pB0+10+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+12+KB,reg5);
            vec_mul_mr(pB0+10+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+12,reg6);
            vec_mul_mr(pB0+12,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+12+KB,reg7);
            vec_mul_mr(pB0+12,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+14,reg4);
            vec_mul_mr(pB0+12+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+14+KB,reg5);
            vec_mul_mr(pB0+12+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+14,reg6);
            vec_mul_mr(pB0+14,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+14+KB,reg7);
            vec_mul_mr(pB0+14,reg5);
            vec_add_rr(reg4,reg0);
            vec_mov_mr(pA0+16,reg4);
            vec_mul_mr(pB0+14+KB,reg6);
            vec_add_rr(reg5,reg1);
            vec_mov_mr(pA0+16+KB,reg5);
            vec_mul_mr(pB0+14+KB,reg7);
            vec_add_rr(reg6,reg2);
            vec_mov_mr(pA0+16,reg6);
            vec_mul_mr(pB0+16,reg4);
            vec_add_rr(reg7,reg3);
            vec_mov_mr(pA0+16+KB,reg7);
            vec_mul_mr(pB0+16,reg5);

            pA0 += 16;
            pB0 += 16;
         }
         vec_add_rr(reg4,reg0);
         vec_add_rr(reg5,reg1);
         vec_mul_mr(pB0+KB,reg6);
         vec_add_rr(reg6,reg2);
         vec_mul_mr(pB0+KB,reg7);
         vec_add_rr(reg7,reg3);
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC1);
         vec_mov_rm_1(reg3,pC1+(1 SHIFT));
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
         pC1 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
      pC1 += incCn;
   }
   while(pB0 != stN);

   vec_exit();
}

#else
   #error Unsupported KB!!
#endif
