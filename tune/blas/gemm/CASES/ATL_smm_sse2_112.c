/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2000 Peter Soendergaard
 */

#ifndef ATL_SSE1
   #error "This routine requires SSE1!"
#endif
#define SSE
#include "SSE3Dnow.h"
#include "atlas_misc.h"

void ATL_USERMM
(const int M, const int N, const int K, const float alpha, const float *A, const int lda, const float *B, const int ldb, const float beta, float *C, const int ldc)

{
   /*--- program info ---*/
   /*  $Revision: 1.3 $  */
   /*  loadfirst = 'a'  */
   /*  nu = 1  */
   /*  k_loop = None  */
   /*  problem = 'gemm'  */
   /*  nregs = 8  */
   /*  split = 1  */
   /*  n_cleanup = {}  */
   /*  mu = 4  */
   /*  ku = 16  */
   /*  rev = '$Revision: 1.3 $'  */
   /*  applyfilter = 0  */
   /*  align_jumps = 0  */
   /*  prec = 'single'  */
   /*  m_cleanup = {}  */
   /*  used_outside_len = 112  */
   /*  k_cleanup = {}  */
   /*  outputdir = 'Linux_P4/'  */
   /*  arch = 'sse'  */
   /*  pipelength = 3  */
   /*  atlasname = 'Linux_P4'  */
   /*  method = 'acc'  */
   /*  used_lastuse = None  */
   /*  used_directload_a = 1  */
   /*  outside_len = 112  */
   /*  veclen = 4  */
   /*  sched = ['spread', 'fuse']  */
   /*  used_directload_b = 0  */
   /*  lastuse = 0  */

   /*--- achitecture specific declarations ---*/

   /*--- program specific declarations ---*/
   int i, j, k;
   vector betavec;
   vector zerovec = {0.0,0.0,0.0,0.0};
   const float *pA0 = A;
   const float *pB0 = B;
   float *pC0 = C;
   const float *stM = A + MB*KB;
   const float *stN = B + NB*KB;
   const int incAm = 4*KB-KB+112;
   const int incBm = -KB+112;
   const int incCm = (4 SHIFT);
   const int incAn = -MB*KB;
   const int incBn = 1*KB;
   const int incCn = ((ldc*1-MB) SHIFT);

   /*--- initial arhitecture specific statements ---*/

   /*--- main program statements ---*/
   vec_mov_mr_1(&beta,reg0);
   vec_mov_rm(reg0,betavec);
   do /* N-loop */
   {
      do /* M-loop */
      {
#ifdef BETA0
         vec_mov_mr(zerovec,reg7);
         vec_mov_rr(reg7,reg0);
         vec_mov_rr(reg7,reg1);
         vec_mov_rr(reg7,reg2);
         vec_mov_rr(reg7,reg3);
#elif defined(BETA1)
         vec_mov_mr_1(pC0,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mov_mr_1(pC0+(2 SHIFT),reg2);
         vec_mov_mr_1(pC0+(3 SHIFT),reg3);
#else
         vec_mov_mr(betavec,reg7);
         vec_mov_mr_1(pC0,reg0);
         vec_mul_rr(reg7,reg0);
         vec_mov_mr_1(pC0+(1 SHIFT),reg1);
         vec_mul_rr(reg7,reg1);
         vec_mov_mr_1(pC0+(2 SHIFT),reg2);
         vec_mul_rr(reg7,reg2);
         vec_mov_mr_1(pC0+(3 SHIFT),reg3);
         vec_mul_rr(reg7,reg3);
#endif
         vec_mov_mr_a(pB0,reg7);
         vec_mov_mr_a(pA0,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+KB,reg5);
         vec_mov_mr_a(pA0+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+4,reg7);
         vec_mov_mr_a(pA0+4,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+4+KB,reg6);
         vec_mov_mr_a(pA0+4+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+4+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+8,reg7);
         vec_mov_mr_a(pA0+8,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+8+KB,reg4);
         vec_mov_mr_a(pA0+8+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+8+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+12,reg7);
         vec_mov_mr_a(pA0+12,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+12+KB,reg5);
         vec_mov_mr_a(pA0+12+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+12+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+16,reg7);
         vec_mov_mr_a(pA0+16,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+16+KB,reg6);
         vec_mov_mr_a(pA0+16+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+16+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+20,reg7);
         vec_mov_mr_a(pA0+20,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+20+KB,reg4);
         vec_mov_mr_a(pA0+20+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+20+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+24,reg7);
         vec_mov_mr_a(pA0+24,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+24+KB,reg5);
         vec_mov_mr_a(pA0+24+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+24+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+28,reg7);
         vec_mov_mr_a(pA0+28,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+28+KB,reg6);
         vec_mov_mr_a(pA0+28+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+28+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+32,reg7);
         vec_mov_mr_a(pA0+32,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+32+KB,reg4);
         vec_mov_mr_a(pA0+32+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+32+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+36,reg7);
         vec_mov_mr_a(pA0+36,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+36+KB,reg5);
         vec_mov_mr_a(pA0+36+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+36+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+40,reg7);
         vec_mov_mr_a(pA0+40,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+40+KB,reg6);
         vec_mov_mr_a(pA0+40+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+40+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+44,reg7);
         vec_mov_mr_a(pA0+44,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+44+KB,reg4);
         vec_mov_mr_a(pA0+44+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+44+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+48,reg7);
         vec_mov_mr_a(pA0+48,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+48+KB,reg5);
         vec_mov_mr_a(pA0+48+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+48+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+52,reg7);
         vec_mov_mr_a(pA0+52,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+52+KB,reg6);
         vec_mov_mr_a(pA0+52+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+52+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+56,reg7);
         vec_mov_mr_a(pA0+56,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+56+KB,reg4);
         vec_mov_mr_a(pA0+56+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+56+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+60,reg7);
         vec_mov_mr_a(pA0+60,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+60+KB,reg5);
         vec_mov_mr_a(pA0+60+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+60+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+64,reg7);
         vec_mov_mr_a(pA0+64,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+64+KB,reg6);
         vec_mov_mr_a(pA0+64+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+64+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+68,reg7);
         vec_mov_mr_a(pA0+68,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+68+KB,reg4);
         vec_mov_mr_a(pA0+68+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+68+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+72,reg7);
         vec_mov_mr_a(pA0+72,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+72+KB,reg5);
         vec_mov_mr_a(pA0+72+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+72+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+76,reg7);
         vec_mov_mr_a(pA0+76,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+76+KB,reg6);
         vec_mov_mr_a(pA0+76+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+76+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+80,reg7);
         vec_mov_mr_a(pA0+80,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+80+KB,reg4);
         vec_mov_mr_a(pA0+80+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+80+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+84,reg7);
         vec_mov_mr_a(pA0+84,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+84+KB,reg5);
         vec_mov_mr_a(pA0+84+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+84+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+88,reg7);
         vec_mov_mr_a(pA0+88,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+88+KB,reg6);
         vec_mov_mr_a(pA0+88+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+88+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+92,reg7);
         vec_mov_mr_a(pA0+92,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+92+KB,reg4);
         vec_mov_mr_a(pA0+92+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+92+3*KB,reg6);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg3);
         vec_mov_mr_a(pB0+96,reg7);
         vec_mov_mr_a(pA0+96,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+96+KB,reg5);
         vec_mov_mr_a(pA0+96+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+96+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+100,reg7);
         vec_mov_mr_a(pA0+100,reg5);
         vec_mul_rr(reg7,reg5);
         vec_mov_mr_a(pA0+100+KB,reg6);
         vec_mov_mr_a(pA0+100+2*KB,reg4);
         vec_add_rr(reg5,reg0);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg1);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg2);
         vec_mov_mr_a(pA0+100+3*KB,reg5);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg3);
         vec_mov_mr_a(pB0+104,reg7);
         vec_mov_mr_a(pA0+104,reg6);
         vec_mul_rr(reg7,reg6);
         vec_mov_mr_a(pA0+104+KB,reg4);
         vec_mov_mr_a(pA0+104+2*KB,reg5);
         vec_add_rr(reg6,reg0);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg1);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg2);
         vec_mov_mr_a(pA0+104+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
         vec_mov_mr_a(pB0+108,reg7);
         vec_mov_mr_a(pA0+108,reg4);
         vec_mul_rr(reg7,reg4);
         vec_mov_mr_a(pA0+108+KB,reg5);
         vec_mov_mr_a(pA0+108+2*KB,reg6);
         vec_add_rr(reg4,reg0);
         vec_mul_rr(reg7,reg5);
         vec_add_rr(reg5,reg1);
         vec_mul_rr(reg7,reg6);
         vec_add_rr(reg6,reg2);
         vec_mov_mr_a(pA0+108+3*KB,reg4);
         vec_mul_rr(reg7,reg4);
         vec_add_rr(reg4,reg3);
#ifndef TCPLX
         vec_sum_full(reg0,reg1,reg2,reg3,reg5,reg6,reg7);
         vec_mov_rm(reg5,pC0);
#else
         vec_sum(reg0);
         vec_sum(reg1);
         vec_sum(reg2);
         vec_sum(reg3);
         vec_mov_rm_1(reg0,pC0);
         vec_mov_rm_1(reg1,pC0+(1 SHIFT));
         vec_mov_rm_1(reg2,pC0+(2 SHIFT));
         vec_mov_rm_1(reg3,pC0+(3 SHIFT));
#endif
         pA0 += incAm;
         pB0 += incBm;
         pC0 += incCm;
      }
      while(pA0 != stM);

      pA0 += incAn;
      pB0 += incBn;
      pC0 += incCn;
   }
   while(pB0 != stN);

   }
