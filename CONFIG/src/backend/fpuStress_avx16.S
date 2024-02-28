#define nrep  %rdi
#define dp    %rsi
#ifdef XMM_
   #define zero %xmm0
   #define acc0 %xmm1
   #define acc1 %xmm2
   #define acc2 %xmm3
   #define acc3 %xmm4
   #define acc4 %xmm5
   #define acc5 %xmm6
   #define acc6 %xmm7
   #define acc7 %xmm8
   #define acc8 %xmm9
   #define acc9 %xmm10
   #define accA %xmm11
   #define accB %xmm12
   #define accC %xmm13
   #define accD %xmm14
   #define accE %xmm15
#else
   #define zero %ymm0
   #define acc0 %ymm1
   #define acc1 %ymm2
   #define acc2 %ymm3
   #define acc3 %ymm4
   #define acc4 %ymm5
   #define acc5 %ymm6
   #define acc6 %ymm7
   #define acc7 %ymm8
   #define acc8 %ymm9
   #define acc9 %ymm10
   #define accA %ymm11
   #define accB %ymm12
   #define accC %ymm13
   #define accD %ymm14
   #define accE %ymm15
#endif
/*
 *                               rdi        rsi
 * void fpuStress(unsigned long nrep, double *d);
 */
.text
.globl fpuStress
fpuStress:
   vxorpd zero, zero, zero
   vmovapd zero, accE;
   vmovapd zero, accD;
   vmovapd zero, accC;
   vmovapd zero, accB;
   vmovapd zero, accA;
   vmovapd zero, acc9;
   vmovapd zero, acc8;
   vmovapd zero, acc7;
   vmovapd zero, acc6;
   vmovapd zero, acc5;
   vmovapd zero, acc4;
   vmovapd zero, acc3;
   vmovapd zero, acc2;
   vmovapd zero, acc1;
   vmovapd zero, acc0;
   cmp $0, %rdi
   je RET_LOOP_FLOPS
   .local LOOP
   LOOP:
      vfmadd231pd zero, zero, accE
      vfmadd231pd zero, zero, accD
      vfmadd231pd zero, zero, accC
      vfmadd231pd zero, zero, accB
      vfmadd231pd zero, zero, accA
      vfmadd231pd zero, zero, acc9
      vfmadd231pd zero, zero, acc8
      vfmadd231pd zero, zero, acc7
      vfmadd231pd zero, zero, acc6
      vfmadd231pd zero, zero, acc5
      vfmadd231pd zero, zero, acc4
      vfmadd231pd zero, zero, acc3
      vfmadd231pd zero, zero, acc2
      vfmadd231pd zero, zero, acc1
      vfmadd231pd zero, zero, acc0
      dec nrep
   jnz LOOP
.local RET_LOOP_FLOPS
RET_LOOP_FLOPS:
   fld1             /* ST={1.0} */
   fadd %st,%st     /* ST={2.0} */
   fld1             /* ST={1.0,2.0} */
   fadd  %st(1),%st /* ST={3.0,2.0} */
   fadd %st(1), %st /* ST={5.0,2.0} */
   fmul %st(1), %st /* ST={10.0,2.0} */
   fmul %st, %st(1) /* ST={10.0,20.0} */
   faddp            /* ST={30.0} */
   fadd %st,%st     /* ST={60.0} VLEN=2 */
   #ifndef XMM_
      fadd %st,%st  /* ST={120.0} VLEN=4 */
   #endif
   fstpL (dp)       /* dp[0] = # flops in loop */
   ret
