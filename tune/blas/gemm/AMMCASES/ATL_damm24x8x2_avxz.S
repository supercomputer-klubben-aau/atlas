/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2014 R. Clint Whaley
 */
#ifndef ATL_ARCH_XeonPHI
   #define vprefetch0 prefetcht0
   #define vprefetch1 prefetcht1
   #define vprefetch2 prefetcht2
   #define vprefetchwt1 prefetchw
   #define vprefetche0 prefetchw
   #define vprefetche1 prefetchw
   #define vprefetche2 prefetchw
   #define vprefetchw prefetchw
#endif
#include "atlas_asm.h"

#if !defined(KB) || KB == 0
   #ifdef KB
      #undef KB
   #endif
   #if defined(NB) && NB >= 24
      #define KB NB
   #elif defined(MB) && MB >= 24
      #define KB MB
   #else
      #define KB 48
   #endif
#endif
#if defined(SREAL) || defined(SCPLX)
   #define SZ 4
   #define movapd movaps
   #define vmulpd vmulps
   #define vsubpd vsubps
   #define vaddpd vaddps
   #define vfmadd231pd vfmadd231ps
   #define vbroadcastsd vbroadcastss
   #if KB > 48 && defined(ATL_ARCH_XeonPHI)
      #define MOVC vmovnraps
   #else
      #define MOVC vmovaps
   #endif
   #define vbroadcastf64x4 vbroadcastf32x4
#else
   #define SZ 8
   #if KB > 48 && defined(ATL_ARCH_XeonPHI)
      #define MOVC vmovnrapd
   #else
      #define MOVC vmovapd
   #endif
#endif
#define FMAC vfmadd231pd   /* FMAC m256/r256, rs1, rd */
#if defined(BETAN) || defined(BETAn)
   #define BETAN1
#endif
#ifdef BETAN1
   #define VCOP vsubpd
#else
   #define VCOP vaddpd
#endif
#ifdef ATL_ARCH_XeonPHI
   #ifdef BETA0
      #define prefC(m_) vprefetche0 m_
   #else
      #define prefC(m_) vprefetche0 m_
   #endif
   #define prefN(m_) vprefetch2 m_
#else
   #define prefC(m_) vprefetchw m_
   #define prefN(m_) prefetcht2 m_
#endif

#define rB0     %zmm0
#define rB1     %zmm1
#define rA0     %zmm2
#define rA1     %zmm3
#define rA2     %zmm4
#define ra0     %zmm5
#define ra1     %zmm6
#define ra2     %zmm7
#define rC00    %zmm8
#define rC10    %zmm9
#define rC20    %zmm10
#define rC01    %zmm11
#define rC11    %zmm12
#define rC21    %zmm13
#define rC02    %zmm14
#define rC12    %zmm15
#define rC22    %zmm16
#define rC03    %zmm17
#define rC13    %zmm18
#define rC23    %zmm19
#define rC04    %zmm20
#define rC14    %zmm21
#define rC24    %zmm22
#define rC05    %zmm23
#define rC15    %zmm24
#define rC25    %zmm25
#define rC06    %zmm26
#define rC16    %zmm27
#define rC26    %zmm28
#define rC07    %zmm29
#define rC17    %zmm30
#define rC27    %zmm31

/*
 * rbp, rbx, r12 & r15 callee-saved, so use them last
 * 15-4 = 11 regs can be used w/o save
 */

/*
 * ABI-mandated registers: 6 registers --> 5 regs left
 */
#define nmu %rdi
#define nnu %rsi
#define KK  %rdx
#define pA  %rcx
#define pB  %r8
#define pC  %r9
/*
 * regs of choice: 5 regs -->
 */
#define KK0 %r10
#define nnu0 %r11
#define incAm %r13
#define pB0 %rax
#define pA0 %r14
/*
 * Regs that must be saved and restored:rbp,rbx,r12,r15
 */
#define incBn %r15
#define incK  %r12
#define incC  %rbx
#define pfA   %rbp
/*
                    rdi      rsi    rdx        rcx         r8        r9
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                  8(%rsp)    16(%rsp)     24(%rsp)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
.text
ALIGN16
.globl ATL_asmdecor(ATL_USERMM)
ATL_asmdecor(ATL_USERMM):
      vprefetch0 (pB)
   movq %r15, -8(%rsp)
      vprefetch0 64(pB)
   movq %r14, -16(%rsp)
      vprefetch0 128(pB)
   movq %r13, -24(%rsp)
      vprefetch0 192(pB)
   movq %r12, -32(%rsp)
      vprefetch0 256(pB)
   movq %rbp, -40(%rsp)
      vprefetch0 320(pB)
   movq %rbx, -48(%rsp)
      vprefetch0 384(pB)
/*
 * incAm = MU*sizeof*K = 24*8*K = 64*3*K
 */
   lea (KK, KK,2), incAm   /* incAm = 3*K */
      vprefetch0 448(pB)
   shl $6, incAm           /* incAm = 64*3*K */
      vprefetch0 512(pB)
   mov nnu, nnu0
      vprefetch0 576(pB)
   mov pA, pA0
      vprefetch0 640(pB)
   movq 8(%rsp), pfA
/*
 * incBn = NU*SZ*K = 8*SZ*K
 */
   mov KK, incBn
      vprefetch0 704(pB)
   #if defined(SREAL) || defined(SCPLX)
      shl $5, incBn
   #else
      shl $6, incBn
   #endif
      vprefetch0 768(pB)
   mov $1536, incC
      vprefetch0 789(pB)
/*
 * Adjust pB and loop bounds so that we can use loop bound to index pB
 */
   mov incBn, KK
      vprefetch0 (pA)
   add incBn, pB
      vprefetch0 64(pA)
   neg KK
      vprefetch0 128(pA)
   mov pB, pB0
      vprefetch0 3*64(pA)
   mov KK, KK0
      vprefetch0 4*64(pA)
   mov $8*2*SZ, incK  /* incK = NU*KU*SZ */
      vprefetch0 5*64(pA)

#define PFDA0 384
#define PFDA2 8*PFDA0
#define PFDA3 2*PFDA2

#define PFDB0 256
#define PFDB2 8*PFDB0 /* KB*SZ*8 */
#define PFDB3 2*PFDB2
   MLOOP:
         vbroadcastf64x4 (pB,KK), rB0    /* {b3,b2,b1,b0,b3,b2,b1,b0} */
         vbroadcastf64x4 4*SZ(pB,KK), rB1  /* {b7,b6,b5,b4,b7,b6,b5,b4} */
      NLOOP:
         vmovapd (pA), rA0
         vmovapd 64(pA), rA1
         vmovapd 128(pA), rA2
/*
 *       First K-iteration peeled to allow zeroing of rCxx
 */
         vmulpd rB0{aaaa}, rA0, rC00
            vmovapd 192(pA), ra0
         vmulpd rB0{aaaa}, rA1, rC10
            vmovapd 256(pA), ra1
         vmulpd rB0{aaaa}, rA2, rC20
            vmovapd 320(pA), ra2
         vmulpd rB0{bbbb}, rA0, rC01
            vprefetch1 PFDB2(pB,KK)
         vmulpd rB0{bbbb}, rA1, rC11
            vprefetch1 64+PFDB2(pB,KK)
         vmulpd rB0{bbbb}, rA2, rC21
            vprefetch0 PFDB0(pB,KK)
         vmulpd rB0{cccc}, rA0, rC02
            vprefetch0 64+PFDB0(pB,KK)
         vmulpd rB0{cccc}, rA1, rC12
            vprefetch1 PFDA2(pA)
         vmulpd rB0{cccc}, rA2, rC22
            vprefetch1 64+PFDA2(pA)
         vmulpd rB0{dddd}, rA0, rC03
            vprefetch1 2*64+PFDA2(pA)
         vmulpd rB0{dddd}, rA1, rC13
            vprefetch1 3*64+PFDA2(pA)
         vmulpd rB0{dddd}, rA2, rC23
            vbroadcastf64x4 8*SZ(pB,KK), rB0  /* {b3,b2,b1,b0,b3,b2,b1,b0} */
         vmulpd rB1{aaaa}, rA0, rC04
            vprefetch1 4*64+PFDA2(pA)
         vmulpd rB1{aaaa}, rA1, rC14
            vprefetch1 5*64+PFDA2(pA)
         vmulpd rB1{aaaa}, rA2, rC24
            vprefetch0 PFDA0(pA)
         vmulpd rB1{bbbb}, rA0, rC05
            vprefetch0 64+PFDA0(pA)
         vmulpd rB1{bbbb}, rA1, rC15
            vprefetch0 2*64+PFDA0(pA)
         vmulpd rB1{bbbb}, rA2, rC25
            vprefetch0 3*64+PFDA0(pA)
         vmulpd rB1{cccc}, rA0, rC06
            vprefetch0 4*64+PFDA0(pA)
         vmulpd rB1{cccc}, rA1, rC16
            vprefetch0 5*64+PFDA0(pA)
         vmulpd rB1{cccc}, rA2, rC26
         vmulpd rB1{dddd}, rA0, rC07
            vmovapd 384(pA), rA0
         vmulpd rB1{dddd}, rA1, rC17
            vmovapd 448(pA), rA1
         vmulpd rB1{dddd}, rA2, rC27
            vbroadcastf64x4 12*SZ(pB,KK), rB1/* {b11,b10,b9,b8,b11,b10,b9,b8} */
/*
 *       2nd peeled iteration
 */
         FMAC rB0{aaaa}, ra0, rC00
            vmovapd 512(pA), rA2
         FMAC rB0{aaaa}, ra1, rC10
         FMAC rB0{aaaa}, ra2, rC20
         FMAC rB0{bbbb}, ra0, rC01
         FMAC rB0{bbbb}, ra1, rC11
         FMAC rB0{bbbb}, ra2, rC21
         FMAC rB0{cccc}, ra0, rC02
         FMAC rB0{cccc}, ra1, rC12
         FMAC rB0{cccc}, ra2, rC22
         FMAC rB0{dddd}, ra0, rC03
         FMAC rB0{dddd}, ra1, rC13
         FMAC rB0{dddd}, ra2, rC23
            vbroadcastf64x4 16*SZ(pB,KK), rB0  /* {b3,b2,b1,b0,b3,b2,b1,b0} */
         FMAC rB1{aaaa}, ra0, rC04
         FMAC rB1{aaaa}, ra1, rC14
         FMAC rB1{aaaa}, ra2, rC24
         FMAC rB1{bbbb}, ra0, rC05
         FMAC rB1{bbbb}, ra1, rC15
         FMAC rB1{bbbb}, ra2, rC25
         FMAC rB1{cccc}, ra0, rC06
         FMAC rB1{cccc}, ra1, rC16
            add $768, pA
         FMAC rB1{cccc}, ra2, rC26
            add incK, KK
         FMAC rB1{dddd}, ra0, rC07
            vmovapd -192(pA), ra0
         FMAC rB1{dddd}, ra1, rC17
            vmovapd -128(pA), ra1
         FMAC rB1{dddd}, ra2, rC27
            vbroadcastf64x4 4*SZ(pB,KK), rB1  /* {b7,b6,b5,b4,b7,b6,b5,b4} */
         jz KDONE
         KLOOP:
               vmovapd -64(pA), ra2
            FMAC rB0{aaaa}, rA0, rC00
               vprefetch1 PFDB2(pB,KK)
            FMAC rB0{aaaa}, rA1, rC10
               vprefetch1 64+PFDB2(pB,KK)
            FMAC rB0{aaaa}, rA2, rC20
               vprefetch0 PFDB0(pB,KK)
            FMAC rB0{bbbb}, rA0, rC01
               vprefetch0 64+PFDB0(pB,KK)
            FMAC rB0{bbbb}, rA1, rC11
               vprefetch1 PFDA2(pA)
            FMAC rB0{bbbb}, rA2, rC21
               vprefetch1 64+PFDA2(pA)
            FMAC rB0{cccc}, rA0, rC02
               vprefetch1 2*64+PFDA2(pA)
            FMAC rB0{cccc}, rA1, rC12
               vprefetch1 3*64+PFDA2(pA)
            FMAC rB0{cccc}, rA2, rC22
               vprefetch1 4*64+PFDA2(pA)
            FMAC rB0{dddd}, rA0, rC03
               vprefetch1 5*64+PFDA2(pA)
            FMAC rB0{dddd}, rA1, rC13
            FMAC rB0{dddd}, rA2, rC23
               vbroadcastf64x4 8*SZ(pB,KK), rB0  /* {b3,b2,b1,b0,b3,b2,b1,b0} */
            FMAC rB1{aaaa}, rA0, rC04
               vprefetch0 PFDA0(pA)
            FMAC rB1{aaaa}, rA1, rC14
               vprefetch0 64+PFDA0(pA)
            FMAC rB1{aaaa}, rA2, rC24
               vprefetch0 2*64+PFDA0(pA)
            FMAC rB1{bbbb}, rA0, rC05
               vprefetch0 3*64+PFDA0(pA)
            FMAC rB1{bbbb}, rA1, rC15
               vprefetch0 4*64+PFDA0(pA)
            FMAC rB1{bbbb}, rA2, rC25
               vprefetch0 5*64+PFDA0(pA)
            FMAC rB1{cccc}, rA0, rC06
               vprefetch2 PFDA3(pA)
            FMAC rB1{cccc}, rA1, rC16
               vprefetch2 64+PFDA3(pA)
            FMAC rB1{cccc}, rA2, rC26
               vprefetch2 2*64+PFDA3(pA)
            FMAC rB1{dddd}, rA0, rC07
               vmovapd (pA), rA0
            FMAC rB1{dddd}, rA1, rC17
               vmovapd 64(pA), rA1
            FMAC rB1{dddd}, rA2, rC27
               vbroadcastf64x4 12*SZ(pB,KK), rB1/* {b11,b10,b9,b8,b11,b10,b9,b8} */
/*
 *          2nd iteration
 */
            FMAC rB0{aaaa}, ra0, rC00
               vmovapd 128(pA), rA2
            FMAC rB0{aaaa}, ra1, rC10
            FMAC rB0{aaaa}, ra2, rC20
            FMAC rB0{bbbb}, ra0, rC01
            FMAC rB0{bbbb}, ra1, rC11
            FMAC rB0{bbbb}, ra2, rC21
            FMAC rB0{cccc}, ra0, rC02
            FMAC rB0{cccc}, ra1, rC12
            FMAC rB0{cccc}, ra2, rC22
            FMAC rB0{dddd}, ra0, rC03
            FMAC rB0{dddd}, ra1, rC13
            FMAC rB0{dddd}, ra2, rC23
               vbroadcastf64x4 16*SZ(pB,KK), rB0  /* {b3,b2,b1,b0,b3,b2,b1,b0} */
            FMAC rB1{aaaa}, ra0, rC04
            FMAC rB1{aaaa}, ra1, rC14
            FMAC rB1{aaaa}, ra2, rC24
            FMAC rB1{bbbb}, ra0, rC05
            FMAC rB1{bbbb}, ra1, rC15
            FMAC rB1{bbbb}, ra2, rC25
            FMAC rB1{cccc}, ra0, rC06
            FMAC rB1{cccc}, ra1, rC16
               add $384, pA
            FMAC rB1{cccc}, ra2, rC26
               add incK, KK
            FMAC rB1{dddd}, ra0, rC07
               vmovapd 192-384(pA), ra0
            FMAC rB1{dddd}, ra1, rC17
               vmovapd 256-384(pA), ra1
            FMAC rB1{dddd}, ra2, rC27
               vbroadcastf64x4 4*SZ(pB,KK), rB1  /* {b7,b6,b5,b4,b7,b6,b5,b4} */
         jnz KLOOP
         KDONE:
         #ifndef BETA0
            VCOP (pC), rC00, rC00
            VCOP 64(pC), rC10, rC10
            VCOP 128(pC), rC20, rC20
            VCOP 192(pC), rC01, rC01
            VCOP 256(pC), rC11, rC11
            VCOP 320(pC), rC21, rC21
            VCOP 384(pC), rC02, rC02
            VCOP 448(pC), rC12, rC12
         #endif
         MOVC rC00, (pC)
            prefC(1536(pC))
         MOVC rC10, 64(pC)
            prefC(64+1536(pC))
         MOVC rC20, 2*64(pC)
            prefC(2*64+1536(pC))
         MOVC rC01, 3*64(pC)
            prefC(3*64+1536(pC))
         MOVC rC11, 4*64(pC)
            prefC(4*64+1536(pC))
         MOVC rC21, 5*64(pC)
            prefC(5*64+1536(pC))
         MOVC rC02, 6*64(pC)
            prefC(6*64+1536(pC))
         MOVC rC12, 7*64(pC)
            prefC(7*64+1536(pC))
         #ifndef BETA0
            VCOP 512(pC), rC22, rC22
            VCOP 576(pC), rC03, rC03
            VCOP 640(pC), rC13, rC13
            VCOP 704(pC), rC23, rC23
            VCOP 768(pC), rC04, rC04
            VCOP 832(pC), rC14, rC14
            VCOP 896(pC), rC24, rC24
            VCOP 960(pC), rC05, rC05
         #endif
         MOVC rC22, 8*64(pC)
            prefC(8*64+1536(pC))
         MOVC rC03, 9*64(pC)
            prefC(9*64+1536(pC))
         MOVC rC13, 10*64(pC)
            prefC(10*64+1536(pC))
         MOVC rC23, 11*64(pC)
            prefC(11*64+1536(pC))
         MOVC rC04, 12*64(pC)
            prefC(12*64+1536(pC))
         MOVC rC14, 13*64(pC)
            prefC(13*64+1536(pC))
         MOVC rC24, 14*64(pC)
            prefC(14*64+1536(pC))
         MOVC rC05, 15*64(pC)
            prefC(15*64+1536(pC))
         #ifndef BETA0
            VCOP 1024(pC), rC15, rC15
            VCOP 1088(pC), rC25, rC25
            mov KK0, KK
            VCOP 1152(pC), rC06, rC06
            add incBn, pB
            VCOP 1216(pC), rC16, rC16
            mov pA0, pA
            VCOP 1280(pC), rC26, rC26
            VCOP 1344(pC), rC07, rC07
            VCOP 1408(pC), rC17, rC17
            VCOP 1472(pC), rC27, rC27
         #endif
         MOVC rC15, 16*64(pC)
            prefC(16*64+1536(pC))
         MOVC rC25, 17*64(pC)
            prefC(17*64+1536(pC))
         MOVC rC06, 18*64(pC)
            prefC(18*64+1536(pC))
         MOVC rC16, 19*64(pC)
            prefC(19*64+1536(pC))
         MOVC rC26, 20*64(pC)
            prefC(20*64+1536(pC))
         MOVC rC07, 21*64(pC)
            prefC(21*64+1536(pC))
         MOVC rC17, 22*64(pC)
            prefC(22*64+1536(pC))
         MOVC rC27, 23*64(pC)
            prefC(23*64+1536(pC))
         #ifdef BETA0
            add incC, pC
            mov KK0, KK
            add incBn, pB
            mov pA0, pA
         #else
            add incC, pC
         #endif
         dec nnu
      jnz NLOOP
      mov pB0, pB
      add incAm, pA0
      dec nmu
      mov pA0, pA
      mov nnu0, nnu
   jnz MLOOP


DONE:
   movq  -8(%rsp), %r15
   movq -16(%rsp), %r14
   movq -24(%rsp), %r13
   movq -32(%rsp), %r12
   movq -40(%rsp), %rbp
   movq -48(%rsp), %rbx
   ret
