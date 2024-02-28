/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_asm.h"
/*
 * This kernel optimized for POWER9, specifically early versions that have
 * generic functional units (sometimes called slices) that do all integer &
 * floating point operations using the same generic units.  This means that
 * doing loop & ptr arithmetic reduces your floating-point peak, since they
 * use the same functional unit!
 * Note that these generic FUs also do all vector permuations, which means
 * we must avoid any type of permute in the loop.  Since m-vectorized kerns
 * require broadcast of B elements (splat for PPC), we use k-vectorized instead
 * which allows us to use straight loads in the innerloop.
 * The system has 4
 * load/store units, so this is fine.
 */

#define rC00 v0
#define rC10 v1
#define rC20 v2
#define rC30 v3
#define rC40 v4
#define rC50 v5
#define rC60 v6
#define rC70 v7
#define rC01 v8
#define rC11 v9
#define rC21 v10
#define rC31 v11
#define rC41 v12
#define rC51 v13
#define rC61 v14
#define rC71 v15
#define rC02 v16
#define rC12 v17
#define rC22 v18
#define rC32 v19
#define rC42 v20
#define rC52 v21
#define rC62 v22
#define rC72 v23
#define rC03 v24
#define rC13 v25
#define rC23 v26
#define rC33 v27
#define rC43 v28
#define rC53 v29
#define rC63 v30
#define rC73 v31
#define rC04 v32
#define rC14 v33
#define rC24 v34
#define rC34 v35
#define rC44 v36
#define rC54 v37
#define rC64 v38
#define rC74 v39
#define rC05 v40
#define rC15 v41
#define rC25 v42
#define rC35 v43
#define rC45 v44
#define rC55 v45
#define rC65 v46
#define rC75 v47

#define rA0 v63
   #define lA0 31
#define rA1 v62
   #define lA1 30
#define rA2 v61
   #define lA2 29
#define rA3 v60
   #define lA3 28
#define rA4 v59
   #define lA4 27
#define rA5 v58
   #define lA5 26
#define rA6 v57
   #define lA6 25
#define rA7 v56
   #define lA7 24

#define rB0 v55
   #define lB0 23
#define rB1 v54
   #define lB1 22
#define rB2 v53
   #define lB2 21
#define rB3 v52
   #define lB3 20
#define rB4 v51
   #define lB4 19
#define rB5 v50
   #define lB5 18

#if defined(SCPLX) && !defined(SREAL)
   #define SREAL 1
#endif
#ifdef SREAL
   #define ELTSZ 4
   #define lxvd2x lxvw4x
   #define xvsubdp xvsubsp
   #define xvadddp xvaddsp
   #define xvmaddadp xvmaddasp
   #define xvmuldp xvmulsp
   #define KMUL *2
   #if KB > 516
   #endif
#else
   #define ELTSZ 8
   #define KMUL
   #if KB > 258
   #endif
#endif
#ifdef BETA1
   #define VADDC xvadddp
#elif !defined(BETA0)
   #define VADDC xvsubdp
#endif
#define VECSZ 16
/*
 * NOTE: 1st param at 96(r1), because caller must reserve space for 8 reg-passed
 * params in param passing space on stack.  ABI says parm passing space at
 * 32(r1) -> 32+64 -> 96(r1).
 * param save space can be written to by callee.
 */
/*
                     r3       r4     r5         r6         r7        r8
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                       r9          r10     96(r1)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
#define pB0 r0
#define NMU r3
#define NNU r4
#define pA  r6
#define pB  r7
#define pC  r8
#define pfA r9
#define pfB r10
#define ADDA addi
#define ADDB addi
#define incAm KB*8*ELTSZ    /* KB*mu*sizeof */
#define incBn KB*6*ELTSZ    /* KB*nu*sizeof */
#define i128_0 0
#define i128_1 r14
#define i128_2 r15
#if defined(__BIG_ENDIAN__)
   #define BIGENDIEN 1
#elif defined(__ORDER_BIG_ENDIAN__) && defined(__BYTE_ORDER__)
   #if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
      #define BIGENDIEN 1
   #endif
#endif
#if defined(SREAL) || defined(SCPLX)
   #ifdef BIGENDIEN
      #define ATL_vvrsum4(s0_, s1_, s2_, s3_, t_, h_) \
         xxmrghw t_, s0_, s1_ ; \
         xxmrglw s0_, s0_, s1_ ; \
         xvaddsp s0_, s0_, t_ ; \
         xxmrghw h_, s2_, s3_ ; \
         xxmrglw s2_, s2_, s3_ ; \
         xvaddsp s2_, s2_, h_ ; \
         xxpermdi t_, s0_, s2_, 0 ; \
         xxpermdi s0_, s0_, s2_, 3 ; \
         xvaddsp s0_, s0_, t_
    #else
      #define ATL_vvrsum4(s3_, s2_, s1_, s0_, t_, h_) \
         xxmrghw t_, s0_, s1_ ; \
         xxmrglw s0_, s0_, s1_ ; \
         xvaddsp s0_, s0_, t_ ; \
         xxmrghw h_, s2_, s3_ ; \
         xxmrglw s2_, s2_, s3_ ; \
         xvaddsp s2_, s2_, h_ ; \
         xxpermdi t_, s0_, s2_, 0 ; \
         xxpermdi s0_, s0_, s2_, 3 ; \
         xvaddsp s3_, s0_, t_
    #endif
#else
   #ifdef BIGENDIEN
      #define ATL_vvrsum2(s0_, s1_, t_) \
         xxpermdi t_, s0_, s1_, 0 ; \
         xxpermdi s0_, s0_, s1_, 3 ; \
         xvadddp s0_, s0_, t_
   #else
      #define ATL_vvrsum2(s1_, s0_, t_) \
         xxpermdi t_, s0_, s1_, 0 ; \
         xxpermdi s0_, s0_, s1_, 3 ; \
         xvadddp s1_, s0_, t_
   #endif
#endif
.local MLOOP
.local NLOOP
.globl ATL_asmdecor(ATL_USERMM)
ATL_asmdecor(ATL_USERMM):
/*
 * Must save & restore:
 *    r1&r2 -> in sframe tail
 *    r14-r31 -> 18*8 = 144 bytes (if used)
 *    f14-f31 -> vsr14-vsr31 -> 18*8 = 144 bytes
 *    v20-v31 -> vsr52-vsr63 -> 12*16 = 192 bytes
 */
#define ISZ 2*8
#define FSZ 18*8
#define VSZ 12*16
#define LSZ 0
#define MSZ 32      /* mandatory 32 bytes for stack frame tail */
#define VOFF (MSZ+LSZ)
#define FOFF (VOFF+VSZ)
#define IOFF (FOFF+FSZ)
#define FSIZE (MSZ+LSZ+VSZ+FSZ+ISZ)
/*
 * Prefetch B using stream 3; pref A using stream 2
 */
   li r9, -1
   rldicr r9, r9, 0, 2
   dcbt r9, pB, 0x8
   li r9, -1
   #if 1
   rldicr r9, r9, 0, 0
   dcbt r9, pA, 0x8
   #endif
/*
 * Allocate local stack frame in order to save & restore
 */
   stdu r1, -FSIZE(r1)   /* allocate stack & store old stack ptr */
/*
 * Save iregs, and initialize fixed 16-increment indices for pC
 */
   std r2, 24(r1)        /* TOC reg always at 24(sp) */
   mfcr r0
   stw r0, 8(r1)         /* CR must be 4-byte saved 8(sp) */
   std r14, IOFF+8*0(r1)
   std r15, IOFF+8*1(r1)
   li  i128_1, 128*1
   li  i128_2, 128*2
/*
 * Save vec regs using AltiVec instructions
 */
   addi pfA, r1, VOFF   /* pfA points to start of vector save area */
   stvx 20, r0, pfA
      addi pfA, pfA, 16
   stvx 21, r0, pfA
      addi pfA, pfA, 16
   stvx 22, r0, pfA
      addi pfA, pfA, 16
   stvx 23, r0, pfA
      addi pfA, pfA, 16
   stvx 24, r0, pfA
      addi pfA, pfA, 16
   stvx 25, r0, pfA
      addi pfA, pfA, 16
   stvx 26, r0, pfA
      addi pfA, pfA, 16
   stvx 27, r0, pfA
      addi pfA, pfA, 16
   stvx 28, r0, pfA
      addi pfA, pfA, 16
   stvx 29, r0, pfA
      addi pfA, pfA, 16
   stvx 30, r0, pfA
      addi pfA, pfA, 16
   stvx 31, r0, pfA
/*
 * f14-31 aliased with VSX 14-31, so they must be saved!
 */
   stfd f14, FOFF+0*8(r1)
   stfd f15, FOFF+1*8(r1)
   stfd f16, FOFF+2*8(r1)
   stfd f17, FOFF+3*8(r1)
   stfd f18, FOFF+4*8(r1)
   stfd f19, FOFF+5*8(r1)
   stfd f20, FOFF+6*8(r1)
   stfd f21, FOFF+7*8(r1)
   stfd f22, FOFF+8*8(r1)
   stfd f23, FOFF+9*8(r1)
   stfd f24, FOFF+10*8(r1)
   stfd f25, FOFF+11*8(r1)
   stfd f26, FOFF+12*8(r1)
   stfd f27, FOFF+13*8(r1)
   stfd f28, FOFF+14*8(r1)
   stfd f29, FOFF+15*8(r1)
   stfd f30, FOFF+16*8(r1)
   stfd f31, FOFF+17*8(r1)
/*
 * PWR addresses take a 16 bit signed offset, so we can double the length
 * of our possible unrolling by adding 2^15-1 to ptrs, and counting offsets
 * from negatives.
 */
#define POFF -32768  /* -2^15 */
   li pB0, POFF
   sub pB, pB, pB0
   sub pA, pA, pB0
   #if defined(pA1) || defined(pB1)
      add pB0, pB0, pB0  /* pB0 = 2*(-32768) */
   #endif
   #if defined(pA1)
      sub pA1, pA, pB0  /* pA1 = pA + 65536 */
   #endif
   #if defined(pB1)
      sub pB1, pB, pB0  /* pB1 = pB + 65536 */
   #endif
   mr pB0, pB
   mtctr NNU

   ALIGN16
   MLOOP:
         lxv rA0, POFF+0(pA)
         lxv rB0, POFF+0(pB)
         lxv rA1, POFF+VECSZ*1(pA)
         lxv rA2, POFF+VECSZ*2(pA)

         lxv rA3, POFF+VECSZ*3(pA)
         lxv rA4, POFF+VECSZ*4(pA)
         lxv rA5, POFF+VECSZ*5(pA)
         lxv rA6, POFF+VECSZ*6(pA)
         lxv rA7, POFF+VECSZ*7(pA)

         lxv rB1, POFF+VECSZ*1(pB)
         lxv rB2, POFF+VECSZ*2(pB)
         lxv rB3, POFF+VECSZ*3(pB)
         lxv rB4, POFF+VECSZ*4(pB)
/*
 *       peeled unroll K=0,1 iteration starts Cregs at zero
 */
         xvmuldp rC00, rA0, rB0
            lxv rB5, POFF+VECSZ*5(pB)
         xvmuldp rC10, rA1, rB0
         xvmuldp rC20, rA2, rB0
            dcbtstt 0, pC
         xvmuldp rC30, rA3, rB0
            dcbtstt i128_1, pC
         xvmuldp rC40, rA4, rB0
         xvmuldp rC50, rA5, rB0
         xvmuldp rC60, rA6, rB0
         xvmuldp rC70, rA7, rB0
         #if KB > 2 KMUL
            lxv rB0, POFF+VECSZ*6(pB)
         #endif
         xvmuldp rC01, rA0, rB1
         xvmuldp rC11, rA1, rB1
         #ifdef SREAL
            dcbtstt i128_2, pC
         #endif
         xvmuldp rC21, rA2, rB1
         xvmuldp rC31, rA3, rB1
         xvmuldp rC41, rA4, rB1
         xvmuldp rC51, rA5, rB1
         xvmuldp rC61, rA6, rB1
         xvmuldp rC71, rA7, rB1
         #if KB > 2 KMUL
            lxv rB1, POFF+VECSZ*7(pB)
         #endif
         xvmuldp rC02, rA0, rB2
         xvmuldp rC12, rA1, rB2
         xvmuldp rC22, rA2, rB2
         xvmuldp rC32, rA3, rB2
         xvmuldp rC42, rA4, rB2
         xvmuldp rC52, rA5, rB2
         xvmuldp rC62, rA6, rB2
         xvmuldp rC72, rA7, rB2
         #if KB > 2 KMUL
            lxv rB2, POFF+VECSZ*8(pB)
         #endif
         xvmuldp rC03, rA0, rB3
         xvmuldp rC13, rA1, rB3
         xvmuldp rC23, rA2, rB3
         xvmuldp rC33, rA3, rB3
         xvmuldp rC43, rA4, rB3
         xvmuldp rC53, rA5, rB3
         xvmuldp rC63, rA6, rB3
         xvmuldp rC73, rA7, rB3
         #if KB > 2 KMUL
            lxv rB3, POFF+VECSZ*9(pB)
         #endif
         xvmuldp rC04, rA0, rB4
         xvmuldp rC14, rA1, rB4
         xvmuldp rC24, rA2, rB4
         xvmuldp rC34, rA3, rB4
         xvmuldp rC44, rA4, rB4
         xvmuldp rC54, rA5, rB4
         xvmuldp rC64, rA6, rB4
         xvmuldp rC74, rA7, rB4
         #if KB > 2 KMUL
            lxv rB4, POFF+VECSZ*10(pB)
         #endif
         xvmuldp rC05, rA0, rB5
         #if KB > 2 KMUL
            lxv rA0, POFF+VECSZ*8(pA)
         #endif
         xvmuldp rC15, rA1, rB5
         #if KB > 2 KMUL
            lxv rA1, POFF+VECSZ*9(pA)
         #endif
         xvmuldp rC25, rA2, rB5
         #if KB > 2 KMUL
            lxv rA2, POFF+VECSZ*10(pA)
         #endif
         xvmuldp rC35, rA3, rB5
         #if KB > 2 KMUL
            lxv rA3, POFF+VECSZ*11(pA)
         #endif
         xvmuldp rC45, rA4, rB5
         #if KB > 2 KMUL
            lxv rA4, POFF+VECSZ*12(pA)
         #endif
         xvmuldp rC55, rA5, rB5
         #if KB > 2 KMUL
            lxv rA5, POFF+VECSZ*13(pA)
         #endif
         xvmuldp rC65, rA6, rB5
         #if KB > 2 KMUL
            lxv rA6, POFF+VECSZ*14(pA)
         #endif
         xvmuldp rC75, rA7, rB5
         #if KB > 2 KMUL
            lxv rA7, POFF+VECSZ*15(pA)
         #endif

/*
 *       K-loop fully unrolled to avoid iops that use VPU/slice
 */
         #if KB >= 4 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*11(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 4 KMUL
               lxv rB0, POFF+VECSZ*12(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 4 KMUL
               lxv rB1, POFF+VECSZ*13(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 4 KMUL
               lxv rB2, POFF+VECSZ*14(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 4 KMUL
               lxv rB3, POFF+VECSZ*15(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 4 KMUL
               lxv rB4, POFF+VECSZ*16(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 4 KMUL
               lxv rA0, POFF+VECSZ*16(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 4 KMUL
               lxv rA1, POFF+VECSZ*17(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 4 KMUL
               lxv rA2, POFF+VECSZ*18(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 4 KMUL
               lxv rA3, POFF+VECSZ*19(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 4 KMUL
               lxv rA4, POFF+VECSZ*20(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 4 KMUL
               lxv rA5, POFF+VECSZ*21(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 4 KMUL
               lxv rA6, POFF+VECSZ*22(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 4 KMUL
               lxv rA7, POFF+VECSZ*23(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 6 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*17(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 6 KMUL
               lxv rB0, POFF+VECSZ*18(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 6 KMUL
               lxv rB1, POFF+VECSZ*19(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 6 KMUL
               lxv rB2, POFF+VECSZ*20(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 6 KMUL
               lxv rB3, POFF+VECSZ*21(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 6 KMUL
               lxv rB4, POFF+VECSZ*22(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 6 KMUL
               lxv rA0, POFF+VECSZ*24(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 6 KMUL
               lxv rA1, POFF+VECSZ*25(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 6 KMUL
               lxv rA2, POFF+VECSZ*26(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 6 KMUL
               lxv rA3, POFF+VECSZ*27(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 6 KMUL
               lxv rA4, POFF+VECSZ*28(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 6 KMUL
               lxv rA5, POFF+VECSZ*29(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 6 KMUL
               lxv rA6, POFF+VECSZ*30(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 6 KMUL
               lxv rA7, POFF+VECSZ*31(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 8 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*23(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 8 KMUL
               lxv rB0, POFF+VECSZ*24(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 8 KMUL
               lxv rB1, POFF+VECSZ*25(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 8 KMUL
               lxv rB2, POFF+VECSZ*26(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 8 KMUL
               lxv rB3, POFF+VECSZ*27(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 8 KMUL
               lxv rB4, POFF+VECSZ*28(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 8 KMUL
               lxv rA0, POFF+VECSZ*32(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 8 KMUL
               lxv rA1, POFF+VECSZ*33(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 8 KMUL
               lxv rA2, POFF+VECSZ*34(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 8 KMUL
               lxv rA3, POFF+VECSZ*35(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 8 KMUL
               lxv rA4, POFF+VECSZ*36(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 8 KMUL
               lxv rA5, POFF+VECSZ*37(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 8 KMUL
               lxv rA6, POFF+VECSZ*38(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 8 KMUL
               lxv rA7, POFF+VECSZ*39(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 10 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*29(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 10 KMUL
               lxv rB0, POFF+VECSZ*30(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 10 KMUL
               lxv rB1, POFF+VECSZ*31(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 10 KMUL
               lxv rB2, POFF+VECSZ*32(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 10 KMUL
               lxv rB3, POFF+VECSZ*33(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 10 KMUL
               lxv rB4, POFF+VECSZ*34(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 10 KMUL
               lxv rA0, POFF+VECSZ*40(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 10 KMUL
               lxv rA1, POFF+VECSZ*41(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 10 KMUL
               lxv rA2, POFF+VECSZ*42(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 10 KMUL
               lxv rA3, POFF+VECSZ*43(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 10 KMUL
               lxv rA4, POFF+VECSZ*44(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 10 KMUL
               lxv rA5, POFF+VECSZ*45(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 10 KMUL
               lxv rA6, POFF+VECSZ*46(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 10 KMUL
               lxv rA7, POFF+VECSZ*47(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 12 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*35(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 12 KMUL
               lxv rB0, POFF+VECSZ*36(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 12 KMUL
               lxv rB1, POFF+VECSZ*37(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 12 KMUL
               lxv rB2, POFF+VECSZ*38(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 12 KMUL
               lxv rB3, POFF+VECSZ*39(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 12 KMUL
               lxv rB4, POFF+VECSZ*40(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 12 KMUL
               lxv rA0, POFF+VECSZ*48(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 12 KMUL
               lxv rA1, POFF+VECSZ*49(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 12 KMUL
               lxv rA2, POFF+VECSZ*50(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 12 KMUL
               lxv rA3, POFF+VECSZ*51(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 12 KMUL
               lxv rA4, POFF+VECSZ*52(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 12 KMUL
               lxv rA5, POFF+VECSZ*53(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 12 KMUL
               lxv rA6, POFF+VECSZ*54(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 12 KMUL
               lxv rA7, POFF+VECSZ*55(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 14 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*41(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 14 KMUL
               lxv rB0, POFF+VECSZ*42(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 14 KMUL
               lxv rB1, POFF+VECSZ*43(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 14 KMUL
               lxv rB2, POFF+VECSZ*44(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 14 KMUL
               lxv rB3, POFF+VECSZ*45(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 14 KMUL
               lxv rB4, POFF+VECSZ*46(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 14 KMUL
               lxv rA0, POFF+VECSZ*56(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 14 KMUL
               lxv rA1, POFF+VECSZ*57(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 14 KMUL
               lxv rA2, POFF+VECSZ*58(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 14 KMUL
               lxv rA3, POFF+VECSZ*59(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 14 KMUL
               lxv rA4, POFF+VECSZ*60(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 14 KMUL
               lxv rA5, POFF+VECSZ*61(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 14 KMUL
               lxv rA6, POFF+VECSZ*62(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 14 KMUL
               lxv rA7, POFF+VECSZ*63(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 16 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*47(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 16 KMUL
               lxv rB0, POFF+VECSZ*48(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 16 KMUL
               lxv rB1, POFF+VECSZ*49(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 16 KMUL
               lxv rB2, POFF+VECSZ*50(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 16 KMUL
               lxv rB3, POFF+VECSZ*51(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 16 KMUL
               lxv rB4, POFF+VECSZ*52(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 16 KMUL
               lxv rA0, POFF+VECSZ*64(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 16 KMUL
               lxv rA1, POFF+VECSZ*65(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 16 KMUL
               lxv rA2, POFF+VECSZ*66(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 16 KMUL
               lxv rA3, POFF+VECSZ*67(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 16 KMUL
               lxv rA4, POFF+VECSZ*68(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 16 KMUL
               lxv rA5, POFF+VECSZ*69(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 16 KMUL
               lxv rA6, POFF+VECSZ*70(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 16 KMUL
               lxv rA7, POFF+VECSZ*71(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 18 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*53(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 18 KMUL
               lxv rB0, POFF+VECSZ*54(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 18 KMUL
               lxv rB1, POFF+VECSZ*55(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 18 KMUL
               lxv rB2, POFF+VECSZ*56(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 18 KMUL
               lxv rB3, POFF+VECSZ*57(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 18 KMUL
               lxv rB4, POFF+VECSZ*58(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 18 KMUL
               lxv rA0, POFF+VECSZ*72(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 18 KMUL
               lxv rA1, POFF+VECSZ*73(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 18 KMUL
               lxv rA2, POFF+VECSZ*74(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 18 KMUL
               lxv rA3, POFF+VECSZ*75(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 18 KMUL
               lxv rA4, POFF+VECSZ*76(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 18 KMUL
               lxv rA5, POFF+VECSZ*77(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 18 KMUL
               lxv rA6, POFF+VECSZ*78(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 18 KMUL
               lxv rA7, POFF+VECSZ*79(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 20 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*59(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 20 KMUL
               lxv rB0, POFF+VECSZ*60(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 20 KMUL
               lxv rB1, POFF+VECSZ*61(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 20 KMUL
               lxv rB2, POFF+VECSZ*62(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 20 KMUL
               lxv rB3, POFF+VECSZ*63(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 20 KMUL
               lxv rB4, POFF+VECSZ*64(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 20 KMUL
               lxv rA0, POFF+VECSZ*80(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 20 KMUL
               lxv rA1, POFF+VECSZ*81(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 20 KMUL
               lxv rA2, POFF+VECSZ*82(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 20 KMUL
               lxv rA3, POFF+VECSZ*83(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 20 KMUL
               lxv rA4, POFF+VECSZ*84(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 20 KMUL
               lxv rA5, POFF+VECSZ*85(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 20 KMUL
               lxv rA6, POFF+VECSZ*86(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 20 KMUL
               lxv rA7, POFF+VECSZ*87(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 22 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*65(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 22 KMUL
               lxv rB0, POFF+VECSZ*66(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 22 KMUL
               lxv rB1, POFF+VECSZ*67(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 22 KMUL
               lxv rB2, POFF+VECSZ*68(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 22 KMUL
               lxv rB3, POFF+VECSZ*69(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 22 KMUL
               lxv rB4, POFF+VECSZ*70(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 22 KMUL
               lxv rA0, POFF+VECSZ*88(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 22 KMUL
               lxv rA1, POFF+VECSZ*89(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 22 KMUL
               lxv rA2, POFF+VECSZ*90(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 22 KMUL
               lxv rA3, POFF+VECSZ*91(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 22 KMUL
               lxv rA4, POFF+VECSZ*92(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 22 KMUL
               lxv rA5, POFF+VECSZ*93(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 22 KMUL
               lxv rA6, POFF+VECSZ*94(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 22 KMUL
               lxv rA7, POFF+VECSZ*95(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 24 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*71(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 24 KMUL
               lxv rB0, POFF+VECSZ*72(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 24 KMUL
               lxv rB1, POFF+VECSZ*73(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 24 KMUL
               lxv rB2, POFF+VECSZ*74(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 24 KMUL
               lxv rB3, POFF+VECSZ*75(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 24 KMUL
               lxv rB4, POFF+VECSZ*76(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 24 KMUL
               lxv rA0, POFF+VECSZ*96(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 24 KMUL
               lxv rA1, POFF+VECSZ*97(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 24 KMUL
               lxv rA2, POFF+VECSZ*98(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 24 KMUL
               lxv rA3, POFF+VECSZ*99(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 24 KMUL
               lxv rA4, POFF+VECSZ*100(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 24 KMUL
               lxv rA5, POFF+VECSZ*101(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 24 KMUL
               lxv rA6, POFF+VECSZ*102(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 24 KMUL
               lxv rA7, POFF+VECSZ*103(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 26 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*77(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 26 KMUL
               lxv rB0, POFF+VECSZ*78(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 26 KMUL
               lxv rB1, POFF+VECSZ*79(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 26 KMUL
               lxv rB2, POFF+VECSZ*80(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 26 KMUL
               lxv rB3, POFF+VECSZ*81(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 26 KMUL
               lxv rB4, POFF+VECSZ*82(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 26 KMUL
               lxv rA0, POFF+VECSZ*104(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 26 KMUL
               lxv rA1, POFF+VECSZ*105(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 26 KMUL
               lxv rA2, POFF+VECSZ*106(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 26 KMUL
               lxv rA3, POFF+VECSZ*107(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 26 KMUL
               lxv rA4, POFF+VECSZ*108(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 26 KMUL
               lxv rA5, POFF+VECSZ*109(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 26 KMUL
               lxv rA6, POFF+VECSZ*110(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 26 KMUL
               lxv rA7, POFF+VECSZ*111(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 28 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*83(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 28 KMUL
               lxv rB0, POFF+VECSZ*84(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 28 KMUL
               lxv rB1, POFF+VECSZ*85(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 28 KMUL
               lxv rB2, POFF+VECSZ*86(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 28 KMUL
               lxv rB3, POFF+VECSZ*87(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 28 KMUL
               lxv rB4, POFF+VECSZ*88(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 28 KMUL
               lxv rA0, POFF+VECSZ*112(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 28 KMUL
               lxv rA1, POFF+VECSZ*113(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 28 KMUL
               lxv rA2, POFF+VECSZ*114(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 28 KMUL
               lxv rA3, POFF+VECSZ*115(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 28 KMUL
               lxv rA4, POFF+VECSZ*116(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 28 KMUL
               lxv rA5, POFF+VECSZ*117(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 28 KMUL
               lxv rA6, POFF+VECSZ*118(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 28 KMUL
               lxv rA7, POFF+VECSZ*119(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 30 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*89(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 30 KMUL
               lxv rB0, POFF+VECSZ*90(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 30 KMUL
               lxv rB1, POFF+VECSZ*91(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 30 KMUL
               lxv rB2, POFF+VECSZ*92(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 30 KMUL
               lxv rB3, POFF+VECSZ*93(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 30 KMUL
               lxv rB4, POFF+VECSZ*94(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 30 KMUL
               lxv rA0, POFF+VECSZ*120(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 30 KMUL
               lxv rA1, POFF+VECSZ*121(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 30 KMUL
               lxv rA2, POFF+VECSZ*122(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 30 KMUL
               lxv rA3, POFF+VECSZ*123(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 30 KMUL
               lxv rA4, POFF+VECSZ*124(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 30 KMUL
               lxv rA5, POFF+VECSZ*125(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 30 KMUL
               lxv rA6, POFF+VECSZ*126(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 30 KMUL
               lxv rA7, POFF+VECSZ*127(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 32 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*95(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 32 KMUL
               lxv rB0, POFF+VECSZ*96(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 32 KMUL
               lxv rB1, POFF+VECSZ*97(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 32 KMUL
               lxv rB2, POFF+VECSZ*98(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 32 KMUL
               lxv rB3, POFF+VECSZ*99(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 32 KMUL
               lxv rB4, POFF+VECSZ*100(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 32 KMUL
               lxv rA0, POFF+VECSZ*128(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 32 KMUL
               lxv rA1, POFF+VECSZ*129(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 32 KMUL
               lxv rA2, POFF+VECSZ*130(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 32 KMUL
               lxv rA3, POFF+VECSZ*131(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 32 KMUL
               lxv rA4, POFF+VECSZ*132(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 32 KMUL
               lxv rA5, POFF+VECSZ*133(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 32 KMUL
               lxv rA6, POFF+VECSZ*134(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 32 KMUL
               lxv rA7, POFF+VECSZ*135(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 34 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*101(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 34 KMUL
               lxv rB0, POFF+VECSZ*102(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 34 KMUL
               lxv rB1, POFF+VECSZ*103(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 34 KMUL
               lxv rB2, POFF+VECSZ*104(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 34 KMUL
               lxv rB3, POFF+VECSZ*105(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 34 KMUL
               lxv rB4, POFF+VECSZ*106(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 34 KMUL
               lxv rA0, POFF+VECSZ*136(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 34 KMUL
               lxv rA1, POFF+VECSZ*137(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 34 KMUL
               lxv rA2, POFF+VECSZ*138(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 34 KMUL
               lxv rA3, POFF+VECSZ*139(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 34 KMUL
               lxv rA4, POFF+VECSZ*140(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 34 KMUL
               lxv rA5, POFF+VECSZ*141(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 34 KMUL
               lxv rA6, POFF+VECSZ*142(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 34 KMUL
               lxv rA7, POFF+VECSZ*143(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 36 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*107(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 36 KMUL
               lxv rB0, POFF+VECSZ*108(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 36 KMUL
               lxv rB1, POFF+VECSZ*109(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 36 KMUL
               lxv rB2, POFF+VECSZ*110(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 36 KMUL
               lxv rB3, POFF+VECSZ*111(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 36 KMUL
               lxv rB4, POFF+VECSZ*112(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 36 KMUL
               lxv rA0, POFF+VECSZ*144(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 36 KMUL
               lxv rA1, POFF+VECSZ*145(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 36 KMUL
               lxv rA2, POFF+VECSZ*146(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 36 KMUL
               lxv rA3, POFF+VECSZ*147(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 36 KMUL
               lxv rA4, POFF+VECSZ*148(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 36 KMUL
               lxv rA5, POFF+VECSZ*149(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 36 KMUL
               lxv rA6, POFF+VECSZ*150(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 36 KMUL
               lxv rA7, POFF+VECSZ*151(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 38 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*113(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 38 KMUL
               lxv rB0, POFF+VECSZ*114(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 38 KMUL
               lxv rB1, POFF+VECSZ*115(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 38 KMUL
               lxv rB2, POFF+VECSZ*116(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 38 KMUL
               lxv rB3, POFF+VECSZ*117(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 38 KMUL
               lxv rB4, POFF+VECSZ*118(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 38 KMUL
               lxv rA0, POFF+VECSZ*152(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 38 KMUL
               lxv rA1, POFF+VECSZ*153(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 38 KMUL
               lxv rA2, POFF+VECSZ*154(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 38 KMUL
               lxv rA3, POFF+VECSZ*155(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 38 KMUL
               lxv rA4, POFF+VECSZ*156(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 38 KMUL
               lxv rA5, POFF+VECSZ*157(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 38 KMUL
               lxv rA6, POFF+VECSZ*158(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 38 KMUL
               lxv rA7, POFF+VECSZ*159(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 40 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*119(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 40 KMUL
               lxv rB0, POFF+VECSZ*120(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 40 KMUL
               lxv rB1, POFF+VECSZ*121(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 40 KMUL
               lxv rB2, POFF+VECSZ*122(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 40 KMUL
               lxv rB3, POFF+VECSZ*123(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 40 KMUL
               lxv rB4, POFF+VECSZ*124(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 40 KMUL
               lxv rA0, POFF+VECSZ*160(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 40 KMUL
               lxv rA1, POFF+VECSZ*161(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 40 KMUL
               lxv rA2, POFF+VECSZ*162(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 40 KMUL
               lxv rA3, POFF+VECSZ*163(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 40 KMUL
               lxv rA4, POFF+VECSZ*164(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 40 KMUL
               lxv rA5, POFF+VECSZ*165(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 40 KMUL
               lxv rA6, POFF+VECSZ*166(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 40 KMUL
               lxv rA7, POFF+VECSZ*167(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 42 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*125(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 42 KMUL
               lxv rB0, POFF+VECSZ*126(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 42 KMUL
               lxv rB1, POFF+VECSZ*127(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 42 KMUL
               lxv rB2, POFF+VECSZ*128(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 42 KMUL
               lxv rB3, POFF+VECSZ*129(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 42 KMUL
               lxv rB4, POFF+VECSZ*130(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 42 KMUL
               lxv rA0, POFF+VECSZ*168(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 42 KMUL
               lxv rA1, POFF+VECSZ*169(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 42 KMUL
               lxv rA2, POFF+VECSZ*170(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 42 KMUL
               lxv rA3, POFF+VECSZ*171(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 42 KMUL
               lxv rA4, POFF+VECSZ*172(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 42 KMUL
               lxv rA5, POFF+VECSZ*173(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 42 KMUL
               lxv rA6, POFF+VECSZ*174(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 42 KMUL
               lxv rA7, POFF+VECSZ*175(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 44 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*131(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 44 KMUL
               lxv rB0, POFF+VECSZ*132(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 44 KMUL
               lxv rB1, POFF+VECSZ*133(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 44 KMUL
               lxv rB2, POFF+VECSZ*134(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 44 KMUL
               lxv rB3, POFF+VECSZ*135(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 44 KMUL
               lxv rB4, POFF+VECSZ*136(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 44 KMUL
               lxv rA0, POFF+VECSZ*176(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 44 KMUL
               lxv rA1, POFF+VECSZ*177(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 44 KMUL
               lxv rA2, POFF+VECSZ*178(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 44 KMUL
               lxv rA3, POFF+VECSZ*179(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 44 KMUL
               lxv rA4, POFF+VECSZ*180(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 44 KMUL
               lxv rA5, POFF+VECSZ*181(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 44 KMUL
               lxv rA6, POFF+VECSZ*182(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 44 KMUL
               lxv rA7, POFF+VECSZ*183(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 46 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*137(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 46 KMUL
               lxv rB0, POFF+VECSZ*138(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 46 KMUL
               lxv rB1, POFF+VECSZ*139(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 46 KMUL
               lxv rB2, POFF+VECSZ*140(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 46 KMUL
               lxv rB3, POFF+VECSZ*141(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 46 KMUL
               lxv rB4, POFF+VECSZ*142(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 46 KMUL
               lxv rA0, POFF+VECSZ*184(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 46 KMUL
               lxv rA1, POFF+VECSZ*185(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 46 KMUL
               lxv rA2, POFF+VECSZ*186(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 46 KMUL
               lxv rA3, POFF+VECSZ*187(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 46 KMUL
               lxv rA4, POFF+VECSZ*188(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 46 KMUL
               lxv rA5, POFF+VECSZ*189(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 46 KMUL
               lxv rA6, POFF+VECSZ*190(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 46 KMUL
               lxv rA7, POFF+VECSZ*191(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 48 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*143(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 48 KMUL
               lxv rB0, POFF+VECSZ*144(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 48 KMUL
               lxv rB1, POFF+VECSZ*145(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 48 KMUL
               lxv rB2, POFF+VECSZ*146(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 48 KMUL
               lxv rB3, POFF+VECSZ*147(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 48 KMUL
               lxv rB4, POFF+VECSZ*148(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 48 KMUL
               lxv rA0, POFF+VECSZ*192(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 48 KMUL
               lxv rA1, POFF+VECSZ*193(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 48 KMUL
               lxv rA2, POFF+VECSZ*194(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 48 KMUL
               lxv rA3, POFF+VECSZ*195(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 48 KMUL
               lxv rA4, POFF+VECSZ*196(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 48 KMUL
               lxv rA5, POFF+VECSZ*197(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 48 KMUL
               lxv rA6, POFF+VECSZ*198(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 48 KMUL
               lxv rA7, POFF+VECSZ*199(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 50 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*149(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 50 KMUL
               lxv rB0, POFF+VECSZ*150(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 50 KMUL
               lxv rB1, POFF+VECSZ*151(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 50 KMUL
               lxv rB2, POFF+VECSZ*152(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 50 KMUL
               lxv rB3, POFF+VECSZ*153(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 50 KMUL
               lxv rB4, POFF+VECSZ*154(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 50 KMUL
               lxv rA0, POFF+VECSZ*200(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 50 KMUL
               lxv rA1, POFF+VECSZ*201(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 50 KMUL
               lxv rA2, POFF+VECSZ*202(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 50 KMUL
               lxv rA3, POFF+VECSZ*203(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 50 KMUL
               lxv rA4, POFF+VECSZ*204(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 50 KMUL
               lxv rA5, POFF+VECSZ*205(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 50 KMUL
               lxv rA6, POFF+VECSZ*206(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 50 KMUL
               lxv rA7, POFF+VECSZ*207(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 52 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*155(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 52 KMUL
               lxv rB0, POFF+VECSZ*156(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 52 KMUL
               lxv rB1, POFF+VECSZ*157(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 52 KMUL
               lxv rB2, POFF+VECSZ*158(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 52 KMUL
               lxv rB3, POFF+VECSZ*159(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 52 KMUL
               lxv rB4, POFF+VECSZ*160(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 52 KMUL
               lxv rA0, POFF+VECSZ*208(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 52 KMUL
               lxv rA1, POFF+VECSZ*209(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 52 KMUL
               lxv rA2, POFF+VECSZ*210(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 52 KMUL
               lxv rA3, POFF+VECSZ*211(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 52 KMUL
               lxv rA4, POFF+VECSZ*212(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 52 KMUL
               lxv rA5, POFF+VECSZ*213(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 52 KMUL
               lxv rA6, POFF+VECSZ*214(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 52 KMUL
               lxv rA7, POFF+VECSZ*215(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 54 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*161(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 54 KMUL
               lxv rB0, POFF+VECSZ*162(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 54 KMUL
               lxv rB1, POFF+VECSZ*163(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 54 KMUL
               lxv rB2, POFF+VECSZ*164(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 54 KMUL
               lxv rB3, POFF+VECSZ*165(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 54 KMUL
               lxv rB4, POFF+VECSZ*166(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 54 KMUL
               lxv rA0, POFF+VECSZ*216(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 54 KMUL
               lxv rA1, POFF+VECSZ*217(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 54 KMUL
               lxv rA2, POFF+VECSZ*218(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 54 KMUL
               lxv rA3, POFF+VECSZ*219(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 54 KMUL
               lxv rA4, POFF+VECSZ*220(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 54 KMUL
               lxv rA5, POFF+VECSZ*221(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 54 KMUL
               lxv rA6, POFF+VECSZ*222(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 54 KMUL
               lxv rA7, POFF+VECSZ*223(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 56 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*167(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 56 KMUL
               lxv rB0, POFF+VECSZ*168(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 56 KMUL
               lxv rB1, POFF+VECSZ*169(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 56 KMUL
               lxv rB2, POFF+VECSZ*170(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 56 KMUL
               lxv rB3, POFF+VECSZ*171(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 56 KMUL
               lxv rB4, POFF+VECSZ*172(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 56 KMUL
               lxv rA0, POFF+VECSZ*224(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 56 KMUL
               lxv rA1, POFF+VECSZ*225(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 56 KMUL
               lxv rA2, POFF+VECSZ*226(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 56 KMUL
               lxv rA3, POFF+VECSZ*227(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 56 KMUL
               lxv rA4, POFF+VECSZ*228(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 56 KMUL
               lxv rA5, POFF+VECSZ*229(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 56 KMUL
               lxv rA6, POFF+VECSZ*230(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 56 KMUL
               lxv rA7, POFF+VECSZ*231(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 58 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*173(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 58 KMUL
               lxv rB0, POFF+VECSZ*174(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 58 KMUL
               lxv rB1, POFF+VECSZ*175(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 58 KMUL
               lxv rB2, POFF+VECSZ*176(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 58 KMUL
               lxv rB3, POFF+VECSZ*177(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 58 KMUL
               lxv rB4, POFF+VECSZ*178(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 58 KMUL
               lxv rA0, POFF+VECSZ*232(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 58 KMUL
               lxv rA1, POFF+VECSZ*233(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 58 KMUL
               lxv rA2, POFF+VECSZ*234(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 58 KMUL
               lxv rA3, POFF+VECSZ*235(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 58 KMUL
               lxv rA4, POFF+VECSZ*236(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 58 KMUL
               lxv rA5, POFF+VECSZ*237(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 58 KMUL
               lxv rA6, POFF+VECSZ*238(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 58 KMUL
               lxv rA7, POFF+VECSZ*239(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 60 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*179(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 60 KMUL
               lxv rB0, POFF+VECSZ*180(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 60 KMUL
               lxv rB1, POFF+VECSZ*181(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 60 KMUL
               lxv rB2, POFF+VECSZ*182(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 60 KMUL
               lxv rB3, POFF+VECSZ*183(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 60 KMUL
               lxv rB4, POFF+VECSZ*184(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 60 KMUL
               lxv rA0, POFF+VECSZ*240(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 60 KMUL
               lxv rA1, POFF+VECSZ*241(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 60 KMUL
               lxv rA2, POFF+VECSZ*242(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 60 KMUL
               lxv rA3, POFF+VECSZ*243(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 60 KMUL
               lxv rA4, POFF+VECSZ*244(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 60 KMUL
               lxv rA5, POFF+VECSZ*245(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 60 KMUL
               lxv rA6, POFF+VECSZ*246(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 60 KMUL
               lxv rA7, POFF+VECSZ*247(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 62 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*185(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 62 KMUL
               lxv rB0, POFF+VECSZ*186(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 62 KMUL
               lxv rB1, POFF+VECSZ*187(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 62 KMUL
               lxv rB2, POFF+VECSZ*188(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 62 KMUL
               lxv rB3, POFF+VECSZ*189(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 62 KMUL
               lxv rB4, POFF+VECSZ*190(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 62 KMUL
               lxv rA0, POFF+VECSZ*248(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 62 KMUL
               lxv rA1, POFF+VECSZ*249(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 62 KMUL
               lxv rA2, POFF+VECSZ*250(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 62 KMUL
               lxv rA3, POFF+VECSZ*251(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 62 KMUL
               lxv rA4, POFF+VECSZ*252(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 62 KMUL
               lxv rA5, POFF+VECSZ*253(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 62 KMUL
               lxv rA6, POFF+VECSZ*254(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 62 KMUL
               lxv rA7, POFF+VECSZ*255(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 64 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*191(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 64 KMUL
               lxv rB0, POFF+VECSZ*192(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 64 KMUL
               lxv rB1, POFF+VECSZ*193(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 64 KMUL
               lxv rB2, POFF+VECSZ*194(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 64 KMUL
               lxv rB3, POFF+VECSZ*195(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 64 KMUL
               lxv rB4, POFF+VECSZ*196(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 64 KMUL
               lxv rA0, POFF+VECSZ*256(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 64 KMUL
               lxv rA1, POFF+VECSZ*257(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 64 KMUL
               lxv rA2, POFF+VECSZ*258(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 64 KMUL
               lxv rA3, POFF+VECSZ*259(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 64 KMUL
               lxv rA4, POFF+VECSZ*260(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 64 KMUL
               lxv rA5, POFF+VECSZ*261(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 64 KMUL
               lxv rA6, POFF+VECSZ*262(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 64 KMUL
               lxv rA7, POFF+VECSZ*263(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 66 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*197(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 66 KMUL
               lxv rB0, POFF+VECSZ*198(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 66 KMUL
               lxv rB1, POFF+VECSZ*199(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 66 KMUL
               lxv rB2, POFF+VECSZ*200(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 66 KMUL
               lxv rB3, POFF+VECSZ*201(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 66 KMUL
               lxv rB4, POFF+VECSZ*202(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 66 KMUL
               lxv rA0, POFF+VECSZ*264(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 66 KMUL
               lxv rA1, POFF+VECSZ*265(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 66 KMUL
               lxv rA2, POFF+VECSZ*266(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 66 KMUL
               lxv rA3, POFF+VECSZ*267(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 66 KMUL
               lxv rA4, POFF+VECSZ*268(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 66 KMUL
               lxv rA5, POFF+VECSZ*269(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 66 KMUL
               lxv rA6, POFF+VECSZ*270(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 66 KMUL
               lxv rA7, POFF+VECSZ*271(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 68 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*203(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 68 KMUL
               lxv rB0, POFF+VECSZ*204(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 68 KMUL
               lxv rB1, POFF+VECSZ*205(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 68 KMUL
               lxv rB2, POFF+VECSZ*206(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 68 KMUL
               lxv rB3, POFF+VECSZ*207(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 68 KMUL
               lxv rB4, POFF+VECSZ*208(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 68 KMUL
               lxv rA0, POFF+VECSZ*272(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 68 KMUL
               lxv rA1, POFF+VECSZ*273(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 68 KMUL
               lxv rA2, POFF+VECSZ*274(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 68 KMUL
               lxv rA3, POFF+VECSZ*275(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 68 KMUL
               lxv rA4, POFF+VECSZ*276(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 68 KMUL
               lxv rA5, POFF+VECSZ*277(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 68 KMUL
               lxv rA6, POFF+VECSZ*278(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 68 KMUL
               lxv rA7, POFF+VECSZ*279(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 70 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*209(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 70 KMUL
               lxv rB0, POFF+VECSZ*210(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 70 KMUL
               lxv rB1, POFF+VECSZ*211(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 70 KMUL
               lxv rB2, POFF+VECSZ*212(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 70 KMUL
               lxv rB3, POFF+VECSZ*213(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 70 KMUL
               lxv rB4, POFF+VECSZ*214(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 70 KMUL
               lxv rA0, POFF+VECSZ*280(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 70 KMUL
               lxv rA1, POFF+VECSZ*281(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 70 KMUL
               lxv rA2, POFF+VECSZ*282(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 70 KMUL
               lxv rA3, POFF+VECSZ*283(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 70 KMUL
               lxv rA4, POFF+VECSZ*284(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 70 KMUL
               lxv rA5, POFF+VECSZ*285(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 70 KMUL
               lxv rA6, POFF+VECSZ*286(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 70 KMUL
               lxv rA7, POFF+VECSZ*287(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 72 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*215(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 72 KMUL
               lxv rB0, POFF+VECSZ*216(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 72 KMUL
               lxv rB1, POFF+VECSZ*217(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 72 KMUL
               lxv rB2, POFF+VECSZ*218(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 72 KMUL
               lxv rB3, POFF+VECSZ*219(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 72 KMUL
               lxv rB4, POFF+VECSZ*220(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 72 KMUL
               lxv rA0, POFF+VECSZ*288(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 72 KMUL
               lxv rA1, POFF+VECSZ*289(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 72 KMUL
               lxv rA2, POFF+VECSZ*290(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 72 KMUL
               lxv rA3, POFF+VECSZ*291(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 72 KMUL
               lxv rA4, POFF+VECSZ*292(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 72 KMUL
               lxv rA5, POFF+VECSZ*293(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 72 KMUL
               lxv rA6, POFF+VECSZ*294(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 72 KMUL
               lxv rA7, POFF+VECSZ*295(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 74 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*221(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 74 KMUL
               lxv rB0, POFF+VECSZ*222(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 74 KMUL
               lxv rB1, POFF+VECSZ*223(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 74 KMUL
               lxv rB2, POFF+VECSZ*224(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 74 KMUL
               lxv rB3, POFF+VECSZ*225(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 74 KMUL
               lxv rB4, POFF+VECSZ*226(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 74 KMUL
               lxv rA0, POFF+VECSZ*296(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 74 KMUL
               lxv rA1, POFF+VECSZ*297(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 74 KMUL
               lxv rA2, POFF+VECSZ*298(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 74 KMUL
               lxv rA3, POFF+VECSZ*299(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 74 KMUL
               lxv rA4, POFF+VECSZ*300(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 74 KMUL
               lxv rA5, POFF+VECSZ*301(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 74 KMUL
               lxv rA6, POFF+VECSZ*302(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 74 KMUL
               lxv rA7, POFF+VECSZ*303(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 76 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*227(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 76 KMUL
               lxv rB0, POFF+VECSZ*228(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 76 KMUL
               lxv rB1, POFF+VECSZ*229(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 76 KMUL
               lxv rB2, POFF+VECSZ*230(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 76 KMUL
               lxv rB3, POFF+VECSZ*231(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 76 KMUL
               lxv rB4, POFF+VECSZ*232(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 76 KMUL
               lxv rA0, POFF+VECSZ*304(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 76 KMUL
               lxv rA1, POFF+VECSZ*305(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 76 KMUL
               lxv rA2, POFF+VECSZ*306(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 76 KMUL
               lxv rA3, POFF+VECSZ*307(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 76 KMUL
               lxv rA4, POFF+VECSZ*308(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 76 KMUL
               lxv rA5, POFF+VECSZ*309(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 76 KMUL
               lxv rA6, POFF+VECSZ*310(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 76 KMUL
               lxv rA7, POFF+VECSZ*311(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 78 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*233(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 78 KMUL
               lxv rB0, POFF+VECSZ*234(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 78 KMUL
               lxv rB1, POFF+VECSZ*235(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 78 KMUL
               lxv rB2, POFF+VECSZ*236(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 78 KMUL
               lxv rB3, POFF+VECSZ*237(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 78 KMUL
               lxv rB4, POFF+VECSZ*238(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 78 KMUL
               lxv rA0, POFF+VECSZ*312(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 78 KMUL
               lxv rA1, POFF+VECSZ*313(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 78 KMUL
               lxv rA2, POFF+VECSZ*314(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 78 KMUL
               lxv rA3, POFF+VECSZ*315(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 78 KMUL
               lxv rA4, POFF+VECSZ*316(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 78 KMUL
               lxv rA5, POFF+VECSZ*317(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 78 KMUL
               lxv rA6, POFF+VECSZ*318(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 78 KMUL
               lxv rA7, POFF+VECSZ*319(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 80 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*239(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 80 KMUL
               lxv rB0, POFF+VECSZ*240(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 80 KMUL
               lxv rB1, POFF+VECSZ*241(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 80 KMUL
               lxv rB2, POFF+VECSZ*242(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 80 KMUL
               lxv rB3, POFF+VECSZ*243(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 80 KMUL
               lxv rB4, POFF+VECSZ*244(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 80 KMUL
               lxv rA0, POFF+VECSZ*320(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 80 KMUL
               lxv rA1, POFF+VECSZ*321(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 80 KMUL
               lxv rA2, POFF+VECSZ*322(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 80 KMUL
               lxv rA3, POFF+VECSZ*323(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 80 KMUL
               lxv rA4, POFF+VECSZ*324(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 80 KMUL
               lxv rA5, POFF+VECSZ*325(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 80 KMUL
               lxv rA6, POFF+VECSZ*326(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 80 KMUL
               lxv rA7, POFF+VECSZ*327(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 82 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*245(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 82 KMUL
               lxv rB0, POFF+VECSZ*246(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 82 KMUL
               lxv rB1, POFF+VECSZ*247(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 82 KMUL
               lxv rB2, POFF+VECSZ*248(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 82 KMUL
               lxv rB3, POFF+VECSZ*249(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 82 KMUL
               lxv rB4, POFF+VECSZ*250(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 82 KMUL
               lxv rA0, POFF+VECSZ*328(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 82 KMUL
               lxv rA1, POFF+VECSZ*329(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 82 KMUL
               lxv rA2, POFF+VECSZ*330(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 82 KMUL
               lxv rA3, POFF+VECSZ*331(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 82 KMUL
               lxv rA4, POFF+VECSZ*332(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 82 KMUL
               lxv rA5, POFF+VECSZ*333(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 82 KMUL
               lxv rA6, POFF+VECSZ*334(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 82 KMUL
               lxv rA7, POFF+VECSZ*335(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 84 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*251(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 84 KMUL
               lxv rB0, POFF+VECSZ*252(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 84 KMUL
               lxv rB1, POFF+VECSZ*253(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 84 KMUL
               lxv rB2, POFF+VECSZ*254(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 84 KMUL
               lxv rB3, POFF+VECSZ*255(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 84 KMUL
               lxv rB4, POFF+VECSZ*256(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 84 KMUL
               lxv rA0, POFF+VECSZ*336(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 84 KMUL
               lxv rA1, POFF+VECSZ*337(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 84 KMUL
               lxv rA2, POFF+VECSZ*338(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 84 KMUL
               lxv rA3, POFF+VECSZ*339(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 84 KMUL
               lxv rA4, POFF+VECSZ*340(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 84 KMUL
               lxv rA5, POFF+VECSZ*341(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 84 KMUL
               lxv rA6, POFF+VECSZ*342(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 84 KMUL
               lxv rA7, POFF+VECSZ*343(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 86 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*257(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 86 KMUL
               lxv rB0, POFF+VECSZ*258(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 86 KMUL
               lxv rB1, POFF+VECSZ*259(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 86 KMUL
               lxv rB2, POFF+VECSZ*260(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 86 KMUL
               lxv rB3, POFF+VECSZ*261(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 86 KMUL
               lxv rB4, POFF+VECSZ*262(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 86 KMUL
               lxv rA0, POFF+VECSZ*344(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 86 KMUL
               lxv rA1, POFF+VECSZ*345(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 86 KMUL
               lxv rA2, POFF+VECSZ*346(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 86 KMUL
               lxv rA3, POFF+VECSZ*347(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 86 KMUL
               lxv rA4, POFF+VECSZ*348(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 86 KMUL
               lxv rA5, POFF+VECSZ*349(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 86 KMUL
               lxv rA6, POFF+VECSZ*350(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 86 KMUL
               lxv rA7, POFF+VECSZ*351(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 88 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*263(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 88 KMUL
               lxv rB0, POFF+VECSZ*264(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 88 KMUL
               lxv rB1, POFF+VECSZ*265(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 88 KMUL
               lxv rB2, POFF+VECSZ*266(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 88 KMUL
               lxv rB3, POFF+VECSZ*267(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 88 KMUL
               lxv rB4, POFF+VECSZ*268(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 88 KMUL
               lxv rA0, POFF+VECSZ*352(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 88 KMUL
               lxv rA1, POFF+VECSZ*353(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 88 KMUL
               lxv rA2, POFF+VECSZ*354(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 88 KMUL
               lxv rA3, POFF+VECSZ*355(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 88 KMUL
               lxv rA4, POFF+VECSZ*356(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 88 KMUL
               lxv rA5, POFF+VECSZ*357(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 88 KMUL
               lxv rA6, POFF+VECSZ*358(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 88 KMUL
               lxv rA7, POFF+VECSZ*359(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 90 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*269(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 90 KMUL
               lxv rB0, POFF+VECSZ*270(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 90 KMUL
               lxv rB1, POFF+VECSZ*271(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 90 KMUL
               lxv rB2, POFF+VECSZ*272(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 90 KMUL
               lxv rB3, POFF+VECSZ*273(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 90 KMUL
               lxv rB4, POFF+VECSZ*274(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 90 KMUL
               lxv rA0, POFF+VECSZ*360(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 90 KMUL
               lxv rA1, POFF+VECSZ*361(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 90 KMUL
               lxv rA2, POFF+VECSZ*362(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 90 KMUL
               lxv rA3, POFF+VECSZ*363(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 90 KMUL
               lxv rA4, POFF+VECSZ*364(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 90 KMUL
               lxv rA5, POFF+VECSZ*365(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 90 KMUL
               lxv rA6, POFF+VECSZ*366(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 90 KMUL
               lxv rA7, POFF+VECSZ*367(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 92 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*275(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 92 KMUL
               lxv rB0, POFF+VECSZ*276(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 92 KMUL
               lxv rB1, POFF+VECSZ*277(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 92 KMUL
               lxv rB2, POFF+VECSZ*278(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 92 KMUL
               lxv rB3, POFF+VECSZ*279(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 92 KMUL
               lxv rB4, POFF+VECSZ*280(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 92 KMUL
               lxv rA0, POFF+VECSZ*368(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 92 KMUL
               lxv rA1, POFF+VECSZ*369(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 92 KMUL
               lxv rA2, POFF+VECSZ*370(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 92 KMUL
               lxv rA3, POFF+VECSZ*371(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 92 KMUL
               lxv rA4, POFF+VECSZ*372(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 92 KMUL
               lxv rA5, POFF+VECSZ*373(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 92 KMUL
               lxv rA6, POFF+VECSZ*374(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 92 KMUL
               lxv rA7, POFF+VECSZ*375(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 94 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*281(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 94 KMUL
               lxv rB0, POFF+VECSZ*282(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 94 KMUL
               lxv rB1, POFF+VECSZ*283(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 94 KMUL
               lxv rB2, POFF+VECSZ*284(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 94 KMUL
               lxv rB3, POFF+VECSZ*285(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 94 KMUL
               lxv rB4, POFF+VECSZ*286(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 94 KMUL
               lxv rA0, POFF+VECSZ*376(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 94 KMUL
               lxv rA1, POFF+VECSZ*377(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 94 KMUL
               lxv rA2, POFF+VECSZ*378(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 94 KMUL
               lxv rA3, POFF+VECSZ*379(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 94 KMUL
               lxv rA4, POFF+VECSZ*380(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 94 KMUL
               lxv rA5, POFF+VECSZ*381(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 94 KMUL
               lxv rA6, POFF+VECSZ*382(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 94 KMUL
               lxv rA7, POFF+VECSZ*383(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 96 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*287(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 96 KMUL
               lxv rB0, POFF+VECSZ*288(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 96 KMUL
               lxv rB1, POFF+VECSZ*289(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 96 KMUL
               lxv rB2, POFF+VECSZ*290(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 96 KMUL
               lxv rB3, POFF+VECSZ*291(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 96 KMUL
               lxv rB4, POFF+VECSZ*292(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 96 KMUL
               lxv rA0, POFF+VECSZ*384(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 96 KMUL
               lxv rA1, POFF+VECSZ*385(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 96 KMUL
               lxv rA2, POFF+VECSZ*386(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 96 KMUL
               lxv rA3, POFF+VECSZ*387(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 96 KMUL
               lxv rA4, POFF+VECSZ*388(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 96 KMUL
               lxv rA5, POFF+VECSZ*389(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 96 KMUL
               lxv rA6, POFF+VECSZ*390(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 96 KMUL
               lxv rA7, POFF+VECSZ*391(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 98 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*293(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 98 KMUL
               lxv rB0, POFF+VECSZ*294(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 98 KMUL
               lxv rB1, POFF+VECSZ*295(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 98 KMUL
               lxv rB2, POFF+VECSZ*296(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 98 KMUL
               lxv rB3, POFF+VECSZ*297(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 98 KMUL
               lxv rB4, POFF+VECSZ*298(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 98 KMUL
               lxv rA0, POFF+VECSZ*392(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 98 KMUL
               lxv rA1, POFF+VECSZ*393(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 98 KMUL
               lxv rA2, POFF+VECSZ*394(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 98 KMUL
               lxv rA3, POFF+VECSZ*395(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 98 KMUL
               lxv rA4, POFF+VECSZ*396(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 98 KMUL
               lxv rA5, POFF+VECSZ*397(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 98 KMUL
               lxv rA6, POFF+VECSZ*398(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 98 KMUL
               lxv rA7, POFF+VECSZ*399(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 100 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*299(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 100 KMUL
               lxv rB0, POFF+VECSZ*300(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 100 KMUL
               lxv rB1, POFF+VECSZ*301(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 100 KMUL
               lxv rB2, POFF+VECSZ*302(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 100 KMUL
               lxv rB3, POFF+VECSZ*303(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 100 KMUL
               lxv rB4, POFF+VECSZ*304(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 100 KMUL
               lxv rA0, POFF+VECSZ*400(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 100 KMUL
               lxv rA1, POFF+VECSZ*401(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 100 KMUL
               lxv rA2, POFF+VECSZ*402(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 100 KMUL
               lxv rA3, POFF+VECSZ*403(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 100 KMUL
               lxv rA4, POFF+VECSZ*404(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 100 KMUL
               lxv rA5, POFF+VECSZ*405(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 100 KMUL
               lxv rA6, POFF+VECSZ*406(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 100 KMUL
               lxv rA7, POFF+VECSZ*407(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 102 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*305(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 102 KMUL
               lxv rB0, POFF+VECSZ*306(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 102 KMUL
               lxv rB1, POFF+VECSZ*307(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 102 KMUL
               lxv rB2, POFF+VECSZ*308(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 102 KMUL
               lxv rB3, POFF+VECSZ*309(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 102 KMUL
               lxv rB4, POFF+VECSZ*310(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 102 KMUL
               lxv rA0, POFF+VECSZ*408(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 102 KMUL
               lxv rA1, POFF+VECSZ*409(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 102 KMUL
               lxv rA2, POFF+VECSZ*410(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 102 KMUL
               lxv rA3, POFF+VECSZ*411(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 102 KMUL
               lxv rA4, POFF+VECSZ*412(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 102 KMUL
               lxv rA5, POFF+VECSZ*413(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 102 KMUL
               lxv rA6, POFF+VECSZ*414(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 102 KMUL
               lxv rA7, POFF+VECSZ*415(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 104 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*311(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 104 KMUL
               lxv rB0, POFF+VECSZ*312(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 104 KMUL
               lxv rB1, POFF+VECSZ*313(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 104 KMUL
               lxv rB2, POFF+VECSZ*314(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 104 KMUL
               lxv rB3, POFF+VECSZ*315(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 104 KMUL
               lxv rB4, POFF+VECSZ*316(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 104 KMUL
               lxv rA0, POFF+VECSZ*416(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 104 KMUL
               lxv rA1, POFF+VECSZ*417(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 104 KMUL
               lxv rA2, POFF+VECSZ*418(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 104 KMUL
               lxv rA3, POFF+VECSZ*419(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 104 KMUL
               lxv rA4, POFF+VECSZ*420(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 104 KMUL
               lxv rA5, POFF+VECSZ*421(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 104 KMUL
               lxv rA6, POFF+VECSZ*422(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 104 KMUL
               lxv rA7, POFF+VECSZ*423(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 106 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*317(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 106 KMUL
               lxv rB0, POFF+VECSZ*318(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 106 KMUL
               lxv rB1, POFF+VECSZ*319(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 106 KMUL
               lxv rB2, POFF+VECSZ*320(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 106 KMUL
               lxv rB3, POFF+VECSZ*321(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 106 KMUL
               lxv rB4, POFF+VECSZ*322(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 106 KMUL
               lxv rA0, POFF+VECSZ*424(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 106 KMUL
               lxv rA1, POFF+VECSZ*425(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 106 KMUL
               lxv rA2, POFF+VECSZ*426(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 106 KMUL
               lxv rA3, POFF+VECSZ*427(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 106 KMUL
               lxv rA4, POFF+VECSZ*428(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 106 KMUL
               lxv rA5, POFF+VECSZ*429(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 106 KMUL
               lxv rA6, POFF+VECSZ*430(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 106 KMUL
               lxv rA7, POFF+VECSZ*431(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 108 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*323(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 108 KMUL
               lxv rB0, POFF+VECSZ*324(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 108 KMUL
               lxv rB1, POFF+VECSZ*325(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 108 KMUL
               lxv rB2, POFF+VECSZ*326(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 108 KMUL
               lxv rB3, POFF+VECSZ*327(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 108 KMUL
               lxv rB4, POFF+VECSZ*328(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 108 KMUL
               lxv rA0, POFF+VECSZ*432(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 108 KMUL
               lxv rA1, POFF+VECSZ*433(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 108 KMUL
               lxv rA2, POFF+VECSZ*434(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 108 KMUL
               lxv rA3, POFF+VECSZ*435(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 108 KMUL
               lxv rA4, POFF+VECSZ*436(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 108 KMUL
               lxv rA5, POFF+VECSZ*437(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 108 KMUL
               lxv rA6, POFF+VECSZ*438(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 108 KMUL
               lxv rA7, POFF+VECSZ*439(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 110 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*329(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 110 KMUL
               lxv rB0, POFF+VECSZ*330(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 110 KMUL
               lxv rB1, POFF+VECSZ*331(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 110 KMUL
               lxv rB2, POFF+VECSZ*332(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 110 KMUL
               lxv rB3, POFF+VECSZ*333(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 110 KMUL
               lxv rB4, POFF+VECSZ*334(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 110 KMUL
               lxv rA0, POFF+VECSZ*440(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 110 KMUL
               lxv rA1, POFF+VECSZ*441(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 110 KMUL
               lxv rA2, POFF+VECSZ*442(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 110 KMUL
               lxv rA3, POFF+VECSZ*443(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 110 KMUL
               lxv rA4, POFF+VECSZ*444(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 110 KMUL
               lxv rA5, POFF+VECSZ*445(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 110 KMUL
               lxv rA6, POFF+VECSZ*446(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 110 KMUL
               lxv rA7, POFF+VECSZ*447(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 112 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*335(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 112 KMUL
               lxv rB0, POFF+VECSZ*336(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 112 KMUL
               lxv rB1, POFF+VECSZ*337(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 112 KMUL
               lxv rB2, POFF+VECSZ*338(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 112 KMUL
               lxv rB3, POFF+VECSZ*339(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 112 KMUL
               lxv rB4, POFF+VECSZ*340(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 112 KMUL
               lxv rA0, POFF+VECSZ*448(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 112 KMUL
               lxv rA1, POFF+VECSZ*449(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 112 KMUL
               lxv rA2, POFF+VECSZ*450(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 112 KMUL
               lxv rA3, POFF+VECSZ*451(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 112 KMUL
               lxv rA4, POFF+VECSZ*452(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 112 KMUL
               lxv rA5, POFF+VECSZ*453(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 112 KMUL
               lxv rA6, POFF+VECSZ*454(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 112 KMUL
               lxv rA7, POFF+VECSZ*455(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 114 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*341(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 114 KMUL
               lxv rB0, POFF+VECSZ*342(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 114 KMUL
               lxv rB1, POFF+VECSZ*343(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 114 KMUL
               lxv rB2, POFF+VECSZ*344(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 114 KMUL
               lxv rB3, POFF+VECSZ*345(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 114 KMUL
               lxv rB4, POFF+VECSZ*346(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 114 KMUL
               lxv rA0, POFF+VECSZ*456(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 114 KMUL
               lxv rA1, POFF+VECSZ*457(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 114 KMUL
               lxv rA2, POFF+VECSZ*458(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 114 KMUL
               lxv rA3, POFF+VECSZ*459(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 114 KMUL
               lxv rA4, POFF+VECSZ*460(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 114 KMUL
               lxv rA5, POFF+VECSZ*461(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 114 KMUL
               lxv rA6, POFF+VECSZ*462(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 114 KMUL
               lxv rA7, POFF+VECSZ*463(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 116 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*347(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 116 KMUL
               lxv rB0, POFF+VECSZ*348(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 116 KMUL
               lxv rB1, POFF+VECSZ*349(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 116 KMUL
               lxv rB2, POFF+VECSZ*350(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 116 KMUL
               lxv rB3, POFF+VECSZ*351(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 116 KMUL
               lxv rB4, POFF+VECSZ*352(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 116 KMUL
               lxv rA0, POFF+VECSZ*464(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 116 KMUL
               lxv rA1, POFF+VECSZ*465(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 116 KMUL
               lxv rA2, POFF+VECSZ*466(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 116 KMUL
               lxv rA3, POFF+VECSZ*467(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 116 KMUL
               lxv rA4, POFF+VECSZ*468(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 116 KMUL
               lxv rA5, POFF+VECSZ*469(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 116 KMUL
               lxv rA6, POFF+VECSZ*470(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 116 KMUL
               lxv rA7, POFF+VECSZ*471(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 118 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*353(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 118 KMUL
               lxv rB0, POFF+VECSZ*354(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 118 KMUL
               lxv rB1, POFF+VECSZ*355(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 118 KMUL
               lxv rB2, POFF+VECSZ*356(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 118 KMUL
               lxv rB3, POFF+VECSZ*357(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 118 KMUL
               lxv rB4, POFF+VECSZ*358(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 118 KMUL
               lxv rA0, POFF+VECSZ*472(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 118 KMUL
               lxv rA1, POFF+VECSZ*473(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 118 KMUL
               lxv rA2, POFF+VECSZ*474(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 118 KMUL
               lxv rA3, POFF+VECSZ*475(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 118 KMUL
               lxv rA4, POFF+VECSZ*476(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 118 KMUL
               lxv rA5, POFF+VECSZ*477(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 118 KMUL
               lxv rA6, POFF+VECSZ*478(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 118 KMUL
               lxv rA7, POFF+VECSZ*479(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 120 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*359(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 120 KMUL
               lxv rB0, POFF+VECSZ*360(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 120 KMUL
               lxv rB1, POFF+VECSZ*361(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 120 KMUL
               lxv rB2, POFF+VECSZ*362(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 120 KMUL
               lxv rB3, POFF+VECSZ*363(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 120 KMUL
               lxv rB4, POFF+VECSZ*364(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 120 KMUL
               lxv rA0, POFF+VECSZ*480(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 120 KMUL
               lxv rA1, POFF+VECSZ*481(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 120 KMUL
               lxv rA2, POFF+VECSZ*482(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 120 KMUL
               lxv rA3, POFF+VECSZ*483(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 120 KMUL
               lxv rA4, POFF+VECSZ*484(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 120 KMUL
               lxv rA5, POFF+VECSZ*485(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 120 KMUL
               lxv rA6, POFF+VECSZ*486(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 120 KMUL
               lxv rA7, POFF+VECSZ*487(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 122 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*365(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 122 KMUL
               lxv rB0, POFF+VECSZ*366(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 122 KMUL
               lxv rB1, POFF+VECSZ*367(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 122 KMUL
               lxv rB2, POFF+VECSZ*368(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 122 KMUL
               lxv rB3, POFF+VECSZ*369(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 122 KMUL
               lxv rB4, POFF+VECSZ*370(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 122 KMUL
               lxv rA0, POFF+VECSZ*488(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 122 KMUL
               lxv rA1, POFF+VECSZ*489(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 122 KMUL
               lxv rA2, POFF+VECSZ*490(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 122 KMUL
               lxv rA3, POFF+VECSZ*491(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 122 KMUL
               lxv rA4, POFF+VECSZ*492(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 122 KMUL
               lxv rA5, POFF+VECSZ*493(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 122 KMUL
               lxv rA6, POFF+VECSZ*494(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 122 KMUL
               lxv rA7, POFF+VECSZ*495(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 124 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*371(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 124 KMUL
               lxv rB0, POFF+VECSZ*372(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 124 KMUL
               lxv rB1, POFF+VECSZ*373(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 124 KMUL
               lxv rB2, POFF+VECSZ*374(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 124 KMUL
               lxv rB3, POFF+VECSZ*375(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 124 KMUL
               lxv rB4, POFF+VECSZ*376(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 124 KMUL
               lxv rA0, POFF+VECSZ*496(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 124 KMUL
               lxv rA1, POFF+VECSZ*497(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 124 KMUL
               lxv rA2, POFF+VECSZ*498(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 124 KMUL
               lxv rA3, POFF+VECSZ*499(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 124 KMUL
               lxv rA4, POFF+VECSZ*500(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 124 KMUL
               lxv rA5, POFF+VECSZ*501(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 124 KMUL
               lxv rA6, POFF+VECSZ*502(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 124 KMUL
               lxv rA7, POFF+VECSZ*503(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 126 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*377(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 126 KMUL
               lxv rB0, POFF+VECSZ*378(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 126 KMUL
               lxv rB1, POFF+VECSZ*379(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 126 KMUL
               lxv rB2, POFF+VECSZ*380(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 126 KMUL
               lxv rB3, POFF+VECSZ*381(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 126 KMUL
               lxv rB4, POFF+VECSZ*382(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 126 KMUL
               lxv rA0, POFF+VECSZ*504(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 126 KMUL
               lxv rA1, POFF+VECSZ*505(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 126 KMUL
               lxv rA2, POFF+VECSZ*506(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 126 KMUL
               lxv rA3, POFF+VECSZ*507(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 126 KMUL
               lxv rA4, POFF+VECSZ*508(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 126 KMUL
               lxv rA5, POFF+VECSZ*509(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 126 KMUL
               lxv rA6, POFF+VECSZ*510(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 126 KMUL
               lxv rA7, POFF+VECSZ*511(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 128 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*383(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 128 KMUL
               lxv rB0, POFF+VECSZ*384(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 128 KMUL
               lxv rB1, POFF+VECSZ*385(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 128 KMUL
               lxv rB2, POFF+VECSZ*386(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 128 KMUL
               lxv rB3, POFF+VECSZ*387(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 128 KMUL
               lxv rB4, POFF+VECSZ*388(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 128 KMUL
               lxv rA0, POFF+VECSZ*512(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 128 KMUL
               lxv rA1, POFF+VECSZ*513(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 128 KMUL
               lxv rA2, POFF+VECSZ*514(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 128 KMUL
               lxv rA3, POFF+VECSZ*515(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 128 KMUL
               lxv rA4, POFF+VECSZ*516(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 128 KMUL
               lxv rA5, POFF+VECSZ*517(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 128 KMUL
               lxv rA6, POFF+VECSZ*518(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 128 KMUL
               lxv rA7, POFF+VECSZ*519(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 130 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*389(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 130 KMUL
               lxv rB0, POFF+VECSZ*390(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 130 KMUL
               lxv rB1, POFF+VECSZ*391(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 130 KMUL
               lxv rB2, POFF+VECSZ*392(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 130 KMUL
               lxv rB3, POFF+VECSZ*393(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 130 KMUL
               lxv rB4, POFF+VECSZ*394(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 130 KMUL
               lxv rA0, POFF+VECSZ*520(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 130 KMUL
               lxv rA1, POFF+VECSZ*521(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 130 KMUL
               lxv rA2, POFF+VECSZ*522(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 130 KMUL
               lxv rA3, POFF+VECSZ*523(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 130 KMUL
               lxv rA4, POFF+VECSZ*524(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 130 KMUL
               lxv rA5, POFF+VECSZ*525(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 130 KMUL
               lxv rA6, POFF+VECSZ*526(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 130 KMUL
               lxv rA7, POFF+VECSZ*527(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 132 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*395(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 132 KMUL
               lxv rB0, POFF+VECSZ*396(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 132 KMUL
               lxv rB1, POFF+VECSZ*397(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 132 KMUL
               lxv rB2, POFF+VECSZ*398(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 132 KMUL
               lxv rB3, POFF+VECSZ*399(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 132 KMUL
               lxv rB4, POFF+VECSZ*400(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 132 KMUL
               lxv rA0, POFF+VECSZ*528(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 132 KMUL
               lxv rA1, POFF+VECSZ*529(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 132 KMUL
               lxv rA2, POFF+VECSZ*530(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 132 KMUL
               lxv rA3, POFF+VECSZ*531(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 132 KMUL
               lxv rA4, POFF+VECSZ*532(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 132 KMUL
               lxv rA5, POFF+VECSZ*533(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 132 KMUL
               lxv rA6, POFF+VECSZ*534(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 132 KMUL
               lxv rA7, POFF+VECSZ*535(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 134 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*401(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 134 KMUL
               lxv rB0, POFF+VECSZ*402(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 134 KMUL
               lxv rB1, POFF+VECSZ*403(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 134 KMUL
               lxv rB2, POFF+VECSZ*404(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 134 KMUL
               lxv rB3, POFF+VECSZ*405(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 134 KMUL
               lxv rB4, POFF+VECSZ*406(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 134 KMUL
               lxv rA0, POFF+VECSZ*536(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 134 KMUL
               lxv rA1, POFF+VECSZ*537(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 134 KMUL
               lxv rA2, POFF+VECSZ*538(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 134 KMUL
               lxv rA3, POFF+VECSZ*539(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 134 KMUL
               lxv rA4, POFF+VECSZ*540(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 134 KMUL
               lxv rA5, POFF+VECSZ*541(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 134 KMUL
               lxv rA6, POFF+VECSZ*542(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 134 KMUL
               lxv rA7, POFF+VECSZ*543(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 136 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*407(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 136 KMUL
               lxv rB0, POFF+VECSZ*408(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 136 KMUL
               lxv rB1, POFF+VECSZ*409(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 136 KMUL
               lxv rB2, POFF+VECSZ*410(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 136 KMUL
               lxv rB3, POFF+VECSZ*411(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 136 KMUL
               lxv rB4, POFF+VECSZ*412(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 136 KMUL
               lxv rA0, POFF+VECSZ*544(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 136 KMUL
               lxv rA1, POFF+VECSZ*545(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 136 KMUL
               lxv rA2, POFF+VECSZ*546(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 136 KMUL
               lxv rA3, POFF+VECSZ*547(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 136 KMUL
               lxv rA4, POFF+VECSZ*548(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 136 KMUL
               lxv rA5, POFF+VECSZ*549(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 136 KMUL
               lxv rA6, POFF+VECSZ*550(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 136 KMUL
               lxv rA7, POFF+VECSZ*551(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 138 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*413(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 138 KMUL
               lxv rB0, POFF+VECSZ*414(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 138 KMUL
               lxv rB1, POFF+VECSZ*415(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 138 KMUL
               lxv rB2, POFF+VECSZ*416(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 138 KMUL
               lxv rB3, POFF+VECSZ*417(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 138 KMUL
               lxv rB4, POFF+VECSZ*418(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 138 KMUL
               lxv rA0, POFF+VECSZ*552(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 138 KMUL
               lxv rA1, POFF+VECSZ*553(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 138 KMUL
               lxv rA2, POFF+VECSZ*554(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 138 KMUL
               lxv rA3, POFF+VECSZ*555(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 138 KMUL
               lxv rA4, POFF+VECSZ*556(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 138 KMUL
               lxv rA5, POFF+VECSZ*557(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 138 KMUL
               lxv rA6, POFF+VECSZ*558(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 138 KMUL
               lxv rA7, POFF+VECSZ*559(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 140 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*419(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 140 KMUL
               lxv rB0, POFF+VECSZ*420(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 140 KMUL
               lxv rB1, POFF+VECSZ*421(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 140 KMUL
               lxv rB2, POFF+VECSZ*422(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 140 KMUL
               lxv rB3, POFF+VECSZ*423(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 140 KMUL
               lxv rB4, POFF+VECSZ*424(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 140 KMUL
               lxv rA0, POFF+VECSZ*560(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 140 KMUL
               lxv rA1, POFF+VECSZ*561(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 140 KMUL
               lxv rA2, POFF+VECSZ*562(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 140 KMUL
               lxv rA3, POFF+VECSZ*563(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 140 KMUL
               lxv rA4, POFF+VECSZ*564(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 140 KMUL
               lxv rA5, POFF+VECSZ*565(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 140 KMUL
               lxv rA6, POFF+VECSZ*566(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 140 KMUL
               lxv rA7, POFF+VECSZ*567(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 142 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*425(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 142 KMUL
               lxv rB0, POFF+VECSZ*426(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 142 KMUL
               lxv rB1, POFF+VECSZ*427(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 142 KMUL
               lxv rB2, POFF+VECSZ*428(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 142 KMUL
               lxv rB3, POFF+VECSZ*429(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 142 KMUL
               lxv rB4, POFF+VECSZ*430(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 142 KMUL
               lxv rA0, POFF+VECSZ*568(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 142 KMUL
               lxv rA1, POFF+VECSZ*569(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 142 KMUL
               lxv rA2, POFF+VECSZ*570(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 142 KMUL
               lxv rA3, POFF+VECSZ*571(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 142 KMUL
               lxv rA4, POFF+VECSZ*572(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 142 KMUL
               lxv rA5, POFF+VECSZ*573(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 142 KMUL
               lxv rA6, POFF+VECSZ*574(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 142 KMUL
               lxv rA7, POFF+VECSZ*575(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 144 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*431(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 144 KMUL
               lxv rB0, POFF+VECSZ*432(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 144 KMUL
               lxv rB1, POFF+VECSZ*433(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 144 KMUL
               lxv rB2, POFF+VECSZ*434(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 144 KMUL
               lxv rB3, POFF+VECSZ*435(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 144 KMUL
               lxv rB4, POFF+VECSZ*436(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 144 KMUL
               lxv rA0, POFF+VECSZ*576(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 144 KMUL
               lxv rA1, POFF+VECSZ*577(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 144 KMUL
               lxv rA2, POFF+VECSZ*578(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 144 KMUL
               lxv rA3, POFF+VECSZ*579(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 144 KMUL
               lxv rA4, POFF+VECSZ*580(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 144 KMUL
               lxv rA5, POFF+VECSZ*581(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 144 KMUL
               lxv rA6, POFF+VECSZ*582(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 144 KMUL
               lxv rA7, POFF+VECSZ*583(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 146 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*437(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 146 KMUL
               lxv rB0, POFF+VECSZ*438(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 146 KMUL
               lxv rB1, POFF+VECSZ*439(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 146 KMUL
               lxv rB2, POFF+VECSZ*440(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 146 KMUL
               lxv rB3, POFF+VECSZ*441(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 146 KMUL
               lxv rB4, POFF+VECSZ*442(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 146 KMUL
               lxv rA0, POFF+VECSZ*584(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 146 KMUL
               lxv rA1, POFF+VECSZ*585(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 146 KMUL
               lxv rA2, POFF+VECSZ*586(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 146 KMUL
               lxv rA3, POFF+VECSZ*587(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 146 KMUL
               lxv rA4, POFF+VECSZ*588(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 146 KMUL
               lxv rA5, POFF+VECSZ*589(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 146 KMUL
               lxv rA6, POFF+VECSZ*590(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 146 KMUL
               lxv rA7, POFF+VECSZ*591(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 148 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*443(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 148 KMUL
               lxv rB0, POFF+VECSZ*444(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 148 KMUL
               lxv rB1, POFF+VECSZ*445(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 148 KMUL
               lxv rB2, POFF+VECSZ*446(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 148 KMUL
               lxv rB3, POFF+VECSZ*447(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 148 KMUL
               lxv rB4, POFF+VECSZ*448(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 148 KMUL
               lxv rA0, POFF+VECSZ*592(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 148 KMUL
               lxv rA1, POFF+VECSZ*593(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 148 KMUL
               lxv rA2, POFF+VECSZ*594(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 148 KMUL
               lxv rA3, POFF+VECSZ*595(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 148 KMUL
               lxv rA4, POFF+VECSZ*596(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 148 KMUL
               lxv rA5, POFF+VECSZ*597(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 148 KMUL
               lxv rA6, POFF+VECSZ*598(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 148 KMUL
               lxv rA7, POFF+VECSZ*599(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 150 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*449(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 150 KMUL
               lxv rB0, POFF+VECSZ*450(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 150 KMUL
               lxv rB1, POFF+VECSZ*451(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 150 KMUL
               lxv rB2, POFF+VECSZ*452(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 150 KMUL
               lxv rB3, POFF+VECSZ*453(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 150 KMUL
               lxv rB4, POFF+VECSZ*454(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 150 KMUL
               lxv rA0, POFF+VECSZ*600(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 150 KMUL
               lxv rA1, POFF+VECSZ*601(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 150 KMUL
               lxv rA2, POFF+VECSZ*602(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 150 KMUL
               lxv rA3, POFF+VECSZ*603(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 150 KMUL
               lxv rA4, POFF+VECSZ*604(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 150 KMUL
               lxv rA5, POFF+VECSZ*605(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 150 KMUL
               lxv rA6, POFF+VECSZ*606(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 150 KMUL
               lxv rA7, POFF+VECSZ*607(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 152 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*455(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 152 KMUL
               lxv rB0, POFF+VECSZ*456(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 152 KMUL
               lxv rB1, POFF+VECSZ*457(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 152 KMUL
               lxv rB2, POFF+VECSZ*458(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 152 KMUL
               lxv rB3, POFF+VECSZ*459(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 152 KMUL
               lxv rB4, POFF+VECSZ*460(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 152 KMUL
               lxv rA0, POFF+VECSZ*608(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 152 KMUL
               lxv rA1, POFF+VECSZ*609(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 152 KMUL
               lxv rA2, POFF+VECSZ*610(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 152 KMUL
               lxv rA3, POFF+VECSZ*611(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 152 KMUL
               lxv rA4, POFF+VECSZ*612(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 152 KMUL
               lxv rA5, POFF+VECSZ*613(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 152 KMUL
               lxv rA6, POFF+VECSZ*614(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 152 KMUL
               lxv rA7, POFF+VECSZ*615(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 154 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*461(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 154 KMUL
               lxv rB0, POFF+VECSZ*462(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 154 KMUL
               lxv rB1, POFF+VECSZ*463(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 154 KMUL
               lxv rB2, POFF+VECSZ*464(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 154 KMUL
               lxv rB3, POFF+VECSZ*465(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 154 KMUL
               lxv rB4, POFF+VECSZ*466(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 154 KMUL
               lxv rA0, POFF+VECSZ*616(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 154 KMUL
               lxv rA1, POFF+VECSZ*617(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 154 KMUL
               lxv rA2, POFF+VECSZ*618(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 154 KMUL
               lxv rA3, POFF+VECSZ*619(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 154 KMUL
               lxv rA4, POFF+VECSZ*620(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 154 KMUL
               lxv rA5, POFF+VECSZ*621(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 154 KMUL
               lxv rA6, POFF+VECSZ*622(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 154 KMUL
               lxv rA7, POFF+VECSZ*623(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 156 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*467(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 156 KMUL
               lxv rB0, POFF+VECSZ*468(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 156 KMUL
               lxv rB1, POFF+VECSZ*469(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 156 KMUL
               lxv rB2, POFF+VECSZ*470(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 156 KMUL
               lxv rB3, POFF+VECSZ*471(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 156 KMUL
               lxv rB4, POFF+VECSZ*472(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 156 KMUL
               lxv rA0, POFF+VECSZ*624(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 156 KMUL
               lxv rA1, POFF+VECSZ*625(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 156 KMUL
               lxv rA2, POFF+VECSZ*626(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 156 KMUL
               lxv rA3, POFF+VECSZ*627(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 156 KMUL
               lxv rA4, POFF+VECSZ*628(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 156 KMUL
               lxv rA5, POFF+VECSZ*629(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 156 KMUL
               lxv rA6, POFF+VECSZ*630(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 156 KMUL
               lxv rA7, POFF+VECSZ*631(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 158 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*473(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 158 KMUL
               lxv rB0, POFF+VECSZ*474(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 158 KMUL
               lxv rB1, POFF+VECSZ*475(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 158 KMUL
               lxv rB2, POFF+VECSZ*476(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 158 KMUL
               lxv rB3, POFF+VECSZ*477(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 158 KMUL
               lxv rB4, POFF+VECSZ*478(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 158 KMUL
               lxv rA0, POFF+VECSZ*632(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 158 KMUL
               lxv rA1, POFF+VECSZ*633(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 158 KMUL
               lxv rA2, POFF+VECSZ*634(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 158 KMUL
               lxv rA3, POFF+VECSZ*635(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 158 KMUL
               lxv rA4, POFF+VECSZ*636(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 158 KMUL
               lxv rA5, POFF+VECSZ*637(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 158 KMUL
               lxv rA6, POFF+VECSZ*638(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 158 KMUL
               lxv rA7, POFF+VECSZ*639(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 160 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*479(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 160 KMUL
               lxv rB0, POFF+VECSZ*480(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 160 KMUL
               lxv rB1, POFF+VECSZ*481(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 160 KMUL
               lxv rB2, POFF+VECSZ*482(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 160 KMUL
               lxv rB3, POFF+VECSZ*483(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 160 KMUL
               lxv rB4, POFF+VECSZ*484(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 160 KMUL
               lxv rA0, POFF+VECSZ*640(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 160 KMUL
               lxv rA1, POFF+VECSZ*641(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 160 KMUL
               lxv rA2, POFF+VECSZ*642(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 160 KMUL
               lxv rA3, POFF+VECSZ*643(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 160 KMUL
               lxv rA4, POFF+VECSZ*644(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 160 KMUL
               lxv rA5, POFF+VECSZ*645(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 160 KMUL
               lxv rA6, POFF+VECSZ*646(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 160 KMUL
               lxv rA7, POFF+VECSZ*647(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 162 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*485(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 162 KMUL
               lxv rB0, POFF+VECSZ*486(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 162 KMUL
               lxv rB1, POFF+VECSZ*487(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 162 KMUL
               lxv rB2, POFF+VECSZ*488(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 162 KMUL
               lxv rB3, POFF+VECSZ*489(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 162 KMUL
               lxv rB4, POFF+VECSZ*490(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 162 KMUL
               lxv rA0, POFF+VECSZ*648(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 162 KMUL
               lxv rA1, POFF+VECSZ*649(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 162 KMUL
               lxv rA2, POFF+VECSZ*650(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 162 KMUL
               lxv rA3, POFF+VECSZ*651(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 162 KMUL
               lxv rA4, POFF+VECSZ*652(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 162 KMUL
               lxv rA5, POFF+VECSZ*653(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 162 KMUL
               lxv rA6, POFF+VECSZ*654(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 162 KMUL
               lxv rA7, POFF+VECSZ*655(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 164 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*491(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 164 KMUL
               lxv rB0, POFF+VECSZ*492(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 164 KMUL
               lxv rB1, POFF+VECSZ*493(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 164 KMUL
               lxv rB2, POFF+VECSZ*494(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 164 KMUL
               lxv rB3, POFF+VECSZ*495(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 164 KMUL
               lxv rB4, POFF+VECSZ*496(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 164 KMUL
               lxv rA0, POFF+VECSZ*656(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 164 KMUL
               lxv rA1, POFF+VECSZ*657(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 164 KMUL
               lxv rA2, POFF+VECSZ*658(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 164 KMUL
               lxv rA3, POFF+VECSZ*659(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 164 KMUL
               lxv rA4, POFF+VECSZ*660(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 164 KMUL
               lxv rA5, POFF+VECSZ*661(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 164 KMUL
               lxv rA6, POFF+VECSZ*662(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 164 KMUL
               lxv rA7, POFF+VECSZ*663(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 166 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*497(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 166 KMUL
               lxv rB0, POFF+VECSZ*498(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 166 KMUL
               lxv rB1, POFF+VECSZ*499(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 166 KMUL
               lxv rB2, POFF+VECSZ*500(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 166 KMUL
               lxv rB3, POFF+VECSZ*501(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 166 KMUL
               lxv rB4, POFF+VECSZ*502(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 166 KMUL
               lxv rA0, POFF+VECSZ*664(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 166 KMUL
               lxv rA1, POFF+VECSZ*665(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 166 KMUL
               lxv rA2, POFF+VECSZ*666(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 166 KMUL
               lxv rA3, POFF+VECSZ*667(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 166 KMUL
               lxv rA4, POFF+VECSZ*668(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 166 KMUL
               lxv rA5, POFF+VECSZ*669(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 166 KMUL
               lxv rA6, POFF+VECSZ*670(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 166 KMUL
               lxv rA7, POFF+VECSZ*671(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 168 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*503(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 168 KMUL
               lxv rB0, POFF+VECSZ*504(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 168 KMUL
               lxv rB1, POFF+VECSZ*505(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 168 KMUL
               lxv rB2, POFF+VECSZ*506(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 168 KMUL
               lxv rB3, POFF+VECSZ*507(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 168 KMUL
               lxv rB4, POFF+VECSZ*508(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 168 KMUL
               lxv rA0, POFF+VECSZ*672(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 168 KMUL
               lxv rA1, POFF+VECSZ*673(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 168 KMUL
               lxv rA2, POFF+VECSZ*674(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 168 KMUL
               lxv rA3, POFF+VECSZ*675(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 168 KMUL
               lxv rA4, POFF+VECSZ*676(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 168 KMUL
               lxv rA5, POFF+VECSZ*677(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 168 KMUL
               lxv rA6, POFF+VECSZ*678(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 168 KMUL
               lxv rA7, POFF+VECSZ*679(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 170 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*509(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 170 KMUL
               lxv rB0, POFF+VECSZ*510(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 170 KMUL
               lxv rB1, POFF+VECSZ*511(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 170 KMUL
               lxv rB2, POFF+VECSZ*512(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 170 KMUL
               lxv rB3, POFF+VECSZ*513(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 170 KMUL
               lxv rB4, POFF+VECSZ*514(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 170 KMUL
               lxv rA0, POFF+VECSZ*680(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 170 KMUL
               lxv rA1, POFF+VECSZ*681(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 170 KMUL
               lxv rA2, POFF+VECSZ*682(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 170 KMUL
               lxv rA3, POFF+VECSZ*683(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 170 KMUL
               lxv rA4, POFF+VECSZ*684(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 170 KMUL
               lxv rA5, POFF+VECSZ*685(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 170 KMUL
               lxv rA6, POFF+VECSZ*686(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 170 KMUL
               lxv rA7, POFF+VECSZ*687(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 172 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*515(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 172 KMUL
               lxv rB0, POFF+VECSZ*516(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 172 KMUL
               lxv rB1, POFF+VECSZ*517(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 172 KMUL
               lxv rB2, POFF+VECSZ*518(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 172 KMUL
               lxv rB3, POFF+VECSZ*519(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 172 KMUL
               lxv rB4, POFF+VECSZ*520(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 172 KMUL
               lxv rA0, POFF+VECSZ*688(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 172 KMUL
               lxv rA1, POFF+VECSZ*689(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 172 KMUL
               lxv rA2, POFF+VECSZ*690(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 172 KMUL
               lxv rA3, POFF+VECSZ*691(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 172 KMUL
               lxv rA4, POFF+VECSZ*692(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 172 KMUL
               lxv rA5, POFF+VECSZ*693(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 172 KMUL
               lxv rA6, POFF+VECSZ*694(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 172 KMUL
               lxv rA7, POFF+VECSZ*695(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 174 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*521(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 174 KMUL
               lxv rB0, POFF+VECSZ*522(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 174 KMUL
               lxv rB1, POFF+VECSZ*523(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 174 KMUL
               lxv rB2, POFF+VECSZ*524(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 174 KMUL
               lxv rB3, POFF+VECSZ*525(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 174 KMUL
               lxv rB4, POFF+VECSZ*526(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 174 KMUL
               lxv rA0, POFF+VECSZ*696(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 174 KMUL
               lxv rA1, POFF+VECSZ*697(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 174 KMUL
               lxv rA2, POFF+VECSZ*698(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 174 KMUL
               lxv rA3, POFF+VECSZ*699(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 174 KMUL
               lxv rA4, POFF+VECSZ*700(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 174 KMUL
               lxv rA5, POFF+VECSZ*701(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 174 KMUL
               lxv rA6, POFF+VECSZ*702(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 174 KMUL
               lxv rA7, POFF+VECSZ*703(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 176 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*527(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 176 KMUL
               lxv rB0, POFF+VECSZ*528(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 176 KMUL
               lxv rB1, POFF+VECSZ*529(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 176 KMUL
               lxv rB2, POFF+VECSZ*530(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 176 KMUL
               lxv rB3, POFF+VECSZ*531(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 176 KMUL
               lxv rB4, POFF+VECSZ*532(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 176 KMUL
               lxv rA0, POFF+VECSZ*704(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 176 KMUL
               lxv rA1, POFF+VECSZ*705(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 176 KMUL
               lxv rA2, POFF+VECSZ*706(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 176 KMUL
               lxv rA3, POFF+VECSZ*707(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 176 KMUL
               lxv rA4, POFF+VECSZ*708(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 176 KMUL
               lxv rA5, POFF+VECSZ*709(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 176 KMUL
               lxv rA6, POFF+VECSZ*710(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 176 KMUL
               lxv rA7, POFF+VECSZ*711(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 178 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*533(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 178 KMUL
               lxv rB0, POFF+VECSZ*534(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 178 KMUL
               lxv rB1, POFF+VECSZ*535(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 178 KMUL
               lxv rB2, POFF+VECSZ*536(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 178 KMUL
               lxv rB3, POFF+VECSZ*537(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 178 KMUL
               lxv rB4, POFF+VECSZ*538(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 178 KMUL
               lxv rA0, POFF+VECSZ*712(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 178 KMUL
               lxv rA1, POFF+VECSZ*713(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 178 KMUL
               lxv rA2, POFF+VECSZ*714(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 178 KMUL
               lxv rA3, POFF+VECSZ*715(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 178 KMUL
               lxv rA4, POFF+VECSZ*716(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 178 KMUL
               lxv rA5, POFF+VECSZ*717(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 178 KMUL
               lxv rA6, POFF+VECSZ*718(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 178 KMUL
               lxv rA7, POFF+VECSZ*719(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 180 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*539(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 180 KMUL
               lxv rB0, POFF+VECSZ*540(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 180 KMUL
               lxv rB1, POFF+VECSZ*541(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 180 KMUL
               lxv rB2, POFF+VECSZ*542(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 180 KMUL
               lxv rB3, POFF+VECSZ*543(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 180 KMUL
               lxv rB4, POFF+VECSZ*544(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 180 KMUL
               lxv rA0, POFF+VECSZ*720(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 180 KMUL
               lxv rA1, POFF+VECSZ*721(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 180 KMUL
               lxv rA2, POFF+VECSZ*722(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 180 KMUL
               lxv rA3, POFF+VECSZ*723(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 180 KMUL
               lxv rA4, POFF+VECSZ*724(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 180 KMUL
               lxv rA5, POFF+VECSZ*725(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 180 KMUL
               lxv rA6, POFF+VECSZ*726(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 180 KMUL
               lxv rA7, POFF+VECSZ*727(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 182 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*545(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 182 KMUL
               lxv rB0, POFF+VECSZ*546(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 182 KMUL
               lxv rB1, POFF+VECSZ*547(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 182 KMUL
               lxv rB2, POFF+VECSZ*548(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 182 KMUL
               lxv rB3, POFF+VECSZ*549(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 182 KMUL
               lxv rB4, POFF+VECSZ*550(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 182 KMUL
               lxv rA0, POFF+VECSZ*728(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 182 KMUL
               lxv rA1, POFF+VECSZ*729(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 182 KMUL
               lxv rA2, POFF+VECSZ*730(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 182 KMUL
               lxv rA3, POFF+VECSZ*731(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 182 KMUL
               lxv rA4, POFF+VECSZ*732(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 182 KMUL
               lxv rA5, POFF+VECSZ*733(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 182 KMUL
               lxv rA6, POFF+VECSZ*734(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 182 KMUL
               lxv rA7, POFF+VECSZ*735(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 184 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*551(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 184 KMUL
               lxv rB0, POFF+VECSZ*552(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 184 KMUL
               lxv rB1, POFF+VECSZ*553(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 184 KMUL
               lxv rB2, POFF+VECSZ*554(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 184 KMUL
               lxv rB3, POFF+VECSZ*555(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 184 KMUL
               lxv rB4, POFF+VECSZ*556(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 184 KMUL
               lxv rA0, POFF+VECSZ*736(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 184 KMUL
               lxv rA1, POFF+VECSZ*737(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 184 KMUL
               lxv rA2, POFF+VECSZ*738(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 184 KMUL
               lxv rA3, POFF+VECSZ*739(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 184 KMUL
               lxv rA4, POFF+VECSZ*740(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 184 KMUL
               lxv rA5, POFF+VECSZ*741(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 184 KMUL
               lxv rA6, POFF+VECSZ*742(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 184 KMUL
               lxv rA7, POFF+VECSZ*743(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 186 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*557(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 186 KMUL
               lxv rB0, POFF+VECSZ*558(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 186 KMUL
               lxv rB1, POFF+VECSZ*559(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 186 KMUL
               lxv rB2, POFF+VECSZ*560(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 186 KMUL
               lxv rB3, POFF+VECSZ*561(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 186 KMUL
               lxv rB4, POFF+VECSZ*562(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 186 KMUL
               lxv rA0, POFF+VECSZ*744(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 186 KMUL
               lxv rA1, POFF+VECSZ*745(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 186 KMUL
               lxv rA2, POFF+VECSZ*746(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 186 KMUL
               lxv rA3, POFF+VECSZ*747(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 186 KMUL
               lxv rA4, POFF+VECSZ*748(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 186 KMUL
               lxv rA5, POFF+VECSZ*749(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 186 KMUL
               lxv rA6, POFF+VECSZ*750(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 186 KMUL
               lxv rA7, POFF+VECSZ*751(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 188 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*563(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 188 KMUL
               lxv rB0, POFF+VECSZ*564(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 188 KMUL
               lxv rB1, POFF+VECSZ*565(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 188 KMUL
               lxv rB2, POFF+VECSZ*566(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 188 KMUL
               lxv rB3, POFF+VECSZ*567(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 188 KMUL
               lxv rB4, POFF+VECSZ*568(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 188 KMUL
               lxv rA0, POFF+VECSZ*752(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 188 KMUL
               lxv rA1, POFF+VECSZ*753(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 188 KMUL
               lxv rA2, POFF+VECSZ*754(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 188 KMUL
               lxv rA3, POFF+VECSZ*755(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 188 KMUL
               lxv rA4, POFF+VECSZ*756(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 188 KMUL
               lxv rA5, POFF+VECSZ*757(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 188 KMUL
               lxv rA6, POFF+VECSZ*758(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 188 KMUL
               lxv rA7, POFF+VECSZ*759(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 190 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*569(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 190 KMUL
               lxv rB0, POFF+VECSZ*570(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 190 KMUL
               lxv rB1, POFF+VECSZ*571(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 190 KMUL
               lxv rB2, POFF+VECSZ*572(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 190 KMUL
               lxv rB3, POFF+VECSZ*573(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 190 KMUL
               lxv rB4, POFF+VECSZ*574(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 190 KMUL
               lxv rA0, POFF+VECSZ*760(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 190 KMUL
               lxv rA1, POFF+VECSZ*761(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 190 KMUL
               lxv rA2, POFF+VECSZ*762(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 190 KMUL
               lxv rA3, POFF+VECSZ*763(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 190 KMUL
               lxv rA4, POFF+VECSZ*764(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 190 KMUL
               lxv rA5, POFF+VECSZ*765(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 190 KMUL
               lxv rA6, POFF+VECSZ*766(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 190 KMUL
               lxv rA7, POFF+VECSZ*767(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 192 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*575(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 192 KMUL
               lxv rB0, POFF+VECSZ*576(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 192 KMUL
               lxv rB1, POFF+VECSZ*577(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 192 KMUL
               lxv rB2, POFF+VECSZ*578(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 192 KMUL
               lxv rB3, POFF+VECSZ*579(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 192 KMUL
               lxv rB4, POFF+VECSZ*580(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 192 KMUL
               lxv rA0, POFF+VECSZ*768(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 192 KMUL
               lxv rA1, POFF+VECSZ*769(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 192 KMUL
               lxv rA2, POFF+VECSZ*770(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 192 KMUL
               lxv rA3, POFF+VECSZ*771(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 192 KMUL
               lxv rA4, POFF+VECSZ*772(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 192 KMUL
               lxv rA5, POFF+VECSZ*773(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 192 KMUL
               lxv rA6, POFF+VECSZ*774(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 192 KMUL
               lxv rA7, POFF+VECSZ*775(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 194 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*581(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 194 KMUL
               lxv rB0, POFF+VECSZ*582(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 194 KMUL
               lxv rB1, POFF+VECSZ*583(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 194 KMUL
               lxv rB2, POFF+VECSZ*584(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 194 KMUL
               lxv rB3, POFF+VECSZ*585(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 194 KMUL
               lxv rB4, POFF+VECSZ*586(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 194 KMUL
               lxv rA0, POFF+VECSZ*776(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 194 KMUL
               lxv rA1, POFF+VECSZ*777(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 194 KMUL
               lxv rA2, POFF+VECSZ*778(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 194 KMUL
               lxv rA3, POFF+VECSZ*779(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 194 KMUL
               lxv rA4, POFF+VECSZ*780(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 194 KMUL
               lxv rA5, POFF+VECSZ*781(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 194 KMUL
               lxv rA6, POFF+VECSZ*782(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 194 KMUL
               lxv rA7, POFF+VECSZ*783(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 196 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*587(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 196 KMUL
               lxv rB0, POFF+VECSZ*588(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 196 KMUL
               lxv rB1, POFF+VECSZ*589(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 196 KMUL
               lxv rB2, POFF+VECSZ*590(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 196 KMUL
               lxv rB3, POFF+VECSZ*591(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 196 KMUL
               lxv rB4, POFF+VECSZ*592(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 196 KMUL
               lxv rA0, POFF+VECSZ*784(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 196 KMUL
               lxv rA1, POFF+VECSZ*785(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 196 KMUL
               lxv rA2, POFF+VECSZ*786(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 196 KMUL
               lxv rA3, POFF+VECSZ*787(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 196 KMUL
               lxv rA4, POFF+VECSZ*788(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 196 KMUL
               lxv rA5, POFF+VECSZ*789(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 196 KMUL
               lxv rA6, POFF+VECSZ*790(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 196 KMUL
               lxv rA7, POFF+VECSZ*791(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 198 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*593(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 198 KMUL
               lxv rB0, POFF+VECSZ*594(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 198 KMUL
               lxv rB1, POFF+VECSZ*595(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 198 KMUL
               lxv rB2, POFF+VECSZ*596(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 198 KMUL
               lxv rB3, POFF+VECSZ*597(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 198 KMUL
               lxv rB4, POFF+VECSZ*598(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 198 KMUL
               lxv rA0, POFF+VECSZ*792(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 198 KMUL
               lxv rA1, POFF+VECSZ*793(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 198 KMUL
               lxv rA2, POFF+VECSZ*794(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 198 KMUL
               lxv rA3, POFF+VECSZ*795(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 198 KMUL
               lxv rA4, POFF+VECSZ*796(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 198 KMUL
               lxv rA5, POFF+VECSZ*797(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 198 KMUL
               lxv rA6, POFF+VECSZ*798(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 198 KMUL
               lxv rA7, POFF+VECSZ*799(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 200 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*599(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 200 KMUL
               lxv rB0, POFF+VECSZ*600(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 200 KMUL
               lxv rB1, POFF+VECSZ*601(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 200 KMUL
               lxv rB2, POFF+VECSZ*602(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 200 KMUL
               lxv rB3, POFF+VECSZ*603(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 200 KMUL
               lxv rB4, POFF+VECSZ*604(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 200 KMUL
               lxv rA0, POFF+VECSZ*800(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 200 KMUL
               lxv rA1, POFF+VECSZ*801(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 200 KMUL
               lxv rA2, POFF+VECSZ*802(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 200 KMUL
               lxv rA3, POFF+VECSZ*803(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 200 KMUL
               lxv rA4, POFF+VECSZ*804(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 200 KMUL
               lxv rA5, POFF+VECSZ*805(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 200 KMUL
               lxv rA6, POFF+VECSZ*806(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 200 KMUL
               lxv rA7, POFF+VECSZ*807(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 202 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*605(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 202 KMUL
               lxv rB0, POFF+VECSZ*606(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 202 KMUL
               lxv rB1, POFF+VECSZ*607(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 202 KMUL
               lxv rB2, POFF+VECSZ*608(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 202 KMUL
               lxv rB3, POFF+VECSZ*609(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 202 KMUL
               lxv rB4, POFF+VECSZ*610(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 202 KMUL
               lxv rA0, POFF+VECSZ*808(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 202 KMUL
               lxv rA1, POFF+VECSZ*809(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 202 KMUL
               lxv rA2, POFF+VECSZ*810(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 202 KMUL
               lxv rA3, POFF+VECSZ*811(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 202 KMUL
               lxv rA4, POFF+VECSZ*812(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 202 KMUL
               lxv rA5, POFF+VECSZ*813(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 202 KMUL
               lxv rA6, POFF+VECSZ*814(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 202 KMUL
               lxv rA7, POFF+VECSZ*815(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 204 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*611(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 204 KMUL
               lxv rB0, POFF+VECSZ*612(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 204 KMUL
               lxv rB1, POFF+VECSZ*613(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 204 KMUL
               lxv rB2, POFF+VECSZ*614(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 204 KMUL
               lxv rB3, POFF+VECSZ*615(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 204 KMUL
               lxv rB4, POFF+VECSZ*616(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 204 KMUL
               lxv rA0, POFF+VECSZ*816(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 204 KMUL
               lxv rA1, POFF+VECSZ*817(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 204 KMUL
               lxv rA2, POFF+VECSZ*818(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 204 KMUL
               lxv rA3, POFF+VECSZ*819(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 204 KMUL
               lxv rA4, POFF+VECSZ*820(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 204 KMUL
               lxv rA5, POFF+VECSZ*821(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 204 KMUL
               lxv rA6, POFF+VECSZ*822(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 204 KMUL
               lxv rA7, POFF+VECSZ*823(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 206 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*617(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 206 KMUL
               lxv rB0, POFF+VECSZ*618(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 206 KMUL
               lxv rB1, POFF+VECSZ*619(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 206 KMUL
               lxv rB2, POFF+VECSZ*620(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 206 KMUL
               lxv rB3, POFF+VECSZ*621(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 206 KMUL
               lxv rB4, POFF+VECSZ*622(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 206 KMUL
               lxv rA0, POFF+VECSZ*824(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 206 KMUL
               lxv rA1, POFF+VECSZ*825(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 206 KMUL
               lxv rA2, POFF+VECSZ*826(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 206 KMUL
               lxv rA3, POFF+VECSZ*827(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 206 KMUL
               lxv rA4, POFF+VECSZ*828(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 206 KMUL
               lxv rA5, POFF+VECSZ*829(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 206 KMUL
               lxv rA6, POFF+VECSZ*830(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 206 KMUL
               lxv rA7, POFF+VECSZ*831(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 208 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*623(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 208 KMUL
               lxv rB0, POFF+VECSZ*624(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 208 KMUL
               lxv rB1, POFF+VECSZ*625(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 208 KMUL
               lxv rB2, POFF+VECSZ*626(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 208 KMUL
               lxv rB3, POFF+VECSZ*627(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 208 KMUL
               lxv rB4, POFF+VECSZ*628(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 208 KMUL
               lxv rA0, POFF+VECSZ*832(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 208 KMUL
               lxv rA1, POFF+VECSZ*833(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 208 KMUL
               lxv rA2, POFF+VECSZ*834(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 208 KMUL
               lxv rA3, POFF+VECSZ*835(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 208 KMUL
               lxv rA4, POFF+VECSZ*836(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 208 KMUL
               lxv rA5, POFF+VECSZ*837(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 208 KMUL
               lxv rA6, POFF+VECSZ*838(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 208 KMUL
               lxv rA7, POFF+VECSZ*839(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 210 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*629(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 210 KMUL
               lxv rB0, POFF+VECSZ*630(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 210 KMUL
               lxv rB1, POFF+VECSZ*631(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 210 KMUL
               lxv rB2, POFF+VECSZ*632(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 210 KMUL
               lxv rB3, POFF+VECSZ*633(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 210 KMUL
               lxv rB4, POFF+VECSZ*634(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 210 KMUL
               lxv rA0, POFF+VECSZ*840(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 210 KMUL
               lxv rA1, POFF+VECSZ*841(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 210 KMUL
               lxv rA2, POFF+VECSZ*842(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 210 KMUL
               lxv rA3, POFF+VECSZ*843(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 210 KMUL
               lxv rA4, POFF+VECSZ*844(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 210 KMUL
               lxv rA5, POFF+VECSZ*845(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 210 KMUL
               lxv rA6, POFF+VECSZ*846(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 210 KMUL
               lxv rA7, POFF+VECSZ*847(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 212 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*635(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 212 KMUL
               lxv rB0, POFF+VECSZ*636(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 212 KMUL
               lxv rB1, POFF+VECSZ*637(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 212 KMUL
               lxv rB2, POFF+VECSZ*638(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 212 KMUL
               lxv rB3, POFF+VECSZ*639(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 212 KMUL
               lxv rB4, POFF+VECSZ*640(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 212 KMUL
               lxv rA0, POFF+VECSZ*848(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 212 KMUL
               lxv rA1, POFF+VECSZ*849(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 212 KMUL
               lxv rA2, POFF+VECSZ*850(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 212 KMUL
               lxv rA3, POFF+VECSZ*851(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 212 KMUL
               lxv rA4, POFF+VECSZ*852(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 212 KMUL
               lxv rA5, POFF+VECSZ*853(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 212 KMUL
               lxv rA6, POFF+VECSZ*854(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 212 KMUL
               lxv rA7, POFF+VECSZ*855(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 214 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*641(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 214 KMUL
               lxv rB0, POFF+VECSZ*642(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 214 KMUL
               lxv rB1, POFF+VECSZ*643(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 214 KMUL
               lxv rB2, POFF+VECSZ*644(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 214 KMUL
               lxv rB3, POFF+VECSZ*645(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 214 KMUL
               lxv rB4, POFF+VECSZ*646(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 214 KMUL
               lxv rA0, POFF+VECSZ*856(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 214 KMUL
               lxv rA1, POFF+VECSZ*857(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 214 KMUL
               lxv rA2, POFF+VECSZ*858(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 214 KMUL
               lxv rA3, POFF+VECSZ*859(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 214 KMUL
               lxv rA4, POFF+VECSZ*860(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 214 KMUL
               lxv rA5, POFF+VECSZ*861(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 214 KMUL
               lxv rA6, POFF+VECSZ*862(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 214 KMUL
               lxv rA7, POFF+VECSZ*863(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 216 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*647(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 216 KMUL
               lxv rB0, POFF+VECSZ*648(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 216 KMUL
               lxv rB1, POFF+VECSZ*649(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 216 KMUL
               lxv rB2, POFF+VECSZ*650(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 216 KMUL
               lxv rB3, POFF+VECSZ*651(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 216 KMUL
               lxv rB4, POFF+VECSZ*652(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 216 KMUL
               lxv rA0, POFF+VECSZ*864(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 216 KMUL
               lxv rA1, POFF+VECSZ*865(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 216 KMUL
               lxv rA2, POFF+VECSZ*866(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 216 KMUL
               lxv rA3, POFF+VECSZ*867(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 216 KMUL
               lxv rA4, POFF+VECSZ*868(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 216 KMUL
               lxv rA5, POFF+VECSZ*869(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 216 KMUL
               lxv rA6, POFF+VECSZ*870(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 216 KMUL
               lxv rA7, POFF+VECSZ*871(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 218 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*653(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 218 KMUL
               lxv rB0, POFF+VECSZ*654(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 218 KMUL
               lxv rB1, POFF+VECSZ*655(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 218 KMUL
               lxv rB2, POFF+VECSZ*656(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 218 KMUL
               lxv rB3, POFF+VECSZ*657(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 218 KMUL
               lxv rB4, POFF+VECSZ*658(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 218 KMUL
               lxv rA0, POFF+VECSZ*872(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 218 KMUL
               lxv rA1, POFF+VECSZ*873(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 218 KMUL
               lxv rA2, POFF+VECSZ*874(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 218 KMUL
               lxv rA3, POFF+VECSZ*875(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 218 KMUL
               lxv rA4, POFF+VECSZ*876(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 218 KMUL
               lxv rA5, POFF+VECSZ*877(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 218 KMUL
               lxv rA6, POFF+VECSZ*878(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 218 KMUL
               lxv rA7, POFF+VECSZ*879(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 220 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*659(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 220 KMUL
               lxv rB0, POFF+VECSZ*660(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 220 KMUL
               lxv rB1, POFF+VECSZ*661(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 220 KMUL
               lxv rB2, POFF+VECSZ*662(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 220 KMUL
               lxv rB3, POFF+VECSZ*663(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 220 KMUL
               lxv rB4, POFF+VECSZ*664(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 220 KMUL
               lxv rA0, POFF+VECSZ*880(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 220 KMUL
               lxv rA1, POFF+VECSZ*881(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 220 KMUL
               lxv rA2, POFF+VECSZ*882(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 220 KMUL
               lxv rA3, POFF+VECSZ*883(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 220 KMUL
               lxv rA4, POFF+VECSZ*884(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 220 KMUL
               lxv rA5, POFF+VECSZ*885(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 220 KMUL
               lxv rA6, POFF+VECSZ*886(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 220 KMUL
               lxv rA7, POFF+VECSZ*887(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 222 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*665(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 222 KMUL
               lxv rB0, POFF+VECSZ*666(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 222 KMUL
               lxv rB1, POFF+VECSZ*667(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 222 KMUL
               lxv rB2, POFF+VECSZ*668(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 222 KMUL
               lxv rB3, POFF+VECSZ*669(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 222 KMUL
               lxv rB4, POFF+VECSZ*670(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 222 KMUL
               lxv rA0, POFF+VECSZ*888(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 222 KMUL
               lxv rA1, POFF+VECSZ*889(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 222 KMUL
               lxv rA2, POFF+VECSZ*890(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 222 KMUL
               lxv rA3, POFF+VECSZ*891(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 222 KMUL
               lxv rA4, POFF+VECSZ*892(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 222 KMUL
               lxv rA5, POFF+VECSZ*893(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 222 KMUL
               lxv rA6, POFF+VECSZ*894(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 222 KMUL
               lxv rA7, POFF+VECSZ*895(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 224 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*671(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 224 KMUL
               lxv rB0, POFF+VECSZ*672(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 224 KMUL
               lxv rB1, POFF+VECSZ*673(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 224 KMUL
               lxv rB2, POFF+VECSZ*674(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 224 KMUL
               lxv rB3, POFF+VECSZ*675(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 224 KMUL
               lxv rB4, POFF+VECSZ*676(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 224 KMUL
               lxv rA0, POFF+VECSZ*896(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 224 KMUL
               lxv rA1, POFF+VECSZ*897(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 224 KMUL
               lxv rA2, POFF+VECSZ*898(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 224 KMUL
               lxv rA3, POFF+VECSZ*899(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 224 KMUL
               lxv rA4, POFF+VECSZ*900(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 224 KMUL
               lxv rA5, POFF+VECSZ*901(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 224 KMUL
               lxv rA6, POFF+VECSZ*902(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 224 KMUL
               lxv rA7, POFF+VECSZ*903(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 226 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*677(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 226 KMUL
               lxv rB0, POFF+VECSZ*678(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 226 KMUL
               lxv rB1, POFF+VECSZ*679(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 226 KMUL
               lxv rB2, POFF+VECSZ*680(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 226 KMUL
               lxv rB3, POFF+VECSZ*681(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 226 KMUL
               lxv rB4, POFF+VECSZ*682(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 226 KMUL
               lxv rA0, POFF+VECSZ*904(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 226 KMUL
               lxv rA1, POFF+VECSZ*905(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 226 KMUL
               lxv rA2, POFF+VECSZ*906(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 226 KMUL
               lxv rA3, POFF+VECSZ*907(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 226 KMUL
               lxv rA4, POFF+VECSZ*908(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 226 KMUL
               lxv rA5, POFF+VECSZ*909(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 226 KMUL
               lxv rA6, POFF+VECSZ*910(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 226 KMUL
               lxv rA7, POFF+VECSZ*911(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 228 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*683(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 228 KMUL
               lxv rB0, POFF+VECSZ*684(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 228 KMUL
               lxv rB1, POFF+VECSZ*685(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 228 KMUL
               lxv rB2, POFF+VECSZ*686(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 228 KMUL
               lxv rB3, POFF+VECSZ*687(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 228 KMUL
               lxv rB4, POFF+VECSZ*688(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 228 KMUL
               lxv rA0, POFF+VECSZ*912(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 228 KMUL
               lxv rA1, POFF+VECSZ*913(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 228 KMUL
               lxv rA2, POFF+VECSZ*914(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 228 KMUL
               lxv rA3, POFF+VECSZ*915(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 228 KMUL
               lxv rA4, POFF+VECSZ*916(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 228 KMUL
               lxv rA5, POFF+VECSZ*917(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 228 KMUL
               lxv rA6, POFF+VECSZ*918(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 228 KMUL
               lxv rA7, POFF+VECSZ*919(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 230 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*689(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 230 KMUL
               lxv rB0, POFF+VECSZ*690(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 230 KMUL
               lxv rB1, POFF+VECSZ*691(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 230 KMUL
               lxv rB2, POFF+VECSZ*692(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 230 KMUL
               lxv rB3, POFF+VECSZ*693(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 230 KMUL
               lxv rB4, POFF+VECSZ*694(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 230 KMUL
               lxv rA0, POFF+VECSZ*920(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 230 KMUL
               lxv rA1, POFF+VECSZ*921(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 230 KMUL
               lxv rA2, POFF+VECSZ*922(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 230 KMUL
               lxv rA3, POFF+VECSZ*923(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 230 KMUL
               lxv rA4, POFF+VECSZ*924(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 230 KMUL
               lxv rA5, POFF+VECSZ*925(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 230 KMUL
               lxv rA6, POFF+VECSZ*926(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 230 KMUL
               lxv rA7, POFF+VECSZ*927(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 232 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*695(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 232 KMUL
               lxv rB0, POFF+VECSZ*696(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 232 KMUL
               lxv rB1, POFF+VECSZ*697(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 232 KMUL
               lxv rB2, POFF+VECSZ*698(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 232 KMUL
               lxv rB3, POFF+VECSZ*699(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 232 KMUL
               lxv rB4, POFF+VECSZ*700(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 232 KMUL
               lxv rA0, POFF+VECSZ*928(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 232 KMUL
               lxv rA1, POFF+VECSZ*929(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 232 KMUL
               lxv rA2, POFF+VECSZ*930(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 232 KMUL
               lxv rA3, POFF+VECSZ*931(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 232 KMUL
               lxv rA4, POFF+VECSZ*932(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 232 KMUL
               lxv rA5, POFF+VECSZ*933(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 232 KMUL
               lxv rA6, POFF+VECSZ*934(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 232 KMUL
               lxv rA7, POFF+VECSZ*935(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 234 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*701(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 234 KMUL
               lxv rB0, POFF+VECSZ*702(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 234 KMUL
               lxv rB1, POFF+VECSZ*703(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 234 KMUL
               lxv rB2, POFF+VECSZ*704(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 234 KMUL
               lxv rB3, POFF+VECSZ*705(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 234 KMUL
               lxv rB4, POFF+VECSZ*706(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 234 KMUL
               lxv rA0, POFF+VECSZ*936(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 234 KMUL
               lxv rA1, POFF+VECSZ*937(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 234 KMUL
               lxv rA2, POFF+VECSZ*938(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 234 KMUL
               lxv rA3, POFF+VECSZ*939(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 234 KMUL
               lxv rA4, POFF+VECSZ*940(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 234 KMUL
               lxv rA5, POFF+VECSZ*941(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 234 KMUL
               lxv rA6, POFF+VECSZ*942(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 234 KMUL
               lxv rA7, POFF+VECSZ*943(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 236 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*707(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 236 KMUL
               lxv rB0, POFF+VECSZ*708(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 236 KMUL
               lxv rB1, POFF+VECSZ*709(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 236 KMUL
               lxv rB2, POFF+VECSZ*710(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 236 KMUL
               lxv rB3, POFF+VECSZ*711(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 236 KMUL
               lxv rB4, POFF+VECSZ*712(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 236 KMUL
               lxv rA0, POFF+VECSZ*944(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 236 KMUL
               lxv rA1, POFF+VECSZ*945(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 236 KMUL
               lxv rA2, POFF+VECSZ*946(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 236 KMUL
               lxv rA3, POFF+VECSZ*947(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 236 KMUL
               lxv rA4, POFF+VECSZ*948(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 236 KMUL
               lxv rA5, POFF+VECSZ*949(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 236 KMUL
               lxv rA6, POFF+VECSZ*950(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 236 KMUL
               lxv rA7, POFF+VECSZ*951(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 238 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*713(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 238 KMUL
               lxv rB0, POFF+VECSZ*714(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 238 KMUL
               lxv rB1, POFF+VECSZ*715(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 238 KMUL
               lxv rB2, POFF+VECSZ*716(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 238 KMUL
               lxv rB3, POFF+VECSZ*717(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 238 KMUL
               lxv rB4, POFF+VECSZ*718(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 238 KMUL
               lxv rA0, POFF+VECSZ*952(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 238 KMUL
               lxv rA1, POFF+VECSZ*953(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 238 KMUL
               lxv rA2, POFF+VECSZ*954(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 238 KMUL
               lxv rA3, POFF+VECSZ*955(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 238 KMUL
               lxv rA4, POFF+VECSZ*956(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 238 KMUL
               lxv rA5, POFF+VECSZ*957(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 238 KMUL
               lxv rA6, POFF+VECSZ*958(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 238 KMUL
               lxv rA7, POFF+VECSZ*959(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 240 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*719(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 240 KMUL
               lxv rB0, POFF+VECSZ*720(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 240 KMUL
               lxv rB1, POFF+VECSZ*721(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 240 KMUL
               lxv rB2, POFF+VECSZ*722(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 240 KMUL
               lxv rB3, POFF+VECSZ*723(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 240 KMUL
               lxv rB4, POFF+VECSZ*724(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 240 KMUL
               lxv rA0, POFF+VECSZ*960(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 240 KMUL
               lxv rA1, POFF+VECSZ*961(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 240 KMUL
               lxv rA2, POFF+VECSZ*962(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 240 KMUL
               lxv rA3, POFF+VECSZ*963(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 240 KMUL
               lxv rA4, POFF+VECSZ*964(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 240 KMUL
               lxv rA5, POFF+VECSZ*965(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 240 KMUL
               lxv rA6, POFF+VECSZ*966(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 240 KMUL
               lxv rA7, POFF+VECSZ*967(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 242 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*725(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 242 KMUL
               lxv rB0, POFF+VECSZ*726(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 242 KMUL
               lxv rB1, POFF+VECSZ*727(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 242 KMUL
               lxv rB2, POFF+VECSZ*728(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 242 KMUL
               lxv rB3, POFF+VECSZ*729(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 242 KMUL
               lxv rB4, POFF+VECSZ*730(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 242 KMUL
               lxv rA0, POFF+VECSZ*968(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 242 KMUL
               lxv rA1, POFF+VECSZ*969(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 242 KMUL
               lxv rA2, POFF+VECSZ*970(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 242 KMUL
               lxv rA3, POFF+VECSZ*971(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 242 KMUL
               lxv rA4, POFF+VECSZ*972(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 242 KMUL
               lxv rA5, POFF+VECSZ*973(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 242 KMUL
               lxv rA6, POFF+VECSZ*974(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 242 KMUL
               lxv rA7, POFF+VECSZ*975(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 244 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*731(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 244 KMUL
               lxv rB0, POFF+VECSZ*732(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 244 KMUL
               lxv rB1, POFF+VECSZ*733(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 244 KMUL
               lxv rB2, POFF+VECSZ*734(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 244 KMUL
               lxv rB3, POFF+VECSZ*735(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 244 KMUL
               lxv rB4, POFF+VECSZ*736(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 244 KMUL
               lxv rA0, POFF+VECSZ*976(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 244 KMUL
               lxv rA1, POFF+VECSZ*977(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 244 KMUL
               lxv rA2, POFF+VECSZ*978(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 244 KMUL
               lxv rA3, POFF+VECSZ*979(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 244 KMUL
               lxv rA4, POFF+VECSZ*980(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 244 KMUL
               lxv rA5, POFF+VECSZ*981(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 244 KMUL
               lxv rA6, POFF+VECSZ*982(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 244 KMUL
               lxv rA7, POFF+VECSZ*983(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 246 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*737(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 246 KMUL
               lxv rB0, POFF+VECSZ*738(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 246 KMUL
               lxv rB1, POFF+VECSZ*739(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 246 KMUL
               lxv rB2, POFF+VECSZ*740(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 246 KMUL
               lxv rB3, POFF+VECSZ*741(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 246 KMUL
               lxv rB4, POFF+VECSZ*742(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 246 KMUL
               lxv rA0, POFF+VECSZ*984(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 246 KMUL
               lxv rA1, POFF+VECSZ*985(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 246 KMUL
               lxv rA2, POFF+VECSZ*986(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 246 KMUL
               lxv rA3, POFF+VECSZ*987(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 246 KMUL
               lxv rA4, POFF+VECSZ*988(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 246 KMUL
               lxv rA5, POFF+VECSZ*989(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 246 KMUL
               lxv rA6, POFF+VECSZ*990(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 246 KMUL
               lxv rA7, POFF+VECSZ*991(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 248 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*743(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 248 KMUL
               lxv rB0, POFF+VECSZ*744(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 248 KMUL
               lxv rB1, POFF+VECSZ*745(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 248 KMUL
               lxv rB2, POFF+VECSZ*746(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 248 KMUL
               lxv rB3, POFF+VECSZ*747(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 248 KMUL
               lxv rB4, POFF+VECSZ*748(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 248 KMUL
               lxv rA0, POFF+VECSZ*992(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 248 KMUL
               lxv rA1, POFF+VECSZ*993(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 248 KMUL
               lxv rA2, POFF+VECSZ*994(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 248 KMUL
               lxv rA3, POFF+VECSZ*995(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 248 KMUL
               lxv rA4, POFF+VECSZ*996(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 248 KMUL
               lxv rA5, POFF+VECSZ*997(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 248 KMUL
               lxv rA6, POFF+VECSZ*998(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 248 KMUL
               lxv rA7, POFF+VECSZ*999(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 250 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*749(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 250 KMUL
               lxv rB0, POFF+VECSZ*750(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 250 KMUL
               lxv rB1, POFF+VECSZ*751(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 250 KMUL
               lxv rB2, POFF+VECSZ*752(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 250 KMUL
               lxv rB3, POFF+VECSZ*753(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 250 KMUL
               lxv rB4, POFF+VECSZ*754(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 250 KMUL
               lxv rA0, POFF+VECSZ*1000(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 250 KMUL
               lxv rA1, POFF+VECSZ*1001(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 250 KMUL
               lxv rA2, POFF+VECSZ*1002(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 250 KMUL
               lxv rA3, POFF+VECSZ*1003(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 250 KMUL
               lxv rA4, POFF+VECSZ*1004(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 250 KMUL
               lxv rA5, POFF+VECSZ*1005(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 250 KMUL
               lxv rA6, POFF+VECSZ*1006(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 250 KMUL
               lxv rA7, POFF+VECSZ*1007(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 252 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*755(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 252 KMUL
               lxv rB0, POFF+VECSZ*756(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 252 KMUL
               lxv rB1, POFF+VECSZ*757(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 252 KMUL
               lxv rB2, POFF+VECSZ*758(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 252 KMUL
               lxv rB3, POFF+VECSZ*759(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 252 KMUL
               lxv rB4, POFF+VECSZ*760(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 252 KMUL
               lxv rA0, POFF+VECSZ*1008(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 252 KMUL
               lxv rA1, POFF+VECSZ*1009(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 252 KMUL
               lxv rA2, POFF+VECSZ*1010(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 252 KMUL
               lxv rA3, POFF+VECSZ*1011(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 252 KMUL
               lxv rA4, POFF+VECSZ*1012(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 252 KMUL
               lxv rA5, POFF+VECSZ*1013(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 252 KMUL
               lxv rA6, POFF+VECSZ*1014(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 252 KMUL
               lxv rA7, POFF+VECSZ*1015(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 254 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*761(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 254 KMUL
               lxv rB0, POFF+VECSZ*762(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 254 KMUL
               lxv rB1, POFF+VECSZ*763(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 254 KMUL
               lxv rB2, POFF+VECSZ*764(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 254 KMUL
               lxv rB3, POFF+VECSZ*765(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 254 KMUL
               lxv rB4, POFF+VECSZ*766(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 254 KMUL
               lxv rA0, POFF+VECSZ*1016(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 254 KMUL
               lxv rA1, POFF+VECSZ*1017(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 254 KMUL
               lxv rA2, POFF+VECSZ*1018(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 254 KMUL
               lxv rA3, POFF+VECSZ*1019(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 254 KMUL
               lxv rA4, POFF+VECSZ*1020(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 254 KMUL
               lxv rA5, POFF+VECSZ*1021(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 254 KMUL
               lxv rA6, POFF+VECSZ*1022(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 254 KMUL
               lxv rA7, POFF+VECSZ*1023(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 256 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*767(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 256 KMUL
               lxv rB0, POFF+VECSZ*768(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 256 KMUL
               lxv rB1, POFF+VECSZ*769(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 256 KMUL
               lxv rB2, POFF+VECSZ*770(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 256 KMUL
               lxv rB3, POFF+VECSZ*771(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 256 KMUL
               lxv rB4, POFF+VECSZ*772(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 256 KMUL
               lxv rA0, POFF+VECSZ*1024(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 256 KMUL
               lxv rA1, POFF+VECSZ*1025(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 256 KMUL
               lxv rA2, POFF+VECSZ*1026(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 256 KMUL
               lxv rA3, POFF+VECSZ*1027(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 256 KMUL
               lxv rA4, POFF+VECSZ*1028(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 256 KMUL
               lxv rA5, POFF+VECSZ*1029(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 256 KMUL
               lxv rA6, POFF+VECSZ*1030(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 256 KMUL
               lxv rA7, POFF+VECSZ*1031(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
         #if KB >= 258 KMUL
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*773(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
            #if KB > 258 KMUL
               lxv rB0, POFF+VECSZ*774(pB)
            #elif !defined(BETA0)
               lxv rB0, 0(pC)
            #endif
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
            #if KB > 258 KMUL
               lxv rB1, POFF+VECSZ*775(pB)
            #elif !defined(BETA0)
               lxv rB1, 16(pC)
            #endif
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
            #if KB > 258 KMUL
               lxv rB2, POFF+VECSZ*776(pB)
            #elif !defined(BETA0)
               lxv rB2, 32(pC)
            #endif
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
            #if KB > 258 KMUL
               lxv rB3, POFF+VECSZ*777(pB)
            #elif !defined(BETA0)
               lxv rB3, 48(pC)
            #endif
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
            #if KB > 258 KMUL
               lxv rB4, POFF+VECSZ*778(pB)
            #elif !defined(BETA0)
               lxv rB4, 64(pC)
            #endif
            xvmaddadp rC05, rA0, rB5
            #if KB > 258 KMUL
               lxv rA0, POFF+VECSZ*1032(pA)
            #elif !defined(BETA0)
               lxv rA0, 80(pC)
            #endif
            xvmaddadp rC15, rA1, rB5
            #if KB > 258 KMUL
               lxv rA1, POFF+VECSZ*1033(pA)
            #elif !defined(BETA0)
               lxv rA1, 96(pC)
            #endif
            xvmaddadp rC25, rA2, rB5
            #if KB > 258 KMUL
               lxv rA2, POFF+VECSZ*1034(pA)
            #elif !defined(BETA0)
               lxv rA2, 112(pC)
            #endif
            xvmaddadp rC35, rA3, rB5
            #if KB > 258 KMUL
               lxv rA3, POFF+VECSZ*1035(pA)
            #elif !defined(BETA0)
               lxv rA3, 128(pC)
            #endif
            xvmaddadp rC45, rA4, rB5
            #if KB > 258 KMUL
               lxv rA4, POFF+VECSZ*1036(pA)
            #elif !defined(BETA0)
               lxv rA4, 144(pC)
            #endif
            xvmaddadp rC55, rA5, rB5
            #if KB > 258 KMUL
               lxv rA5, POFF+VECSZ*1037(pA)
            #elif !defined(BETA0)
               lxv rA5, 160(pC)
            #endif
            xvmaddadp rC65, rA6, rB5
            #if KB > 258 KMUL
               lxv rA6, POFF+VECSZ*1038(pA)
            #elif !defined(BETA0)
               lxv rA6, 176(pC)
            #endif
            xvmaddadp rC75, rA7, rB5
            #if KB > 258 KMUL
               lxv rA7, POFF+VECSZ*1039(pA)
            #elif !defined(BETA0)
               lxv rA7, 192(pC)
            #endif
         #endif
/*
 *       Write out to C, already loaded: rB0-4=C[0-9], rA0=c[10-11],
 *          rA1=c[12-13], rA2=c[14-15], rA3=c[16-17], rA4=c[18-19],
 *          rA5=c[20-21], rA6=c[22-23]
 */
         #if KB == 2  /* if no unroll, no preload of C! */
            lxv rB0, 0*16(pC)
            lxv rB1, 1*16(pC)
            lxv rB2, 2*16(pC)
            lxv rB3, 3*16(pC)
            lxv rB4, 4*16(pC)
            lxv rA0, 5*16(pC)
            lxv rA1, 6*16(pC)
            lxv rA2, 7*16(pC)
            lxv rA3, 8*16(pC)
            lxv rA4, 9*16(pC)
            lxv rA5, 10*16(pC)
            lxv rA6, 11*16(pC)
         #endif

         #ifdef SREAL
            ATL_vvrsum4(rC00, rC10, rC20, rC30, rB5, rA7);
            #ifndef BETA0
               VADDC rC00, rC00, rB0
            #endif
            ATL_vvrsum4(rC40, rC50, rC60, rC70, rB5, rA7);
            #ifndef BETA0
               VADDC rC40, rC40, rB1
            #endif
            ATL_vvrsum4(rC01, rC11, rC21, rC31, rB5, rA7);
            #ifndef BETA0
               VADDC rC01, rC01, rB2
            #endif
            ATL_vvrsum4(rC41, rC51, rC61, rC71, rB5, rA7);
            #ifndef BETA0
               VADDC rC41, rC41, rB3
            #endif
            stxv rC00, 16*0(pC)
            stxv rC40, 16*1(pC)
            stxv rC01, 16*2(pC)
            stxv rC41, 16*3(pC)
            ATL_vvrsum4(rC02, rC12, rC22, rC32, rB5, rA7);
            #ifndef BETA0
               VADDC rC02, rC02, rB4
            #endif
            ATL_vvrsum4(rC42, rC52, rC62, rC72, rB5, rA7);
            #ifndef BETA0
               VADDC rC42, rC42, rA0
            #endif
            ATL_vvrsum4(rC03, rC13, rC23, rC33, rB5, rA7);
            #ifndef BETA0
               VADDC rC03, rC03, rA1
            #endif
            ATL_vvrsum4(rC43, rC53, rC63, rC73, rB5, rA7);
            #ifndef BETA0
               VADDC rC43, rC43, rA2
            #endif
            stxv rC02, 16*4(pC)
            stxv rC42, 16*5(pC)
            stxv rC03, 16*6(pC)
            stxv rC43, 16*7(pC)
            ATL_vvrsum4(rC04, rC14, rC24, rC34, rB5, rA7);
            #ifndef BETA0
               VADDC rC04, rC04, rA3
            #endif
            ATL_vvrsum4(rC44, rC54, rC64, rC74, rB5, rA7);
            #ifndef BETA0
               VADDC rC44, rC44, rA4
            #endif
            ATL_vvrsum4(rC05, rC15, rC25, rC35, rB5, rA7);
            #ifndef BETA0
               VADDC rC05, rC05, rA5
            #endif
            ATL_vvrsum4(rC45, rC55, rC65, rC75, rB5, rA7);
            #ifndef BETA0
               VADDC rC45, rC45, rA6
            #endif
            stxv rC04, 16*8(pC)
            stxv rC44, 16*9(pC)
            stxv rC05, 16*10(pC)
            stxv rC45, 16*11(pC)
         #else
            ATL_vvrsum2(rC00, rC10, rB5);
            #ifndef BETA0
               VADDC rC00, rC00, rB0
            #endif
            #ifndef BETA0
               lxv rC10, 24*8(pC)   /* rC10 = c[24,25] */
               lxv rB0, 26*8(pC)    /* rB0  = c[26,27] */
            #endif
            ATL_vvrsum2(rC20, rC30, rB5);
            stxv rC00, 8*0(pC)
            #ifndef BETA0
               VADDC rC20, rC20, rB1
            #endif
            #ifndef BETA0
               lxv rC30, 28*8(pC)   /* rC30 = c[28,29] */
               lxv rB1, 30*8(pC)    /* rB1  = c[30,31] */
            #endif
            ATL_vvrsum2(rC40, rC50, rB5);
            stxv rC20, 8*2(pC)
            #ifndef BETA0
               VADDC rC40, rC40, rB2
            #endif
            #ifndef BETA0
               lxv rC50, 32*8(pC)   /* rC50 = c[32,33] */
               lxv rB2, 34*8(pC)    /* rB2  = c[34,35] */
            #endif
            ATL_vvrsum2(rC60, rC70, rB5);
            stxv rC40, 8*4(pC)
            #ifndef BETA0
               VADDC rC60, rC60, rB3
            #endif
            #ifndef BETA0
               lxv rC70, 36*8(pC)   /* rC70 = c[36,37] */
               lxv rB3,  38*8(pC)   /* rB3  = c[38,39] */
            #endif
            ATL_vvrsum2(rC01, rC11, rB5);
            #ifndef BETA0
               lxv rC00, 40*8(pC)   /* rC00 = c[40,41] */
            #endif
            #ifndef BETA0
               lxv rC20, 42*8(pC)   /* rC20 = c[42,43] */
            #endif
            #ifndef BETA0
               lxv rC40, 44*8(pC)   /* rC40 = c[44,45] */
            #endif
            stxv rC60, 16*3(pC)
            #ifndef BETA0
               VADDC rC01, rC01, rB4
            #endif
            #ifndef BETA0
               lxv rC60, 46*8(pC)   /* rC60 = c[46,47] */
            #endif
            ATL_vvrsum2(rC21, rC31, rB5);
            #ifndef BETA0
               VADDC rC21, rC21, rA0
            #endif
            ATL_vvrsum2(rC41, rC51, rB5);
            #ifndef BETA0
               VADDC rC41, rC41, rA1
            #endif
            ATL_vvrsum2(rC61, rC71, rB5);
            #ifndef BETA0
               VADDC rC61, rC61, rA2
            #endif
            stxv rC01, 8*8(pC)
            stxv rC21, 8*10(pC)
            stxv rC41, 8*12(pC)
            stxv rC61, 8*14(pC)
            ATL_vvrsum2(rC02, rC12, rB5);
            #ifndef BETA0
               VADDC rC02, rC02, rA3
            #endif
            ATL_vvrsum2(rC22, rC32, rB5);
            #ifndef BETA0
               VADDC rC22, rC22, rA4
            #endif
            ATL_vvrsum2(rC42, rC52, rB5);
            #ifndef BETA0
               VADDC rC42, rC42, rA5
            #endif
            ATL_vvrsum2(rC62, rC72, rB5);
            #ifndef BETA0
               VADDC rC62, rC62, rA6
            #endif
            stxv rC02, 8*16(pC)
            stxv rC22, 8*18(pC)
            stxv rC42, 8*20(pC)
            stxv rC62, 8*22(pC)
            ATL_vvrsum2(rC03, rC13, rB5);
            #ifndef BETA0
               VADDC rC03, rC03, rC10
            #endif
            dcbf  0, pC, 3
            ATL_vvrsum2(rC23, rC33, rB5);
            #ifndef BETA0
               VADDC rC23, rC23, rB0
            #endif
            ATL_vvrsum2(rC43, rC53, rB5);
            #ifndef BETA0
               VADDC rC43, rC43, rC30
            #endif
            ATL_vvrsum2(rC63, rC73, rB5);
            #ifndef BETA0
               VADDC rC63, rC63, rB1
            #endif
            stxv rC03, 8*24(pC)
            stxv rC23, 8*26(pC)
            stxv rC43, 8*28(pC)
            stxv rC63, 8*30(pC)
            ATL_vvrsum2(rC04, rC14, rB5);
            #ifndef BETA0
               VADDC rC04, rC04, rC50
            #endif
            ATL_vvrsum2(rC24, rC34, rB5);
            stxv rC04, 8*32(pC)
            #ifndef BETA0
               VADDC rC24, rC24, rB2
            #endif
            ATL_vvrsum2(rC44, rC54, rB5);
            stxv rC24, 8*34(pC)
            #ifndef BETA0
               VADDC rC44, rC44, rC70
            #endif
            ATL_vvrsum2(rC64, rC74, rB5);
            #ifndef BETA0
               VADDC rC64, rC64, rB3
            #endif
            stxv rC44, 8*36(pC)
            stxv rC64, 8*38(pC)
            ATL_vvrsum2(rC05, rC15, rB5);
            dcbf  i128_1, pC, 3
            #ifndef BETA0
               VADDC rC05, rC05, rC00
            #endif
            ATL_vvrsum2(rC25, rC35, rB5);
            stxv rC05, 8*40(pC)
            #ifndef BETA0
               VADDC rC25, rC25, rC20
            #endif
            ATL_vvrsum2(rC45, rC55, rB5);
            stxv rC25, 8*42(pC)
            #ifndef BETA0
               VADDC rC45, rC45, rC40
            #endif
            ATL_vvrsum2(rC65, rC75, rB5);
            stxv rC45, 8*44(pC)
            #ifndef BETA0
               VADDC rC65, rC65, rC60
            #endif
            stxv rC65, 8*46(pC)
            dcbf  i128_2, pC, 3
         #endif
         ADDB pB, pB, incBn     /* pB += KB*nu*sizeof */
         addi pC, pC, 8*6*ELTSZ /* pC += mu*nu*sizeof */
      bdnz+ MLOOP               /* really N-loop */
      mtctr NNU
      mr pB, pB0
      ADDA pA, pA, incAm        /* pA += KB*mu*sizeof */
      addic. NMU, NMU, -1
   bne+ MLOOP
/*
 * Restore vec regs using Altivec instructions
 */
   addi pfA, r1, VOFF   /* pfA points to start of vector save area */
   lvx 20, r0, pfA
      addi pfA, pfA, 16
   lvx 21, r0, pfA
      addi pfA, pfA, 16
   lvx 22, r0, pfA
      addi pfA, pfA, 16
   lvx 23, r0, pfA
      addi pfA, pfA, 16
   lvx 24, r0, pfA
      addi pfA, pfA, 16
   lvx 25, r0, pfA
      addi pfA, pfA, 16
   lvx 26, r0, pfA
      addi pfA, pfA, 16
   lvx 27, r0, pfA
      addi pfA, pfA, 16
   lvx 28, r0, pfA
      addi pfA, pfA, 16
   lvx 29, r0, pfA
      addi pfA, pfA, 16
   lvx 30, r0, pfA
      addi pfA, pfA, 16
   lvx 31, r0, pfA
/*
 * restore f14-31
 */
   lfd f14, FOFF+0*8(r1)
   lfd f15, FOFF+1*8(r1)
   lfd f16, FOFF+2*8(r1)
   lfd f17, FOFF+3*8(r1)
   lfd f18, FOFF+4*8(r1)
   lfd f19, FOFF+5*8(r1)
   lfd f20, FOFF+6*8(r1)
   lfd f21, FOFF+7*8(r1)
   lfd f22, FOFF+8*8(r1)
   lfd f23, FOFF+9*8(r1)
   lfd f24, FOFF+10*8(r1)
   lfd f25, FOFF+11*8(r1)
   lfd f26, FOFF+12*8(r1)
   lfd f27, FOFF+13*8(r1)
   lfd f28, FOFF+14*8(r1)
   lfd f29, FOFF+15*8(r1)
   lfd f30, FOFF+16*8(r1)
   lfd f31, FOFF+17*8(r1)
/*
 * Restore iregs
 */
   ld r2, 24(r1)        /* TOC reg always at 24(sp) */
   lwz r0, 8(r1)        /* CR must be 4-byte saved 8(sp) */
   mtcr r0
   ld r14, IOFF+8*0(r1)
   ld r15, IOFF+8*1(r1)
/*
 * Dealloc stack frame & return
 */
   addi r1, r1, FSIZE
   blr

