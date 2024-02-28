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
#define K0 r5
#define KK r11
#define pA0 r12
#define pBn r14
#define incAm r15
#define incBn r16
#define ADDA add
#define ADDB add
#define i128_0 0
#define i128_1 r17
#define i128_2 r18
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
#define ISZ (5+1)*8 /* +1 keeps 16-byte aligned stack frame */
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
   std r16, IOFF+8*2(r1)
   std r17, IOFF+8*3(r1)
   std r18, IOFF+8*4(r1)
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
   #if 0
      mulli incAm, K0, 8*ELTSZ  /* incAm = K*mu*sizeof */
   #else                        /* incAm = K*8*8 = K*64 */
       #ifdef SREAL
         sldi incAm, K0, 5      /* incAm = K*mu*sizeof = K*8*4 = K*32 */
       #else
         sldi incAm, K0, 6      /* incAm = K*mu*sizeof = K*8*8 = K*64 */
       #endif
   #endif
   #if 0
      mulli incBn, K0, 6*ELTSZ  /* incBn = K*nu*sizeof */
   #else                     /* k*nu*sizeof = K*6*sizeof =K*sizeof*2*3*/
      #ifdef SREAL
         sldi pA0, K0, 3     /* pA0 = K*4*2 = K*sizeof*2 */
      #else
         sldi pA0, K0, 4     /* pA0 = K*16 = K*sizeof*2 */
      #endif
      add incBn, pA0, pA0    /* incBn = K*4*sizeof */
      add incBn, incBn, pA0  /* incBn = K*6*sizeof */
   #endif
   mr pA0, pA
   add pBn, pB, incBn
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
            lxv rB0, POFF+VECSZ*6(pB)
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
            lxv rB1, POFF+VECSZ*7(pB)
         xvmuldp rC02, rA0, rB2
         xvmuldp rC12, rA1, rB2
         xvmuldp rC22, rA2, rB2
         xvmuldp rC32, rA3, rB2
         xvmuldp rC42, rA4, rB2
         xvmuldp rC52, rA5, rB2
         xvmuldp rC62, rA6, rB2
         xvmuldp rC72, rA7, rB2
            lxv rB2, POFF+VECSZ*8(pB)
         xvmuldp rC03, rA0, rB3
         xvmuldp rC13, rA1, rB3
         xvmuldp rC23, rA2, rB3
         xvmuldp rC33, rA3, rB3
         xvmuldp rC43, rA4, rB3
         xvmuldp rC53, rA5, rB3
         xvmuldp rC63, rA6, rB3
         xvmuldp rC73, rA7, rB3
            lxv rB3, POFF+VECSZ*9(pB)
         xvmuldp rC04, rA0, rB4
         xvmuldp rC14, rA1, rB4
         xvmuldp rC24, rA2, rB4
         xvmuldp rC34, rA3, rB4
         xvmuldp rC44, rA4, rB4
         xvmuldp rC54, rA5, rB4
         xvmuldp rC64, rA6, rB4
         xvmuldp rC74, rA7, rB4
            lxv rB4, POFF+VECSZ*10(pB)
         xvmuldp rC05, rA0, rB5
            lxv rA0, POFF+VECSZ*8(pA)
         xvmuldp rC15, rA1, rB5
            lxv rA1, POFF+VECSZ*9(pA)
         xvmuldp rC25, rA2, rB5
            lxv rA2, POFF+VECSZ*10(pA)
         xvmuldp rC35, rA3, rB5
            lxv rA3, POFF+VECSZ*11(pA)
         xvmuldp rC45, rA4, rB5
            lxv rA4, POFF+VECSZ*12(pA)
         xvmuldp rC55, rA5, rB5
            lxv rA5, POFF+VECSZ*13(pA)
         xvmuldp rC65, rA6, rB5
            lxv rA6, POFF+VECSZ*14(pA)
         xvmuldp rC75, rA7, rB5
            lxv rA7, POFF+VECSZ*15(pA)

         addic. KK, K0, -2 KMUL  /* finished VLEN K iterations in peel */
         beq KDONE               /* if (K == 2) goto KDONE; */
/*
 *       This loop unrolled to 16 (32) double (single) iterations, but we
 *       make it work like a ku=VLEN loop by constructing a compiler's
 *       superblock/trace: it has only exists jumping down.
 *       We need to unroll the loop because integer operations interfere
 *       with using the slices for FPU.  This code concentrates all iops
 *       at beginning (computing conditions for trace exit), and end
 *       (loop control & ptr updates).  In this way, cost of switching between
 *       iop/FPU only paid once per 16 (32) iterations.
 */
         KLOOP:
            cmpdi cr2, KK, 4 KMUL
            cmpdi cr3, KK, 8 KMUL
            cmpdi cr4, KK, 12 KMUL
            cmpdi cr5, KK, 16 KMUL
/*
 *          Unroll for K=0,1
 */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*11(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*12(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*13(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*14(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*15(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*16(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*16(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*17(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*18(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*19(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*20(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*21(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*22(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*23(pA)
/*
 *          Unroll for K=2,3
 */
            blt cr2, KDONE /* if (KK == 2) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*17(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*18(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*19(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*20(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*21(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*22(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*24(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*25(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*26(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*27(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*28(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*29(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*30(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*31(pA)
/*
 *          Unroll for K=4,5
 */
            beq cr2, KDONE /* if (KK == 4) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*23(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*24(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*25(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*26(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*27(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*28(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*32(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*33(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*34(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*35(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*36(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*37(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*38(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*39(pA)
/*
 *          Unroll for K=6,7
 */
            blt cr3, KDONE /* if (KK == 6) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*29(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*30(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*31(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*32(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*33(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*34(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*40(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*41(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*42(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*43(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*44(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*45(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*46(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*47(pA)
/*
 *          Unroll for K=8,9
 */
            beq cr3, KDONE /* if (KK == 8) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*35(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*36(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*37(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*38(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*39(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*40(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*48(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*49(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*50(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*51(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*52(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*53(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*54(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*55(pA)
/*
 *          Unroll for K=10,11
 */
            blt cr4, KDONE /* if (KK == 10) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*41(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*42(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*43(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*44(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*45(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*46(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*56(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*57(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*58(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*59(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*60(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*61(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*62(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*63(pA)
/*
 *          Unroll for K=12,13
 */
            beq cr4, KDONE /* if (KK == 12) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*47(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*48(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*49(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*50(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*51(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*52(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*64(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*65(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*66(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*67(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*68(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*69(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*70(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*71(pA)
/*
 *          Unroll for K=14,15
 */
            blt cr5, KDONE /* if (KK == 14) break */
            xvmaddadp rC00, rA0, rB0
               lxv rB5, POFF+VECSZ*53(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
            xvmaddadp rC70, rA7, rB0
               lxv rB0, POFF+VECSZ*54(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
            xvmaddadp rC71, rA7, rB1
               lxv rB1, POFF+VECSZ*55(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
            xvmaddadp rC72, rA7, rB2
               lxv rB2, POFF+VECSZ*56(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
            xvmaddadp rC73, rA7, rB3
               lxv rB3, POFF+VECSZ*57(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
            xvmaddadp rC74, rA7, rB4
               lxv rB4, POFF+VECSZ*58(pB)
            xvmaddadp rC05, rA0, rB5
               lxv rA0, POFF+VECSZ*72(pA)
            xvmaddadp rC15, rA1, rB5
               lxv rA1, POFF+VECSZ*73(pA)
            xvmaddadp rC25, rA2, rB5
               lxv rA2, POFF+VECSZ*74(pA)
            xvmaddadp rC35, rA3, rB5
               lxv rA3, POFF+VECSZ*75(pA)
            xvmaddadp rC45, rA4, rB5
               lxv rA4, POFF+VECSZ*76(pA)
            xvmaddadp rC55, rA5, rB5
               lxv rA5, POFF+VECSZ*77(pA)
            xvmaddadp rC65, rA6, rB5
               lxv rA6, POFF+VECSZ*78(pA)
            xvmaddadp rC75, rA7, rB5
               lxv rA7, POFF+VECSZ*79(pA)
         addi pA, pA, 16*8*8  /* pA += ku*mu*sizeof */
         addi pB, pB, 16*6*8  /* pB += ku*nu*sizeof */
         addi KK, KK, -16 KMUL
      bne cr5, KLOOP  /* while (KK != 16) */
KDONE:
/*
 *       Write out to C
 */

         #ifdef SREAL
            #if !defined(BETA0)
               lxv rA0, 16*0(pC)
               lxv rA1, 16*1(pC)
            #endif
            ATL_vvrsum4(rC00, rC10, rC20, rC30, rB0, rB1);
            #ifndef BETA0
               VADDC rC00, rC00, rA0
            #endif
            stxv rC00, 16*0(pC)
            ATL_vvrsum4(rC40, rC50, rC60, rC70, rB0, rB1);
            #ifndef BETA0
               VADDC rC40, rC40, rA1
            #endif
            stxv rC40, 16*1(pC)
            #if !defined(BETA0)
               lxv rA0, 16*2(pC)
               lxv rA1, 16*3(pC)
            #endif
            ATL_vvrsum4(rC01, rC11, rC21, rC31, rB0, rB1);
            #ifndef BETA0
               VADDC rC01, rC01, rA0
            #endif
            stxv rC01, 16*2(pC)
            ATL_vvrsum4(rC41, rC51, rC61, rC71, rB0, rB1);
            #ifndef BETA0
               VADDC rC41, rC41, rA1
            #endif
            stxv rC41, 16*3(pC)
            #if !defined(BETA0)
               lxv rA0, 16*4(pC)
               lxv rA1, 16*5(pC)
            #endif
            ATL_vvrsum4(rC02, rC12, rC22, rC32, rB0, rB1);
            #ifndef BETA0
               VADDC rC02, rC02, rA0
            #endif
            stxv rC02, 16*4(pC)
            ATL_vvrsum4(rC42, rC52, rC62, rC72, rB0, rB1);
            #ifndef BETA0
               VADDC rC42, rC42, rA1
            #endif
            stxv rC42, 16*5(pC)
            #if !defined(BETA0)
               lxv rA0, 16*6(pC)
               lxv rA1, 16*7(pC)
            #endif
            ATL_vvrsum4(rC03, rC13, rC23, rC33, rB0, rB1);
            #ifndef BETA0
               VADDC rC03, rC03, rA0
            #endif
            stxv rC03, 16*6(pC)
            ATL_vvrsum4(rC43, rC53, rC63, rC73, rB0, rB1);
            #ifndef BETA0
               VADDC rC43, rC43, rA1
            #endif
            stxv rC43, 16*7(pC)
            #if !defined(BETA0)
               lxv rA0, 16*8(pC)
               lxv rA1, 16*9(pC)
            #endif
            ATL_vvrsum4(rC04, rC14, rC24, rC34, rB0, rB1);
            #ifndef BETA0
               VADDC rC04, rC04, rA0
            #endif
            stxv rC04, 16*8(pC)
            ATL_vvrsum4(rC44, rC54, rC64, rC74, rB0, rB1);
            #ifndef BETA0
               VADDC rC44, rC44, rA1
            #endif
            stxv rC44, 16*9(pC)
            #if !defined(BETA0)
               lxv rA0, 16*10(pC)
               lxv rA1, 16*11(pC)
            #endif
            ATL_vvrsum4(rC05, rC15, rC25, rC35, rB0, rB1);
            #ifndef BETA0
               VADDC rC05, rC05, rA0
            #endif
            stxv rC05, 16*10(pC)
            ATL_vvrsum4(rC45, rC55, rC65, rC75, rB0, rB1);
            #ifndef BETA0
               VADDC rC45, rC45, rA1
            #endif
            stxv rC45, 16*11(pC)
         #else
            #if !defined(BETA0)
               lxv rA0, 16*0(pC)
               lxv rA1, 16*1(pC)
               lxv rA2, 16*2(pC)
               lxv rA3, 16*3(pC)
            #endif
            ATL_vvrsum2(rC00, rC10, rB0);
            #ifndef BETA0
               VADDC rC00, rC00, rA0
            #endif
            stxv rC00, 16*0(pC)
            ATL_vvrsum2(rC20, rC30, rB0);
            #ifndef BETA0
               VADDC rC20, rC20, rA1
            #endif
            stxv rC20, 16*1(pC)
            ATL_vvrsum2(rC40, rC50, rB0);
            #ifndef BETA0
               VADDC rC40, rC40, rA2
            #endif
            stxv rC40, 16*2(pC)
            ATL_vvrsum2(rC60, rC70, rB0);
            #ifndef BETA0
               VADDC rC60, rC60, rA3
            #endif
            stxv rC60, 16*3(pC)
            #if !defined(BETA0)
               lxv rA0, 16*4(pC)
               lxv rA1, 16*5(pC)
               lxv rA2, 16*6(pC)
               lxv rA3, 16*7(pC)
            #endif
            ATL_vvrsum2(rC01, rC11, rB0);
            #ifndef BETA0
               VADDC rC01, rC01, rA0
            #endif
            stxv rC01, 16*4(pC)
            ATL_vvrsum2(rC21, rC31, rB0);
            #ifndef BETA0
               VADDC rC21, rC21, rA1
            #endif
            stxv rC21, 16*5(pC)
            ATL_vvrsum2(rC41, rC51, rB0);
            #ifndef BETA0
               VADDC rC41, rC41, rA2
            #endif
            stxv rC41, 16*6(pC)
            ATL_vvrsum2(rC61, rC71, rB0);
            #ifndef BETA0
               VADDC rC61, rC61, rA3
            #endif
            stxv rC61, 16*7(pC)
            #if !defined(BETA0)
               lxv rA0, 16*8(pC)
               lxv rA1, 16*9(pC)
               lxv rA2, 16*10(pC)
               lxv rA3, 16*11(pC)
            #endif
            ATL_vvrsum2(rC02, rC12, rB0);
            #ifndef BETA0
               VADDC rC02, rC02, rA0
            #endif
            stxv rC02, 16*8(pC)
            ATL_vvrsum2(rC22, rC32, rB0);
            #ifndef BETA0
               VADDC rC22, rC22, rA1
            #endif
            stxv rC22, 16*9(pC)
            ATL_vvrsum2(rC42, rC52, rB0);
            #ifndef BETA0
               VADDC rC42, rC42, rA2
            #endif
            stxv rC42, 16*10(pC)
            ATL_vvrsum2(rC62, rC72, rB0);
            #ifndef BETA0
               VADDC rC62, rC62, rA3
            #endif
            stxv rC62, 16*11(pC)
            #if !defined(BETA0)
               lxv rA0, 16*12(pC)
               lxv rA1, 16*13(pC)
               lxv rA2, 16*14(pC)
               lxv rA3, 16*15(pC)
            #endif
            ATL_vvrsum2(rC03, rC13, rB0);
            #ifndef BETA0
               VADDC rC03, rC03, rA0
            #endif
            stxv rC03, 16*12(pC)
            ATL_vvrsum2(rC23, rC33, rB0);
            #ifndef BETA0
               VADDC rC23, rC23, rA1
            #endif
            stxv rC23, 16*13(pC)
            ATL_vvrsum2(rC43, rC53, rB0);
            #ifndef BETA0
               VADDC rC43, rC43, rA2
            #endif
            stxv rC43, 16*14(pC)
            ATL_vvrsum2(rC63, rC73, rB0);
            #ifndef BETA0
               VADDC rC63, rC63, rA3
            #endif
            stxv rC63, 16*15(pC)
            #if !defined(BETA0)
               lxv rA0, 16*16(pC)
               lxv rA1, 16*17(pC)
               lxv rA2, 16*18(pC)
               lxv rA3, 16*19(pC)
            #endif
            ATL_vvrsum2(rC04, rC14, rB0);
            #ifndef BETA0
               VADDC rC04, rC04, rA0
            #endif
            stxv rC04, 16*16(pC)
            ATL_vvrsum2(rC24, rC34, rB0);
            #ifndef BETA0
               VADDC rC24, rC24, rA1
            #endif
            stxv rC24, 16*17(pC)
            ATL_vvrsum2(rC44, rC54, rB0);
            #ifndef BETA0
               VADDC rC44, rC44, rA2
            #endif
            stxv rC44, 16*18(pC)
            ATL_vvrsum2(rC64, rC74, rB0);
            #ifndef BETA0
               VADDC rC64, rC64, rA3
            #endif
            stxv rC64, 16*19(pC)
            #if !defined(BETA0)
               lxv rA0, 16*20(pC)
               lxv rA1, 16*21(pC)
               lxv rA2, 16*22(pC)
               lxv rA3, 16*23(pC)
            #endif
            ATL_vvrsum2(rC05, rC15, rB0);
            #ifndef BETA0
               VADDC rC05, rC05, rA0
            #endif
            stxv rC05, 16*20(pC)
            ATL_vvrsum2(rC25, rC35, rB0);
            #ifndef BETA0
               VADDC rC25, rC25, rA1
            #endif
            stxv rC25, 16*21(pC)
            ATL_vvrsum2(rC45, rC55, rB0);
            #ifndef BETA0
               VADDC rC45, rC45, rA2
            #endif
            stxv rC45, 16*22(pC)
            ATL_vvrsum2(rC65, rC75, rB0);
            #ifndef BETA0
               VADDC rC65, rC65, rA3
            #endif
            stxv rC65, 16*23(pC)
         #endif

         mr pB, pBn
         add pBn, pBn, incBn
         mr pA, pA0
         addi pC, pC, 8*6*ELTSZ /* pC += mu*nu*sizeof */
      bdnz+ MLOOP               /* really N-loop */
      mtctr NNU
      mr pB, pB0
      ADDA pA0, pA0, incAm        /* pA += KB*mu*sizeof */
      mr pA, pA0
      add pBn, pB0, incBn
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
   ld r16, IOFF+8*2(r1)
   ld r17, IOFF+8*3(r1)
   ld r18, IOFF+8*4(r1)
/*
 * Dealloc stack frame & return
 */
   addi r1, r1, FSIZE
   blr

