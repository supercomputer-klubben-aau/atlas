/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_asm.h"

#if defined(SCPLX) && !defined(SREAL)
   #define SREAL 1
#endif

#ifdef SREAL
   #define xvsubdp xvsubsp
   #define xvadddp xvaddsp
   #define xvmaddadp xvmaddasp
   #define xvmuldp xvmulsp
   #define ELTSZ 4
   #define SSH   2             /* size shift, 2^2 = 4 */
   #define KMUL *2
   #define UK 32               /* unroll along K dim in elements */
   #define USH 5               /* unroll shift: 2^5 = 32 */
   #define VLEN 4              /* 4 floats in 1 vector */
   #define VSH  2              /* vector shift: 2^2 = 4 */
#else
   #define ELTSZ 8             /* sizeof(double) = 8 */
   #define SSH   3             /* size shift, 2^3 = 8 */
   #define UK 16               /* unroll along K dim in elements */
   #define USH 4               /* unroll shift: 2^4 = 16 */
   #define VLEN 2              /* 2 doubles in 1 vector */
   #define VSH  1              /* vector shift: 2^1 = 2 */
   #define KMUL
#endif
#define VECSZ 16
#ifdef BETA1
   #define ADD2C xvadddp
#elif !defined(BETA0)
   #define ADD2C xvsubdp
#endif
#define rC00 v0
#define rC10 v1
#define rC20 v2
#define rC30 v3
#define rC40 v4
#define rC50 v5
#define rC60 v6
#define rC01 v7
#define rC11 v8
#define rC21 v9
#define rC31 v10
#define rC41 v11
#define rC51 v12
#define rC61 v13
#define rC02 v14
#define rC12 v15
#define rC22 v16
#define rC32 v17
#define rC42 v18
#define rC52 v19
#define rC62 v20
#define rC03 v21
#define rC13 v22
#define rC23 v23
#define rC33 v24
#define rC43 v25
#define rC53 v26
#define rC63 v27
#define rC04 v28
#define rC14 v29
#define rC24 v30
#define rC34 v31
#define rC44 v32
#define rC54 v33
#define rC64 v34
#define rC05 v35
#define rC15 v36
#define rC25 v37
#define rC35 v38
#define rC45 v39
#define rC55 v40
#define rC65 v41
#define rC06 v42
#define rC16 v43
#define rC26 v44
#define rC36 v45
#define rC46 v46
#define rC56 v47
#define rC66 v48
#define rA0 v49
#define rB0 v56
#define rA1 v50
#define rB1 v57
#define rA2 v51
#define rB2 v58
#define rA3 v52
#define rB3 v59
#define rA4 v53
#define rB4 v60
#define rA5 v54
#define rB5 v61
#define rA6 v55
#define rB6 v62
/*
                     r3       r4     r5         r6         r7        r8
void ATL_USERMM(SZT nmu, SZT nnu, SZT K, CTYPE *pA, CTYPE *pB, TYPE *pC,
                       r9          r10     96(r1)
                CTYPE *pAn, CTYPE *pBn, CTYPE *pCn);
 */
#define pB0 r0
#define NMU r3
#define NNU0 r4
#define pA  r6
#define pB  r7
#define pC  r8
#define pfA r9
#define pfB r10
#define K16 r5  /* The number of 16*VLEN K-iterations to do */
#define KR r11  /* The number of VLEN K-its (cleanup) to do */
#define pA0 r12
#define pBn r30
#define incAm r29
#define incBn r28
#define NNU r27
#define i128_0 0
#define i128_1 r14
#define i128_2 r15
/* #define Kb0 r14 */
#if defined(__BIG_ENDIAN__)
   #define BIGENDIEN 1
#elif defined(__ORDER_BIG_ENDIAN__) && defined(__BYTE_ORDER__)
   #if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
      #define BIGENDIEN 1
   #endif
#endif
#ifdef SREAL
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
      #define ATL_vvrsum1(s0_, t_) \
         xxpermdi t_, s0_, s0_, 0 ; \
         xvaddsp s0_, s0_, t_ ; \
         xxspltw  t_, s0_, 1 ; \
         xvaddsp s0_, s0_, t_ ;
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
      #define ATL_vvrsum1(s0_, t_) \
         xxpermdi t_, s0_, s0_, 0 ; \
         xvaddsp s0_, s0_, t_ ; \
         xxspltw  t_, s0_, 2 ; \
         xvaddsp s0_, s0_, t_ ;
    #endif
#else
   #ifdef BIGENDIEN
      #define ATL_vvrsum2(s0_, s1_, t_) \
         xxpermdi t_, s0_, s1_, 0 ; \
         xxpermdi s0_, s0_, s1_, 3 ; \
         xvadddp s0_, s0_, t_
      #define ATL_vvrsum1(s0_, t_) \
         xxpermdi t_, s0_, s0_, 3 ; \
         xvadddp s0_, s0_, t_
   #else
      #define ATL_vvrsum2(s1_, s0_, t_) \
         xxpermdi t_, s0_, s1_, 0 ; \
         xxpermdi s0_, s0_, s1_, 3 ; \
         xvadddp s1_, s0_, t_
      #define ATL_vvrsum1(s0_, t_) \
         xxpermdi t_, s0_, s0_, 0 ; \
         xvadddp s0_, s0_, t_
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
#define ISZ 8*8     /* only need 7, but need 8 for 16-byte stack alignment */
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
   std i128_1, IOFF+8*0(r1)
   li  i128_1, 128*1
   std i128_2, IOFF+8*1(r1)
   li  i128_2, 128*2
   std NNU, IOFF+8*3(r1)
   std incBn, IOFF+8*4(r1)
   std incAm, IOFF+8*5(r1)
   std pBn, IOFF+8*6(r1)
/*
 * Save vec regs using AltiVec instructions
 */
   addi pA0, r1, VOFF   /* pA0 points to start of vector save area */
   stvx 20, r0, pA0
      addi pA0, pA0, 16
   stvx 21, r0, pA0
      addi pA0, pA0, 16
   stvx 22, r0, pA0
      addi pA0, pA0, 16
   stvx 23, r0, pA0
      addi pA0, pA0, 16
   stvx 24, r0, pA0
      addi pA0, pA0, 16
   stvx 25, r0, pA0
      addi pA0, pA0, 16
   stvx 26, r0, pA0
      addi pA0, pA0, 16
   stvx 27, r0, pA0
      addi pA0, pA0, 16
   stvx 28, r0, pA0
      addi pA0, pA0, 16
   stvx 29, r0, pA0
      addi pA0, pA0, 16
   stvx 30, r0, pA0
      addi pA0, pA0, 16
   stvx 31, r0, pA0
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
 * incAm = mu*K*sizeof = 7*(K<<SSH); tmp = K<<SSH; incAm=8*tmp-tmp;
 */
   sldi pA0, K16, SSH     /* pA0 = K*sizeof(TYPE) */
   sldi incAm, K16, 3+SSH /* incAm = 8*K*sizeof(TYPE) */
   sub incAm, incAm, pA0  /* incAm = 7*K*sizeof(TYPE) = mu*K*sizeof(TYPE) */
   addi K16, K16, -VLEN   /* peeled VLEN iterations */
   srdi K16, K16, VSH     /* K16 now has number of vector K after peel */
   andi. KR, K16, 0xF     /* KR = K%16 */
   srdi K16, K16, 4       /* K16 = K/16 */
   cmpdi cr6, KR, 0       /* check if all we do is K=2 peeled iterations */
   cmpdi cr5, K16, 0      /* check if we do any 16-vec unroll its */
   mr pA0, pA
   mr pB0, pB

   MLOOP:
      mr NNU, NNU0
      lxv rB0, 0(pB)
      lxv rB1, 1*VECSZ(pB)
      lxv rB2, 2*VECSZ(pB)
      lxv rB3, 3*VECSZ(pB)
      lxv rB4, 4*VECSZ(pB)
      lxv rB5, 5*VECSZ(pB)
      NLOOP:
         lxv rA0, 0(pA)
         lxv rA1, 1*VECSZ(pA)
         lxv rA2, 2*VECSZ(pA)
         lxv rA3, 3*VECSZ(pA)

         xvmuldp rC00, rA0, rB0
            lxv rA4, 4*VECSZ(pA)
         xvmuldp rC10, rA1, rB0
            lxv rA5, 5*VECSZ(pA)
         xvmuldp rC20, rA2, rB0
            lxv rA6, 6*VECSZ(pA)
         xvmuldp rC30, rA3, rB0
            lxv rB6, 6*VECSZ(pB)
         xvmuldp rC40, rA4, rB0
         xvmuldp rC50, rA5, rB0
         xvmuldp rC60, rA6, rB0
            lxv rB0, 7*VECSZ(pB)
         xvmuldp rC01, rA0, rB1
         xvmuldp rC11, rA1, rB1
         xvmuldp rC21, rA2, rB1
         xvmuldp rC31, rA3, rB1
         xvmuldp rC41, rA4, rB1
         xvmuldp rC51, rA5, rB1
         xvmuldp rC61, rA6, rB1
            lxv rB1, 8*VECSZ(pB)
         xvmuldp rC02, rA0, rB2
         xvmuldp rC12, rA1, rB2
         xvmuldp rC22, rA2, rB2
         xvmuldp rC32, rA3, rB2
         xvmuldp rC42, rA4, rB2
         xvmuldp rC52, rA5, rB2
         xvmuldp rC62, rA6, rB2
            lxv rB2, 9*VECSZ(pB)
         xvmuldp rC03, rA0, rB3
         xvmuldp rC13, rA1, rB3
         xvmuldp rC23, rA2, rB3
         xvmuldp rC33, rA3, rB3
         xvmuldp rC43, rA4, rB3
         xvmuldp rC53, rA5, rB3
         xvmuldp rC63, rA6, rB3
            lxv rB3, 10*VECSZ(pB)
         xvmuldp rC04, rA0, rB4
         xvmuldp rC14, rA1, rB4
         xvmuldp rC24, rA2, rB4
         xvmuldp rC34, rA3, rB4
         xvmuldp rC44, rA4, rB4
         xvmuldp rC54, rA5, rB4
         xvmuldp rC64, rA6, rB4
            lxv rB4, 11*VECSZ(pB)
         xvmuldp rC05, rA0, rB5
         xvmuldp rC15, rA1, rB5
         xvmuldp rC25, rA2, rB5
         xvmuldp rC35, rA3, rB5
         xvmuldp rC45, rA4, rB5
         xvmuldp rC55, rA5, rB5
         xvmuldp rC65, rA6, rB5
            lxv rB5, 12*VECSZ(pB)
         xvmuldp rC06, rA0, rB6
            lxv rA0, 7*VECSZ(pA)
         xvmuldp rC16, rA1, rB6
            lxv rA1, 8*VECSZ(pA)
         xvmuldp rC26, rA2, rB6
            lxv rA2, 9*VECSZ(pA)
         xvmuldp rC36, rA3, rB6
            lxv rA3, 10*VECSZ(pA)
         xvmuldp rC46, rA4, rB6
            lxv rA4, 11*VECSZ(pA)
         xvmuldp rC56, rA5, rB6
            lxv rA5, 12*VECSZ(pA)
         xvmuldp rC66, rA6, rB6
            lxv rA6, 13*VECSZ(pA)
         addi pA, pA, 7*VECSZ       /* pA += mu*sizeof */
         addi pB, pB, 7*VECSZ       /* pB += mu*sizeof */
         beq cr5, KCLEAN            /* if (K16==0) skip to cleanup loop */
         mtctr K16                  /* move # of 16*vlen K its to ctr */
         KLOOP16:
            /* Vector unroll=0 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 6*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 7*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 8*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 9*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 10*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 11*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 12*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 7*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 8*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 9*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 10*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 11*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 12*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 13*VECSZ(pA)

            /* Vector unroll=1 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 13*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 14*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 15*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 16*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 17*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 18*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 19*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 14*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 15*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 16*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 17*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 18*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 19*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 20*VECSZ(pA)

            /* Vector unroll=2 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 20*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 21*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 22*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 23*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 24*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 25*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 26*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 21*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 22*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 23*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 24*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 25*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 26*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 27*VECSZ(pA)

            /* Vector unroll=3 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 27*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 28*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 29*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 30*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 31*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 32*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 33*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 28*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 29*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 30*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 31*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 32*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 33*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 34*VECSZ(pA)

            /* Vector unroll=4 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 34*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 35*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 36*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 37*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 38*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 39*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 40*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 35*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 36*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 37*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 38*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 39*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 40*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 41*VECSZ(pA)

            /* Vector unroll=5 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 41*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 42*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 43*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 44*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 45*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 46*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 47*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 42*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 43*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 44*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 45*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 46*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 47*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 48*VECSZ(pA)

            /* Vector unroll=6 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 48*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 49*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 50*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 51*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 52*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 53*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 54*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 49*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 50*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 51*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 52*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 53*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 54*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 55*VECSZ(pA)

            /* Vector unroll=7 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 55*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 56*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 57*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 58*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 59*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 60*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 61*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 56*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 57*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 58*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 59*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 60*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 61*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 62*VECSZ(pA)

            /* Vector unroll=8 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 62*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 63*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 64*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 65*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 66*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 67*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 68*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 63*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 64*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 65*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 66*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 67*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 68*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 69*VECSZ(pA)

            /* Vector unroll=9 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 69*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 70*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 71*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 72*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 73*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 74*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 75*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 70*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 71*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 72*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 73*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 74*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 75*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 76*VECSZ(pA)

            /* Vector unroll=10 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 76*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 77*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 78*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 79*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 80*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 81*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 82*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 77*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 78*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 79*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 80*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 81*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 82*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 83*VECSZ(pA)

            /* Vector unroll=11 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 83*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 84*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 85*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 86*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 87*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 88*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 89*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 84*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 85*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 86*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 87*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 88*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 89*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 90*VECSZ(pA)

            /* Vector unroll=12 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 90*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 91*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 92*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 93*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 94*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 95*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 96*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 91*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 92*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 93*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 94*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 95*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 96*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 97*VECSZ(pA)

            /* Vector unroll=13 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 97*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 98*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 99*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 100*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 101*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 102*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 103*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 98*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 99*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 100*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 101*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 102*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 103*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 104*VECSZ(pA)

            /* Vector unroll=14 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 104*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 105*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 106*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 107*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 108*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 109*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 110*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 105*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 106*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 107*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 108*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 109*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 110*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 111*VECSZ(pA)

            /* Vector unroll=15 */
            xvmaddadp rC00, rA0, rB0
                     lxv rB6, 111*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
                     lxv rB0, 112*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
                     lxv rB1, 113*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
                     lxv rB2, 114*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
                     lxv rB3, 115*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
                     lxv rB4, 116*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
                     lxv rB5, 117*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
                     lxv rA0, 112*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
                     lxv rA1, 113*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
                     lxv rA2, 114*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
                     lxv rA3, 115*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
                     lxv rA4, 116*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
                     lxv rA5, 117*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
                     lxv rA6, 118*VECSZ(pA)

            addi pA, pA, 16*7*VECSZ       /* pA += 16*mu*sizeof */
            addi pB, pB, 16*7*VECSZ       /* pB += 16*nu*sizeof */
         bdnz KLOOP16

         KCLEAN:
         beq cr6, KDONE             /* if (KR==0) done with K loop */
         mtctr KR  /* move rolled VLEN K its to do to KR */
         KLOOP1:
            xvmaddadp rC00, rA0, rB0
               lxv rB6, 6*VECSZ(pB)
            xvmaddadp rC10, rA1, rB0
            xvmaddadp rC20, rA2, rB0
            xvmaddadp rC30, rA3, rB0
            xvmaddadp rC40, rA4, rB0
            xvmaddadp rC50, rA5, rB0
            xvmaddadp rC60, rA6, rB0
               lxv rB0, 7*VECSZ(pB)
            xvmaddadp rC01, rA0, rB1
            xvmaddadp rC11, rA1, rB1
            xvmaddadp rC21, rA2, rB1
            xvmaddadp rC31, rA3, rB1
            xvmaddadp rC41, rA4, rB1
            xvmaddadp rC51, rA5, rB1
            xvmaddadp rC61, rA6, rB1
               lxv rB1, 8*VECSZ(pB)
            xvmaddadp rC02, rA0, rB2
            xvmaddadp rC12, rA1, rB2
            xvmaddadp rC22, rA2, rB2
            xvmaddadp rC32, rA3, rB2
            xvmaddadp rC42, rA4, rB2
            xvmaddadp rC52, rA5, rB2
            xvmaddadp rC62, rA6, rB2
               lxv rB2, 9*VECSZ(pB)
            xvmaddadp rC03, rA0, rB3
            xvmaddadp rC13, rA1, rB3
            xvmaddadp rC23, rA2, rB3
            xvmaddadp rC33, rA3, rB3
            xvmaddadp rC43, rA4, rB3
            xvmaddadp rC53, rA5, rB3
            xvmaddadp rC63, rA6, rB3
               lxv rB3, 10*VECSZ(pB)
            xvmaddadp rC04, rA0, rB4
            xvmaddadp rC14, rA1, rB4
            xvmaddadp rC24, rA2, rB4
            xvmaddadp rC34, rA3, rB4
            xvmaddadp rC44, rA4, rB4
            xvmaddadp rC54, rA5, rB4
            xvmaddadp rC64, rA6, rB4
               lxv rB4, 11*VECSZ(pB)
            xvmaddadp rC05, rA0, rB5
            xvmaddadp rC15, rA1, rB5
            xvmaddadp rC25, rA2, rB5
            xvmaddadp rC35, rA3, rB5
            xvmaddadp rC45, rA4, rB5
            xvmaddadp rC55, rA5, rB5
            xvmaddadp rC65, rA6, rB5
               lxv rB5, 12*VECSZ(pB)
            xvmaddadp rC06, rA0, rB6
               lxv rA0, 7*VECSZ(pA)
            xvmaddadp rC16, rA1, rB6
               lxv rA1, 8*VECSZ(pA)
            xvmaddadp rC26, rA2, rB6
               lxv rA2, 9*VECSZ(pA)
            xvmaddadp rC36, rA3, rB6
               lxv rA3, 10*VECSZ(pA)
            xvmaddadp rC46, rA4, rB6
               lxv rA4, 11*VECSZ(pA)
            xvmaddadp rC56, rA5, rB6
               lxv rA5, 12*VECSZ(pA)
            xvmaddadp rC66, rA6, rB6
               lxv rA6, 13*VECSZ(pA)
            addi pA, pA, 7*VECSZ       /* pA += mu*sizeof */
            addi pB, pB, 7*VECSZ       /* pB += mu*sizeof */
         bdnz+ KLOOP1
         KDONE:
	 #ifdef SREAL
	    #ifndef BETA0
	       lxv rA0, 0*VECSZ(pC)
	       lxv rA1, 1*VECSZ(pC)
	       lxv rA2, 2*VECSZ(pC)
	       lxv rA3, 3*VECSZ(pC)
	    #endif
	    ATL_vvrsum4(rC00, rC10, rC20, rC30, rA4, rA5)
	    ATL_vvrsum4(rC40, rC50, rC60, rC01, rA4, rA5)
	    #ifndef BETA0
               ADD2C rC00, rC00, rA0
               ADD2C rC40, rC40, rA1
            #endif
	    ATL_vvrsum4(rC11, rC21, rC31, rC41, rA4, rA5)
	    #ifndef BETA0
               ADD2C rC11, rC11, rA2
            #endif
	    stxv rC00, 0*VECSZ(pC)
	    stxv rC40, 1*VECSZ(pC)
	    ATL_vvrsum4(rC51, rC61, rC02, rC12, rA4, rA5)
	    #ifndef BETA0
               ADD2C rC51, rC51, rA3
	    #endif
	    #ifndef BETA0
	       lxv rA0, 4*VECSZ(pC)
	       lxv rA1, 5*VECSZ(pC)
	       lxv rA2, 6*VECSZ(pC)
	       lxv rA3, 7*VECSZ(pC)
	    #endif
	    ATL_vvrsum4(rC22, rC32, rC42, rC52, rA4, rA5)
	    stxv rC11, 2*VECSZ(pC)
	    stxv rC51, 3*VECSZ(pC)
	    ATL_vvrsum4(rC62, rC03, rC13, rC23, rA4, rA5)
            #ifndef BETA0
               ADD2C rC22, rC22, rA0
            #endif
	    ATL_vvrsum4(rC33, rC43, rC53, rC63, rA4, rA5)
            #ifndef BETA0
               ADD2C rC62, rC62, rA1
            #endif
	    stxv rC22, 4*VECSZ(pC)
	    ATL_vvrsum4(rC04, rC14, rC24, rC34, rA4, rA5)
	    stxv rC62, 5*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 8*VECSZ(pC)
               ADD2C rC33, rC33, rA2
	       lxv rA1, 9*VECSZ(pC)
               ADD2C rC04, rC04, rA3
	    #endif
	    stxv rC33, 6*VECSZ(pC)
	    stxv rC04, 7*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA2, 10*VECSZ(pC)
	       lxv rA3, 11*VECSZ(pC)
	    #endif
	    ATL_vvrsum4(rC44, rC54, rC64, rC05, rA4, rA5)
	    ATL_vvrsum4(rC15, rC25, rC35, rC45, rA4, rA5)
	    ATL_vvrsum4(rC55, rC65, rC06, rC16, rA4, rA5)
	    ATL_vvrsum4(rC26, rC36, rC46, rC56, rA4, rA5)
	    #ifndef BETA0
               ADD2C rC44, rC44, rA0
               ADD2C rC15, rC15, rA1
               ADD2C rC55, rC55, rA2
               ADD2C rC26, rC26, rA3
	    #endif
	    stxv rC44, 8*VECSZ(pC)
	    stxv rC15, 9*VECSZ(pC)
	    stxv rC55, 10*VECSZ(pC)
	    stxv rC26, 11*VECSZ(pC)
	 #else   /* double precision code */
	    #ifndef BETA0
	       lxv rA0, 0*VECSZ(pC)
	       lxv rA1, 1*VECSZ(pC)
	       lxv rA2, 2*VECSZ(pC)
	       lxv rA3, 3*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC00, rC10, rA4)
	    ATL_vvrsum2(rC20, rC30, rA6)
	    ATL_vvrsum2(rC40, rC50, rA5)
	    ATL_vvrsum2(rC60, rC01, rA4)
	    #ifndef BETA0
               ADD2C rC00, rC00, rA0
               ADD2C rC20, rC20, rA1
               ADD2C rC40, rC40, rA2
               ADD2C rC60, rC60, rA3
	    #endif
	    stxv rC00, 0*VECSZ(pC)
	    stxv rC20, 1*VECSZ(pC)
	    stxv rC40, 2*VECSZ(pC)
	    stxv rC60, 3*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 4*VECSZ(pC)
	       lxv rA1, 5*VECSZ(pC)
	       lxv rA2, 6*VECSZ(pC)
	       lxv rA3, 7*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC11, rC21, rA4)
	    ATL_vvrsum2(rC31, rC41, rA6)
	    ATL_vvrsum2(rC51, rC61, rA5)
	    ATL_vvrsum2(rC02, rC12, rA4)
	    #ifndef BETA0
               ADD2C rC11, rC11, rA0
               ADD2C rC31, rC31, rA1
               ADD2C rC51, rC51, rA2
               ADD2C rC02, rC02, rA3
	    #endif
	    stxv rC11, 4*VECSZ(pC)
	    stxv rC31, 5*VECSZ(pC)
	    stxv rC51, 6*VECSZ(pC)
	    stxv rC02, 7*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 8*VECSZ(pC)
	       lxv rA1, 9*VECSZ(pC)
	       lxv rA2, 10*VECSZ(pC)
	       lxv rA3, 11*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC22, rC32, rA4)
	    ATL_vvrsum2(rC42, rC52, rA6)
	    ATL_vvrsum2(rC62, rC03, rA5)
	    ATL_vvrsum2(rC13, rC23, rA4)
	    #ifndef BETA0
               ADD2C rC22, rC22, rA0
               ADD2C rC42, rC42, rA1
               ADD2C rC62, rC62, rA2
               ADD2C rC13, rC13, rA3
	    #endif
	    stxv rC22, 8*VECSZ(pC)
	    stxv rC42, 9*VECSZ(pC)
	    stxv rC62, 10*VECSZ(pC)
	    stxv rC13, 11*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 12*VECSZ(pC)
	       lxv rA1, 13*VECSZ(pC)
	       lxv rA2, 14*VECSZ(pC)
	       lxv rA3, 15*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC33, rC43, rA4)
	    ATL_vvrsum2(rC53, rC63, rA6)
	    ATL_vvrsum2(rC04, rC14, rA5)
	    ATL_vvrsum2(rC24, rC34, rA4)
	    #ifndef BETA0
               ADD2C rC33, rC33, rA0
               ADD2C rC53, rC53, rA1
               ADD2C rC04, rC04, rA2
               ADD2C rC24, rC24, rA3
	    #endif
	    stxv rC33, 12*VECSZ(pC)
	    stxv rC53, 13*VECSZ(pC)
	    stxv rC04, 14*VECSZ(pC)
	    stxv rC24, 15*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 16*VECSZ(pC)
	       lxv rA1, 17*VECSZ(pC)
	       lxv rA2, 18*VECSZ(pC)
	       lxv rA3, 19*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC44, rC54, rA4)
	    ATL_vvrsum2(rC64, rC05, rA6)
	    ATL_vvrsum2(rC15, rC25, rA5)
	    ATL_vvrsum2(rC35, rC45, rA4)
	    #ifndef BETA0
               ADD2C rC44, rC44, rA0
               ADD2C rC64, rC64, rA1
               ADD2C rC15, rC15, rA2
               ADD2C rC35, rC35, rA3
	    #endif
	    stxv rC44, 16*VECSZ(pC)
	    stxv rC64, 17*VECSZ(pC)
	    stxv rC15, 18*VECSZ(pC)
	    stxv rC35, 19*VECSZ(pC)
	    #ifndef BETA0
	       lxv rA0, 20*VECSZ(pC)
	       lxv rA1, 21*VECSZ(pC)
	       lxv rA2, 22*VECSZ(pC)
	       lxv rA3, 23*VECSZ(pC)
	    #endif
	    ATL_vvrsum2(rC55, rC65, rA4)
	    ATL_vvrsum2(rC06, rC16, rA6)
	    ATL_vvrsum2(rC26, rC36, rA5)
	    ATL_vvrsum2(rC46, rC56, rA4)
	    #ifndef BETA0
               ADD2C rC55, rC55, rA0
               ADD2C rC06, rC06, rA1
               ADD2C rC26, rC26, rA2
               ADD2C rC46, rC46, rA3
	    #endif
	    stxv rC55, 20*VECSZ(pC)
	    stxv rC06, 21*VECSZ(pC)
	    stxv rC26, 22*VECSZ(pC)
	    stxv rC46, 23*VECSZ(pC)
         #endif
/*
 *       Handle 49th register; could do with scalar ops, but vec allowed
 */
         #ifndef BETA0
            lxv rA0, 48*ELTSZ(pC)
         #endif
         ATL_vvrsum1(rC66, rA4)
         #ifndef BETA0
             ADD2C rC66, rC66, rA0
         #endif
         stxv rC66, 48*ELTSZ(pC)
/*
 *       Go to next loop iteration
 */
         #ifdef SREAL
            addi  pC, pC, 52*ELTSZ
         #else
            addi  pC, pC, 50*ELTSZ
         #endif
         mr pA, pA0
         addic. NNU, NNU, -1
      bne+ NLOOP
      mr pB, pB0
      add pA0, pA0, incAm
      mr pA, pA0
      addic. NMU, NMU, -1
   bne+ MLOOP
/*
 * Restore vec regs using Altivec instructions
 */
   addi pA0, r1, VOFF   /* pA0 points to start of vector save area */
   lvx 20, r0, pA0
      addi pA0, pA0, 16
   lvx 21, r0, pA0
      addi pA0, pA0, 16
   lvx 22, r0, pA0
      addi pA0, pA0, 16
   lvx 23, r0, pA0
      addi pA0, pA0, 16
   lvx 24, r0, pA0
      addi pA0, pA0, 16
   lvx 25, r0, pA0
      addi pA0, pA0, 16
   lvx 26, r0, pA0
      addi pA0, pA0, 16
   lvx 27, r0, pA0
      addi pA0, pA0, 16
   lvx 28, r0, pA0
      addi pA0, pA0, 16
   lvx 29, r0, pA0
      addi pA0, pA0, 16
   lvx 30, r0, pA0
      addi pA0, pA0, 16
   lvx 31, r0, pA0
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
   ld i128_1, IOFF+8*0(r1)
   ld i128_2, IOFF+8*1(r1)
   ld NNU, IOFF+8*3(r1)
   ld incBn, IOFF+8*4(r1)
   ld incAm, IOFF+8*5(r1)
   ld pBn, IOFF+8*6(r1)
/*
 * Dealloc stack frame & return
 */
   addi r1, r1, FSIZE
   blr
