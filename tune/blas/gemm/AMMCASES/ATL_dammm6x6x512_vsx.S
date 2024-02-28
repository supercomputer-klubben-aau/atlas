/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_asm.h"
#define pB0 r0
#define NNU r2
#define NMU r3
#define NN0 r4
#define KK  r5
#define pA  r6
#define pB  r7
#define pC  r8
#define pfA r9
#define pfB r10
#define incA r11
/* #define Kb0 r14 */
#define i16 r12
   #define i16_1 i16
#define i32 r15
   #define i16_2 i32
#define i48 r16
   #define i16_3 i48
#define i64 r17
   #define i16_4 i64
#define i80 r18
   #define i16_5 i80
#define i96 r19
   #define i16_6 i96
#define i112 r20
   #define i16_7 i112
#define i128 r22
   #define i16_8 i128
#define i144 r23
   #define i16_9 i144
#define i160 r24
   #define i16_10 i160
#define i176 r25
   #define i16_11 i176
#define i192 r26
   #define i16_12 i192
#define i208 r27
   #define i16_13 i208
#define i224 r28
   #define i16_14 i224
#define i240 r29
   #define i16_15 i240

#define rA0 v0
#define rA1 v1
#define rA2 v2
#define rA3 v3
#define rA4 v4
#define rA5 v5
#define ra0 v6
#define ra1 v7
#define ra2 v8
#define ra3 v9
#define ra4 v10
#define ra5 v11
#define rB0 v12
#define rB1 v13
#define rB2 v14
#define rB3 v15
#define rB4 v16
#define rB5 v17
#define rb0 v18
#define rb1 v19
#define rb2 v20
#define rb3 v21
#define rb4 v22
#define rb5 v23
#define rC00 v24
#define rC10 v25
#define rC20 v26
#define rC30 v27
#define rC40 v28
#define rC50 v29
#define rC01 v30
#define rC11 v31
#define rC21 v32
#define rC31 v33
#define rC41 v34
#define rC51 v35
#define rC02 v36
#define rC12 v37
#define rC22 v38
#define rC32 v39
#define rC42 v40
#define rC52 v41
#define rC03 v42
#define rC13 v43
#define rC23 v44
#define rC33 v45
#define rC43 v46
#define rC53 v47
#define rC04 v48
#define rC14 v49
#define rC24 v50
#define rC34 v51
#define rC44 v52
#define rC54 v53
#define rC05 v54
#define rC15 v55
#define rC25 v56
#define rC35 v57
#define rC45 v58
#define rC55 v59
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
.local MLOOP
.local NLOOP
.local KLOOP
.globl ATL_USERMM
ATL_USERMM:
/*
 * Must save & restore:
 *    r1&r2 -> in sframe tail
 *    r14-r31 -> 18*8 = 144 bytes (if used)
 *    f14-f31 -> vsr14-vsr31 -> 18*8 = 144 bytes
 *    v20-v31 -> vsr52-vsr63 -> 12*16 = 192 bytes
 */
#define ISZ 0
#define FSZ 144
#define VSZ 192
#define LSZ 0
#define MSZ 32      /* mandatory 32 bytes for stack frame tail */
#define VOFF (MSZ+LSZ)
#define FOFF (VOFF+VSZ)
#define IOFF (FOFF+FSZ)
#define FSIZE (MSZ+LSZ+VSZ+FSZ+ISZ)
/*
 * Allocate local stack frame in order to save & restore
 */
   stdu r1, -FSIZE(r1)   /* allocate stack & store old stack ptr */
/*
 * Save iregs, UNUSED: r14
 */
   std r2, 24(r1)        /* TOC reg always at 24(sp) */
/*
 * Save f14-f31
 */
   stfd f14, FOFF(r1)
   stfd f15, 8+FOFF(r1)
   stfd f16, 16+FOFF(r1)
   stfd f17, 24+FOFF(r1)
   stfd f18, 32+FOFF(r1)
   stfd f19, 40+FOFF(r1)
   stfd f20, 48+FOFF(r1)
   stfd f21, 56+FOFF(r1)
   stfd f22, 64+FOFF(r1)
   stfd f23, 72+FOFF(r1)
   stfd f24, 80+FOFF(r1)
   stfd f25, 88+FOFF(r1)
   stfd f26, 96+FOFF(r1)
   stfd f27, 104+FOFF(r1)
   stfd f28, 112+FOFF(r1)
   stfd f29, 120+FOFF(r1)
   stfd f30, 128+FOFF(r1)
   stfd f31, 136+FOFF(r1)
/*
 * Save vec regs, use VSX so we don't rely on Altivec
 */
   mr pB0, pB
   addi pB, r1, VOFF
   stxv v52, 0(pB)
   stxv v53, 16*1(pB)
   stxv v54, 16*2(pB)
   stxv v55, 16*3(pB)
   stxv v56, 16*4(pB)
   stxv v57, 16*5(pB)
   stxv v58, 16*6(pB)
   stxv v59, 16*7(pB)
   #if 0
   stxv v60, 16*8(pB)
   stxv v61, 16*9(pB)
   stxv v62, 16*10(pB)
   stxv v63, 16*11(pB)
   #endif

   ALIGN16
   MLOOP:
      mr pB, pB0
      mr NNU, NN0
      NLOOP:
/*
 *       Need to peel to handle beta=0,-1, but I'm too lazy to do it yet
 */
         #ifdef BETA0
            xxlxor rC00, rC00, rC00
            xxlxor rC10, rC10, rC10
            xxlxor rC20, rC20, rC20
            xxlxor rC30, rC30, rC30
            xxlxor rC40, rC40, rC40
            xxlxor rC50, rC50, rC50
            xxlxor rC01, rC01, rC01
            xxlxor rC11, rC11, rC11
            xxlxor rC21, rC21, rC21
            xxlxor rC31, rC31, rC31
            xxlxor rC41, rC41, rC41
            xxlxor rC51, rC51, rC51
            xxlxor rC02, rC02, rC02
            xxlxor rC12, rC12, rC12
            xxlxor rC22, rC22, rC22
            xxlxor rC32, rC32, rC32
            xxlxor rC42, rC42, rC42
            xxlxor rC52, rC52, rC52
            xxlxor rC03, rC03, rC03
            xxlxor rC13, rC13, rC13
            xxlxor rC23, rC23, rC23
            xxlxor rC33, rC33, rC33
            xxlxor rC43, rC43, rC43
            xxlxor rC53, rC53, rC53
            xxlxor rC04, rC04, rC04
            xxlxor rC14, rC14, rC14
            xxlxor rC24, rC24, rC24
            xxlxor rC34, rC34, rC34
            xxlxor rC44, rC44, rC44
            xxlxor rC54, rC54, rC54
            xxlxor rC05, rC05, rC05
            xxlxor rC15, rC15, rC15
            xxlxor rC25, rC25, rC25
            xxlxor rC35, rC35, rC35
            xxlxor rC45, rC45, rC45
            xxlxor rC55, rC55, rC55
         #else
/*            lxvd2x rC00, 0, pC */
            #if 0
            lxsd rC00, 0(pC)
            lxsd rC10, 8(pC)
            lxsd rC20, 8*2(pC)
            lxsd rC30, 8*3(pC)
            lxsd rC40, 8*4(pC)
            lxsd rC50, 8*5(pC)
            lxsd rC01, 8*6(pC)
            lxsd rC11, 8*7(pC)
            lxsd rC21, 8*8(pC)
            lxsd rC31, 8*9(pC)
            lxsd rC41, 8*10(pC)
            lxsd rC51, 8*11(pC)
            lxsd rC02, 8*12(pC)
            lxsd rC12, 8*13(pC)
            lxsd rC22, 8*14(pC)
            lxsd rC32, 8*15(pC)
            lxsd rC42, 8*8(pC)
            lxsd rC52, 8*17(pC)
            lxsd rC03, 8*18(pC)
            lxsd rC13, 8*19(pC)
            lxsd rC23, 8*20(pC)
            lxsd rC33, 8*21(pC)
            lxsd rC43, 8*22(pC)
            lxsd rC53, 8*23(pC)
            lxsd rC04, 8*24(pC)
            lxsd rC14, 8*25(pC)
            lxsd rC24, 8*26(pC)
            lxsd rC34, 8*27(pC)
            lxsd rC44, 8*28(pC)
            lxsd rC54, 8*29(pC)
            lxsd rC05, 8*30(pC)
            lxsd rC15, 8*31(pC)
            lxsd rC25, 8*32(pC)
            lxsd rC35, 8*33(pC)
            lxsd rC45, 8*34(pC)
            lxsd rC55, 8*35(pC)
            #else
            lxv rC10, 0(pC)
            xxspltd rC00, rC10, 1
            lxv rC30, 16(pC)
            xxspltd rC20, rC30, 1
            lxv rC50, 32(pC)
            xxspltd rC40, rC50, 1
            lxv rC11, 48(pC)
            xxspltd rC01, rC11, 1
            lxv rC31, 64(pC)
            xxspltd rC21, rC31, 1
            lxv rC51, 80(pC)
            xxspltd rC41, rC51, 1
            lxv rC12, 96(pC)
            xxspltd rC02, rC12, 1
            lxv rC32, 112(pC)
            xxspltd rC22, rC32, 1
            lxv rC52, 128(pC)
            xxspltd rC42, rC52, 1
            lxv rC13, 144(pC)
            xxspltd rC03, rC13, 1
            lxv rC33, 160(pC)
            xxspltd rC23, rC33, 1
            lxv rC53, 176(pC)
            xxspltd rC43, rC53, 1
            lxv rC14, 192(pC)
            xxspltd rC04, rC14, 1
            lxv rC34, 208(pC)
            xxspltd rC24, rC34, 1
            lxv rC54, 224(pC)
            xxspltd rC44, rC54, 1
            lxv rC15, 240(pC)
            xxspltd rC05, rC15, 1
            lxv rC35, 256(pC)
            xxspltd rC25, rC35, 1
            lxv rC55, 272(pC)
            xxspltd rC45, rC55, 1
            #endif

            #ifndef BETA1
               xsnegdp rC00, rC00
               xsnegdp rC10, rC10
               xsnegdp rC20, rC20
               xsnegdp rC30, rC30
               xsnegdp rC40, rC40
               xsnegdp rC50, rC50
               xsnegdp rC01, rC01
               xsnegdp rC11, rC11
               xsnegdp rC21, rC21
               xsnegdp rC31, rC31
               xsnegdp rC41, rC41
               xsnegdp rC51, rC51
               xsnegdp rC02, rC02
               xsnegdp rC12, rC12
               xsnegdp rC22, rC22
               xsnegdp rC32, rC32
               xsnegdp rC42, rC42
               xsnegdp rC52, rC52
               xsnegdp rC03, rC03
               xsnegdp rC13, rC13
               xsnegdp rC23, rC23
               xsnegdp rC33, rC33
               xsnegdp rC43, rC43
               xsnegdp rC53, rC53
               xsnegdp rC04, rC04
               xsnegdp rC14, rC14
               xsnegdp rC24, rC24
               xsnegdp rC34, rC34
               xsnegdp rC44, rC44
               xsnegdp rC54, rC54
               xsnegdp rC05, rC05
               xsnegdp rC15, rC15
               xsnegdp rC25, rC25
               xsnegdp rC35, rC35
               xsnegdp rC45, rC45
               xsnegdp rC55, rC55
            #endif
         #endif
         lxsd rA0, 0(pA)
         lxsd rA1, 8(pA)
         lxsd rA2, 8*2(pA)
         lxsd rA3, 8*3(pA)
         lxsd rA4, 8*4(pA)
         lxsd rA5, 8*5(pA)

         #if KB > 1
         lxsd ra0, 8*6(pA)
         lxsd ra1, 8*7(pA)
         lxsd ra2, 8*8(pA)
         lxsd ra3, 8*9(pA)
         lxsd ra4, 8*10(pA)
         lxsd ra5, 8*11(pA)
         #endif

         lxsd rB0, 0(pB)
         lxsd rB1, 8(pB)
         lxsd rB2, 8*2(pB)
         lxsd rB3, 8*3(pB)
         lxsd rB4, 8*4(pB)
         lxsd rB5, 8*5(pB)

         #if KB > 1
         lxsd rb0, 8*6(pB)
         lxsd rb1, 8*7(pB)
         lxsd rb2, 8*8(pB)
         lxsd rb3, 8*9(pB)
         lxsd rb4, 8*10(pB)
         lxsd rb5, 8*11(pB)
         #endif

         #if KB >= 1
            xsmaddadp rC00, rA0, rB0
            xsmaddadp rC10, rA1, rB0
            xsmaddadp rC20, rA2, rB0
            xsmaddadp rC30, rA3, rB0
            xsmaddadp rC40, rA4, rB0
            xsmaddadp rC50, rA5, rB0
            #if KB > 1
               lxsd rB0, 96(pB)
            #endif

            xsmaddadp rC01, rA0, rB1
            xsmaddadp rC11, rA1, rB1
            xsmaddadp rC21, rA2, rB1
            xsmaddadp rC31, rA3, rB1
            xsmaddadp rC41, rA4, rB1
            xsmaddadp rC51, rA5, rB1
            #if KB > 1
               lxsd rB1, 104(pB)
            #endif

            xsmaddadp rC02, rA0, rB2
            xsmaddadp rC12, rA1, rB2
            xsmaddadp rC22, rA2, rB2
            xsmaddadp rC32, rA3, rB2
            xsmaddadp rC42, rA4, rB2
            xsmaddadp rC52, rA5, rB2
            #if KB > 1
               lxsd rB2, 112(pB)
            #endif

            xsmaddadp rC03, rA0, rB3
            xsmaddadp rC13, rA1, rB3
            xsmaddadp rC23, rA2, rB3
            xsmaddadp rC33, rA3, rB3
            xsmaddadp rC43, rA4, rB3
            xsmaddadp rC53, rA5, rB3
            #if KB > 1
               lxsd rB3, 120(pB)
            #endif

            xsmaddadp rC04, rA0, rB4
            xsmaddadp rC14, rA1, rB4
            xsmaddadp rC24, rA2, rB4
            xsmaddadp rC34, rA3, rB4
            xsmaddadp rC44, rA4, rB4
            xsmaddadp rC54, rA5, rB4
            #if KB > 1
               lxsd rB4, 128(pB)
            #endif

            xsmaddadp rC05, rA0, rB5
            #if KB > 1
               lxsd rA0, 96(pA)
            #endif
            xsmaddadp rC15, rA1, rB5
            #if KB > 1
               lxsd rA1, 104(pA)
            #endif
            xsmaddadp rC25, rA2, rB5
            #if KB > 1
               lxsd rA2, 112(pA)
            #endif
            xsmaddadp rC35, rA3, rB5
            #if KB > 1
               lxsd rA3, 120(pA)
            #endif
            xsmaddadp rC45, rA4, rB5
            #if KB > 1
                lxsd rA4, 128(pA)
            #endif
            xsmaddadp rC55, rA5, rB5
         #endif

         #if KB >= 2
               lxsd rA5, 136(pA)
               lxsd rB5, 136(pB)

            xsmaddadp rC00, ra0, rb0
            xsmaddadp rC10, ra1, rb0
            xsmaddadp rC20, ra2, rb0
            xsmaddadp rC30, ra3, rb0
            xsmaddadp rC40, ra4, rb0
            xsmaddadp rC50, ra5, rb0
            #if KB > 2
               lxsd rb0, 144(pB)
            #endif

            xsmaddadp rC01, ra0, rb1
            xsmaddadp rC11, ra1, rb1
            xsmaddadp rC21, ra2, rb1
            xsmaddadp rC31, ra3, rb1
            xsmaddadp rC41, ra4, rb1
            xsmaddadp rC51, ra5, rb1
            #if KB > 2
               lxsd rb1, 152(pB)
            #endif
            xsmaddadp rC02, ra0, rb2
            xsmaddadp rC12, ra1, rb2
            xsmaddadp rC22, ra2, rb2
            xsmaddadp rC32, ra3, rb2
            xsmaddadp rC42, ra4, rb2
            xsmaddadp rC52, ra5, rb2
            #if KB > 2
               lxsd rb2, 160(pB)
            #endif

            xsmaddadp rC03, ra0, rb3
            xsmaddadp rC13, ra1, rb3
            xsmaddadp rC23, ra2, rb3
            xsmaddadp rC33, ra3, rb3
            xsmaddadp rC43, ra4, rb3
            xsmaddadp rC53, ra5, rb3
            #if KB > 2
               lxsd rb3, 168(pB)
            #endif
            xsmaddadp rC04, ra0, rb4
            xsmaddadp rC14, ra1, rb4
            xsmaddadp rC24, ra2, rb4
            xsmaddadp rC34, ra3, rb4
            xsmaddadp rC44, ra4, rb4
            xsmaddadp rC54, ra5, rb4
            #if KB > 2
               lxsd rb4, 176(pB)
            #endif

            xsmaddadp rC05, ra0, rb5
            #if KB > 2
               lxsd ra0, 144(pA)
            #endif
            xsmaddadp rC15, ra1, rb5
            #if KB > 2
               lxsd ra1, 152(pA)
            #endif
            xsmaddadp rC25, ra2, rb5
            #if KB > 2
               lxsd ra2, 160(pA)
            #endif
            xsmaddadp rC35, ra3, rb5
            #if KB > 2
               lxsd ra3, 168(pA)
            #endif
            xsmaddadp rC45, ra4, rb5
            #if KB > 2
               lxsd ra4, 176(pA)
            #endif
            xsmaddadp rC55, ra5, rb5
            #if KB > 2
               lxsd ra5, 184(pA)
               lxsd rb5, 184(pB)
            #endif
         #endif
         #if KB >= 3
            xsmaddadp rC00, rA0, rB0
            xsmaddadp rC10, rA1, rB0
            xsmaddadp rC20, rA2, rB0
            xsmaddadp rC30, rA3, rB0
            xsmaddadp rC40, rA4, rB0
            xsmaddadp rC50, rA5, rB0
            #if KB > 3
               lxsd rB0, 192(pB)
            #endif

            xsmaddadp rC01, rA0, rB1
            xsmaddadp rC11, rA1, rB1
            xsmaddadp rC21, rA2, rB1
            xsmaddadp rC31, rA3, rB1
            xsmaddadp rC41, rA4, rB1
            xsmaddadp rC51, rA5, rB1
            #if KB > 3
               lxsd rB1, 200(pB)
            #endif

            xsmaddadp rC02, rA0, rB2
            xsmaddadp rC12, rA1, rB2
            xsmaddadp rC22, rA2, rB2
            xsmaddadp rC32, rA3, rB2
            xsmaddadp rC42, rA4, rB2
            xsmaddadp rC52, rA5, rB2
            #if KB > 3
               lxsd rB2, 208(pB)
            #endif

            xsmaddadp rC03, rA0, rB3
            xsmaddadp rC13, rA1, rB3
            xsmaddadp rC23, rA2, rB3
            xsmaddadp rC33, rA3, rB3
            xsmaddadp rC43, rA4, rB3
            xsmaddadp rC53, rA5, rB3
            #if KB > 3
               lxsd rB3, 216(pB)
            #endif

            xsmaddadp rC04, rA0, rB4
            xsmaddadp rC14, rA1, rB4
            xsmaddadp rC24, rA2, rB4
            xsmaddadp rC34, rA3, rB4
            xsmaddadp rC44, rA4, rB4
            xsmaddadp rC54, rA5, rB4
            #if KB > 3
               lxsd rB4, 224(pB)
            #endif

            xsmaddadp rC05, rA0, rB5
            #if KB > 3
               lxsd rA0, 192(pA)
            #endif
            xsmaddadp rC15, rA1, rB5
            #if KB > 3
               lxsd rA1, 200(pA)
            #endif
            xsmaddadp rC25, rA2, rB5
            #if KB > 3
               lxsd rA2, 208(pA)
            #endif
            xsmaddadp rC35, rA3, rB5
            #if KB > 3
               lxsd rA3, 216(pA)
            #endif
            xsmaddadp rC45, rA4, rB5
            #if KB > 3
                lxsd rA4, 224(pA)
            #endif
            xsmaddadp rC55, rA5, rB5
         #endif

         #if KB >= 4
               lxsd rA5, 232(pA)
               lxsd rB5, 232(pB)

            xsmaddadp rC00, ra0, rb0
            xsmaddadp rC10, ra1, rb0
            xsmaddadp rC20, ra2, rb0
            xsmaddadp rC30, ra3, rb0
            xsmaddadp rC40, ra4, rb0
            xsmaddadp rC50, ra5, rb0
            #if KB > 4
               lxsd rb0, 240(pB)
            #endif

            xsmaddadp rC01, ra0, rb1
            xsmaddadp rC11, ra1, rb1
            xsmaddadp rC21, ra2, rb1
            xsmaddadp rC31, ra3, rb1
            xsmaddadp rC41, ra4, rb1
            xsmaddadp rC51, ra5, rb1
            #if KB > 4
               lxsd rb1, 248(pB)
            #endif
            xsmaddadp rC02, ra0, rb2
            xsmaddadp rC12, ra1, rb2
            xsmaddadp rC22, ra2, rb2
            xsmaddadp rC32, ra3, rb2
            xsmaddadp rC42, ra4, rb2
            xsmaddadp rC52, ra5, rb2
            #if KB > 4
               lxsd rb2, 256(pB)
            #endif

            xsmaddadp rC03, ra0, rb3
            xsmaddadp rC13, ra1, rb3
            xsmaddadp rC23, ra2, rb3
            xsmaddadp rC33, ra3, rb3
            xsmaddadp rC43, ra4, rb3
            xsmaddadp rC53, ra5, rb3
            #if KB > 4
               lxsd rb3, 264(pB)
            #endif
            xsmaddadp rC04, ra0, rb4
            xsmaddadp rC14, ra1, rb4
            xsmaddadp rC24, ra2, rb4
            xsmaddadp rC34, ra3, rb4
            xsmaddadp rC44, ra4, rb4
            xsmaddadp rC54, ra5, rb4
            #if KB > 4
               lxsd rb4, 272(pB)
            #endif

            xsmaddadp rC05, ra0, rb5
            #if KB > 4
               lxsd ra0, 240(pA)
            #endif
            xsmaddadp rC15, ra1, rb5
            #if KB > 4
               lxsd ra1, 248(pA)
            #endif
            xsmaddadp rC25, ra2, rb5
            #if KB > 4
               lxsd ra2, 256(pA)
            #endif
            xsmaddadp rC35, ra3, rb5
            #if KB > 4
               lxsd ra3, 264(pA)
            #endif
            xsmaddadp rC45, ra4, rb5
            #if KB > 4
               lxsd ra4, 272(pA)
            #endif
            xsmaddadp rC55, ra5, rb5
            #if KB > 4
               lxsd ra5, 280(pA)
               lxsd rb5, 280(pB)
            #endif
         #endif
/*
 *       Write out to C
 */
         stxsd rC00, 0(pC)

         addi pC, pC, 8*36
         addi pB, pB, KB*6*8    /* pB += KB*nu*sizeof */
         addic. NNU, NNU, -1
      bne+ NLOOP
      addi pA, pA, KB*6*8       /* pA += KB*mu*sizeof */
      addic. NMU, NMU, -1
   bne+ MLOOP
/*
 * Restore fpregs
 */
   lfd f14, FOFF(r1)
   lfd f15, 8+FOFF(r1)
   lfd f16, 16+FOFF(r1)
   lfd f17, 24+FOFF(r1)
   lfd f18, 32+FOFF(r1)
   lfd f19, 40+FOFF(r1)
   lfd f20, 48+FOFF(r1)
   lfd f21, 56+FOFF(r1)
   lfd f22, 64+FOFF(r1)
   lfd f23, 72+FOFF(r1)
   lfd f24, 80+FOFF(r1)
   lfd f25, 88+FOFF(r1)
   lfd f26, 96+FOFF(r1)
   lfd f27, 104+FOFF(r1)
   lfd f28, 112+FOFF(r1)
   lfd f29, 120+FOFF(r1)
   lfd f30, 128+FOFF(r1)
   lfd f31, 136+FOFF(r1)
/*
 * Restore vec regs
 */
   addi pB, r1, VOFF
   lxv v52, 0(pB)
   lxv v53, 16*1(pB)
   lxv v54, 16*2(pB)
   lxv v55, 16*3(pB)
   lxv v56, 16*4(pB)
   lxv v57, 16*5(pB)
   lxv v58, 16*6(pB)
   lxv v59, 16*7(pB)
   #if 0
   lxv v60, 16*8(pB)
   lxv v61, 16*9(pB)
   lxv v62, 16*10(pB)
   lxv v63, 16*11(pB)
   #endif
/*
 * Restore iregs
 */
   ld r2, 24(r1)        /* TOC reg always at 24(sp) */
/*   ld r14, IOFF(r1) */
/*
 * Dealloc stack frame & return
 */
   addi r1, r1, FSIZE
   blr

