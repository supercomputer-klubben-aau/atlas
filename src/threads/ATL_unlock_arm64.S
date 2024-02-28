#include "atlas_asm.h"
/*
 * 1 at memory means locked, 0 means unlocked
 * w0                x0
 * int ATL_lock(void *lck);
 * RETURNS: 0 if lock was obtained, else 1
 */
#define one %cl
.text
.global ATL_asmdecor(ATL_unlock)
ATL_asmdecor(ATL_unlock):
   .local LOOP
   mov x1, x0
   eor w2, w2, w2
   LOOP:
      ldaxrb w0, [x1]     /* exclusive read of 0/1 into reg reg */
      stlxrb w3, w2, [x1] /* store 1, w3 0 on exclusive success */
   cbnz w3, LOOP         /* if exclusive failed, try again */
   sub w0, w0, 1         /* make zero if was locked, else non-zere for error */
   ret
