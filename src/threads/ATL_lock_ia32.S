#include "atlas_asm.h"
/*
 * 1 at memory means locked, 0 means unlocked
 *                   4/rdi
 * void ATL_lock(void *lck);
 * RETURNS: 0 if lock was obtained, else 1
 */
#define one %cl
.text
.global ATL_asmdecor(ATL_lock)
ATL_asmdecor(ATL_lock):
   xor %eax, %eax
   mov $1, one
   .local LOOP
   LOOP:
      mov one, %al
      xchg %al, (%rdi)
      cmp one, %al      /* if (already locked) */
      je  WAIT          /*    kill some time & try again */
   ret
   .local FAIL
WAIT:
   pause
   jmp LOOP
