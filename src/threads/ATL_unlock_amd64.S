#include "atlas_asm.h"
/*
 * 1 at memory means locked, 0 means unlocked
 *  ax                 4/rdi
 * int ATL_unlock(void *lck);
 */
.text
.global ATL_asmdecor(ATL_unlock)
ATL_asmdecor(ATL_unlock):
   xor %rax, %rax
   xchg %al, (%rdi)
   dec %al            /* signal error if it wasn't locked */
   ret
