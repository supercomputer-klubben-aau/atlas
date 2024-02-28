#include "atlas_asm.h"
#include "atlas_asm.h"
#define ROUT ATL_asmdecor(ATL_atmctr_dec)
#ifdef ATL_GAS_x8664
   /* rax                     rdi  */
   /* int ATL_atmctr_dec(void *vp) */
   .text
   ALIGN32
   .global ROUT
   ROUT:
      movq (%rdi), %rax       /* read cnt from memory */
      mov  %rax, %rcx         /* rcx = cnt */
      sub  $1, %rcx           /* rcx = cnt-1 */
      jl ZERO_RET             /* return 0 if count already below 1 */
      lock                    /* make cmpxchg atomic */
      cmpxchg %rcx, (%rdi)    /* st newcnt if mem still == oldcnt in rax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ROUT                   /* ZF=0 means cmpxch failed, try again */
   ZERO_RET:
      xor %rax, %rax
   DONE:
   ret
#elif defined(ATL_GAS_WOW64)
   /* rax                      rcx  */
   /* int ATL_atmctr_dec(void *vp) */
   .text
   ALIGN32
   .global ROUT
   ROUT:
      movq (%rcx), %rax       /* read cnt from memory */
      mov  %rax, %r8          /* r8 = cnt */
      subl $1, %r8            /* r8 = cnt-1 */
      jl ZERO_RET             /* return 0 if count already below 1 */
      add %rdx, %r8           /* r8 = val+cnt */
      lock                    /* make cmpxchg atomic */
      cmpxchg %r8, (%rcx)     /* st newcnt if mem still == oldcnt in rax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ROUT                   /* ZF=0 means cmpxch failed, try again */
   DONE:
   ret
#elif defined(ATL_GAS_x8632)
   /* eax                    4(%esp)  */
   /* int ATL_DecAtomicCount(void *vp) */
   .text
   .global ROUT
   ROUT:
   movl 4(%esp), %edx
   ATOMIC_LOOP:
      movl (%edx), %eax       /* read cnt from memory */
      movl %eax, %ecx         /* ecx = cnt */
      subl $1, %ecx           /* ecx = cnt-1 */
      jl ZERO_RET             /* return 0 if count already below 1 */
      lock                    /* make cmpxchg atomic */
      cmpxchg %ecx, (%edx)    /* put cnt-1 in mem if mem still == cnt in eax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ATOMIC_LOOP            /* ZF=0 means cmpxch failed, try again */

   ZERO_RET:
      xor %eax, %eax
   DONE:
   ret
#elif defined(ATL_GAS_ARM64)
/* x0                        x0
/* long ATL_atmctr_dec(void *ac) */
   .global ROUT
   ROUT:
   mov x3, x0            /* x3 = address of ctr */
   ATOMIC_LOOP:
      ldxr x0, [x3]      /* exclusive read of cnt into return reg (x0) */
      subs x2, x0, 1     /* dec cnt, set cond codes */
      b.lt   ZERO_RET    /* return 0 if count already below 1 */
      stxr w4, x2, [x3]  /* excl store newval, w4 0 on exclusive success */
   cbnz w4, ATOMIC_LOOP  /* if (w4 != 0) try again */
   ret
   ZERO_RET:
      eor w0, w0, w0
   ret
#endif
