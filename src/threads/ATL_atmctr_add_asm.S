#include "atlas_asm.h"
#include "atlas_asm.h"
#define ROUT ATL_asmdecor(ATL_atmctr_add)
#ifdef ATL_GAS_x8664
   .text
   ALIGN32
   .global ROUT
   ROUT:
      movq (%rdi), %rax       /* read cnt from memory */
      mov  %rax, %rcx         /* rcx = cnt */
      add %rsi, %rcx          /* rcx = val+cnt */
      lock                    /* make cmpxchg atomic */
      cmpxchg %rcx, (%rdi)    /* st newcnt if mem still == oldcnt in rax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ROUT                   /* ZF=0 means cmpxch failed, try again */
   ZERO_RET:
      xor %rax, %rax
   DONE:
   ret
#elif defined(ATL_GAS_WOW64)
   /* rax                      rcx                rdx  */
   /* int ATL_atmctr_add(void *vp, unsigned long val) */
   .text
   ALIGN32
   .global ROUT
   ROUT:
      movq (%rcx), %rax       /* read cnt from memory */
      add %rdx, %r8           /* r8 = val+cnt */
      lock                    /* make cmpxchg atomic */
      cmpxchg %r8, (%rcx)     /* st newcnt if mem still == oldcnt in rax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ROUT                   /* ZF=0 means cmpxch failed, try again */
   ZERO_RET:
      xor %rax, %rax
   DONE:
   ret
#elif defined(ATL_GAS_x8632)
   /* eax                  4(%esp)  8(%esp)  */
   /* int ATL_atmctr_add(void *vp, int val) */
   .text
   .global ROUT
   ROUT:
   push %esi
   movl 4+4(%esp), %edx       /* edx = ac */
   movl 4+8(%esp), %esi       /* esi = val */
   ATOMIC_LOOP:
      movl (%edx), %eax       /* read cnt from memory */
      movl %eax, %ecx         /* ecx = cnt */
      add %esi, %ecx          /* ecx = cnt + val */
      lock                    /* make cmpxchg atomic */
      cmpxchg %ecx, (%edx)    /* st newcnt if mem still == oldcnt in eax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ATOMIC_LOOP            /* ZF=0 means cmpxch failed, try again */
   pop %esi

   DONE:
   ret
#elif defined(ATL_GAS_ARM64)
/* x0                         x0                 x1 */
/* long ATL_atmctr_add(void *ac, unsigned long val) */
   .global ROUT
   ROUT:
   mov x3, x0            /* x3 = address of ctr */
   ATOMIC_LOOP:
      ldxr x0, [x3]      /* exclusive read of cnt into return reg (x0) */
      add x2, x0, x1     /* x2 = cnt + val */
      stxr w4, x2, [x3]  /* excl store newval, w4 0 on exclusive success */
   cbnz w4, ATOMIC_LOOP  /* if (w4 != 0) try again */
   ret
   ret
#endif
