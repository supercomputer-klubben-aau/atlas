#include "atlas_asm.h"
#include "atlas_asm.h"
#define ROUT ATL_asmdecor(ATL_atmctr_set)
#ifdef ATL_GAS_x8664
      mov %rsi, %rcx          /* rcx = val */
   /* rax                      rcx                rdx  */
   /* int ATL_atmctr_set(void *vp, unsigned long val) */
   .text
   ALIGN32
   .global ROUT
   ROUT:
      movq (%rcx), %rax       /* read cnt from memory */
   DONE:
   ret
#elif defined(ATL_GAS_x8632)
   /* eax                  4(%esp)  8(%esp)  */
   /* int ATL_atmctr_set(void *vp, int val) */
   .text
   .global ROUT
   ROUT:
   push %esi
   movl 4+4(%esp), %edx       /* edx = ac */
   movl 4+8(%esp), %esi       /* esi = val */
   ATOMIC_LOOP:
      movl (%edx), %eax       /* read cnt from memory */
      movl %eax, %ecx         /* ecx = cnt */
      mov %esi, %ecx          /* ecx = val */
      lock                    /* make cmpxchg atomic */
      cmpxchg %ecx, (%edx)    /* st newcnt if mem still == oldcnt in eax */
      je DONE                 /* ZF set if cmpxchg wrote to mem */
   jmp ATOMIC_LOOP            /* ZF=0 means cmpxch failed, try again */
   pop %esi

/* x0                         x0                 x1 */
/* long ATL_atmctr_set(void *ac, unsigned long val) */
   .global ROUT
   ROUT:
   mov x3, x0            /* x3 = address of ctr */
   ATOMIC_LOOP:
      ldxr x0, [x3]      /* exclusive read of cnt into return reg (x0) */
      stxr w4, x2, [x3]  /* excl store newval, w4 0 on exclusive success */
   cbnz w4, ATOMIC_LOOP  /* if (w4 != 0) try again */
   ret
   ret
#endif
