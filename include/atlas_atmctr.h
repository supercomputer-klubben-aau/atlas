#ifndef ATLAS_ATMCTR_H
   #define ATLAS_ATMCTR_H
   #include "atlas_tprim.h"
/*
 * atmctr data structure has two variants.  The mutex version is of size
 * 2*SAFELS + sizeof(lock).  _new returns a SAFELS-aligned address, and the
 * data structure from this point is: [<cnt><off>][lock].
 *
 * For systems where we have true atomic ctrs, the lock is ommitted, so size
 * 2*SAFELS.
 *
 * If <off> is nonzero it is the bytes to add to ac to get to the
 * original malloc-returned ptr (used only in _free).
 */
   #ifndef ATL_ATM_ASM
      #if defined(ATL_GAS_x8632) || defined(ATL_GAS_x8664) || \
          defined(ATL_GAS_ARM64) || defined(ATL_GAS_WOW64)
         #define ATL_ATM_ASM 1  /* I've got assembly atomic ops */
      #else
         #define ATL_ATM_ASM 0  /* I do not have assembly atomic ops */
      #endif
   #endif

void *ATL_atmctr_init(void *vp, long cnt);
void *ATL_atmctr_new(long cnt);
void  ATL_atmctr_free(void *ac);
#if ! ATL_ATM_ASM
   #define ATL_atmctr_destroy(ac_) ATL_lock_destroy(ATL_IncBySafeLS(ac_))
   #define ATL_atmctr_sizeof (size_t) \
      ATL_AlignSafeLS(sizeof(ATL_lock_t)+ATL_SAFELS)
#else
   #define ATL_atmctr_sizeof ATL_SAFELS
   #define ATL_atmctr_destroy(ac_)
#endif
/* All these functions return prior value */
long ATL_atmctr_set(void *ac, long val);
long ATL_atmctr_dec(void *ac);
long ATL_atmctr_add(void *ac, unsigned long val);

#define ATL_atmctr_get(ac_) (*((volatile long*)(ac_)))
#endif
