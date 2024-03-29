#ifndef ATLAS_TRSMNB_H
#define ATLAS_TRSMNB_H
#include "atlas_misc.h"
#include "atlas_kern3.h"

#ifndef TRSM_NB
   #ifdef TREAL
      #define TRSM_NB TRSM_Xover
/*
 *    For larger NB than 8, the intel compiler screws up the TRSM kernel,
 *    so force 8 as our largest stopping factor.  This is OK performance-
 *    wise, since all non-x86 archs benefit from not using the x86-specific
 *    kernel too much anyway (that's why we mandate 8 for all non-x86 archs)
 *    NOTE: I'm afraid the ATL_GAS_x86* probes might succeed on some
 *          IA64 due to emulation, and that's why they are explicit
 */
      #if defined(ATL_IntelIccBugs) || defined(ATL_ARCH_IA64Itan2) || \
          defined(ATL_ARCH_IA64Itan2) || \
          (!defined(ATL_GAS_x8632) && !defined(ATL_GAS_x8664))
         #undef TRSM_NB
         #define TRSM_NB 8
      #endif
   #else
      #define TRSM_NB 4
   #endif
#endif

#endif
