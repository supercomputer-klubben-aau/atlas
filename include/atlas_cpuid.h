#ifndef ATLAS_CPUID_H
   #define  ATLAS_CPUID_H 1
/*
 * works like cpuid with these reg names.  Note that some inst take eax|ecx
 * as input, while others take only eax as input.  The number of output regs
 * also varies by query.
 */
#define do_cpuid(eax_, ebx_, ecx_, edx_)  __asm__("cpuid" \
           :"=a"(eax_), "=b"(ebx_), "=c"(ecx_), "=d"(edx_) \
           :"a"(eax_), "c"(ecx_) \
           : )
#endif
