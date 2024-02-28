/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 2007 R. Clint Whaley
 */
#include "atlas_misc.h"

#if ATL_LINEFLUSH    /* set in atlas-aux.h, if on one of the below archs */

#if defined(ATL_ARCH_PPCG5) || defined(ATL_ARCH_PPCG4) || defined(ATL_GAS_PPC)
   #define ATL_flushCacheLine(mem) __asm__ __volatile__ \
      ("dcbf 0, %0" :: "r"((void *)(mem)))
#elif defined(ATL_ARCH_IA64Itan) || defined(ATL_ARCH_IA64Itan2)
   #define ATL_flushCacheLine(mem) __asm__ __volatile__ \
      ("fc %0" :: "r"((void *)(mem)))
#elif defined(ATL_SSE2)
   #define ATL_flushCacheLine(mem) __asm__ __volatile__ \
      ("clflush %0" : : "m" (*((char *)(mem))))
#elif defined(ATL_ARCH_ARM64)  /* patch by Dave Nuechterlein */
   #define ATL_flushCacheLine(mem) __asm__ __volatile__ \
      ("dc civac, %[va]" \
         : \
         : [va] "r" ((void *)(mem)))

#else
   #define ATL_flushCacheLine(mem) \
   { \
      fprintf(stderr, "Cannot do cache-line flushing, %d of %s!\n",  \
              __LINE__, __FILE__); \
      exit(-1); \
   }
#endif

#if defined(ATL_ARCH_TI_C66_BM)       /* On the C66, just give cmd. */
   void ATL_flushCacheByAddr(size_t N, void *vp)
   {
      /*---------------------------------------------------------------------*/
      /* Address    Function                 Value                           */
      /* 0x01840000 Level 2 Cache Config     0x==== 0007 (all cache)         */
      /* 0x01845004 L2 Write-back invalidate 0x0000 0001 (flushes L2 cache)  */
      /* 0x01840040 Level 1 Cache Config     0x==== 0007 (all cache)         */
      /* 0x01845044 L1 Write-back invalidate 0x0000 0001 (flushes L1 cache)  */
      /*---------------------------------------------------------------------*/
      volatile unsigned int *L1_WB_INV = (unsigned int*) (0x01845044);
      volatile unsigned int *L2_WB_INV = (unsigned int*) (0x01845004);
      *L1_WB_INV = 1;                        /* Invalidate L1. */
      *L2_WB_INV = 1;                        /* Invalidate L2. */
   }
#else
   void ATL_flushCacheByAddr(size_t N, void *vp)
   {
      double *dp = vp;  /* assume cache line at least 8 bytes long */
      size_t i;
      for (i=0, N /= sizeof(double); i < N; i++)
         ATL_flushCacheLine(dp+i);
   }
#endif

FLSTRUCT *ATL_GetFlushStruct(void *p, size_t length, FLSTRUCT *next)
{
   FLSTRUCT *fp;

   fp = malloc(sizeof(FLSTRUCT));
   ATL_assert(fp);
   fp->p = p;
   fp->length = length;
   fp->next = next;

   return(fp);
}

void ATL_KillAllFlushStructs(FLSTRUCT *p)
{
   FLSTRUCT *kill;
   while (p)
   {
      kill = p;
      p = p->next;
      free(kill);
   }
}

void ATL_FlushAreasByCL(FLSTRUCT *fp)
{
   int i, n;
   char *cp;
   while (fp)
   {
      ATL_flushCacheByAddr(fp->length, fp->p);
      fp = fp->next;
   }
}

#endif   /* end #if ATL_LINEFLUSH */
