/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 1997, 2015 R. Clint Whaley
 */

/*
 * Cycle-accurate timers not quite thread safe, but all cores should
 * almost always be the same in upper 11 bits, which we need to mask
 * of to avoid losing precision going from long long -> double
 */
#ifdef PentiumCPS
   #include <assert.h>
   #ifdef __LP64__
      #define ATL_GetCycleCount(ll_) __asm__ __volatile__ \
      ("xorq %%rax, %%rax ; .byte 0x0f ; .byte 0x31 ; " \
       "shlq $32, %%rdx ; orq %%rdx, %%rax " \
       : "=a"(ll_)      /* output */ \
       :                /* no input params */ \
       : "%rdx"  /* clobber list */ \
      );
      #define MASK11 (((long long)0xFFE0)<<48)
      double ATL_walltime(void)
      {
         const static double mul=1.0e-6 / ((double)PentiumCPS);
         static long long top11bits=1, t11;
         long long cycnt;

         ATL_GetCycleCount(cycnt)
         t11 = cycnt & MASK11;
         if (top11bits != t11)
         {
            assert(top11bits == 1);
            top11bits = t11;
         }
         cycnt &= ~MASK11;
         return(mul*cycnt);
      }
   #else
      #define ATL_GetCycleCount(hi_, lo_) __asm__ __volatile__ \
      (".byte 0x0f ; .byte 0x31 ; " \
       : "=a"(lo_), "=d"(hi_)     /* output */ \
       :                          /* no input params */ \
       :                          /* clobber list */ \
      );
      #define MASK11 (0xFFE0<<16)
      double ATL_walltime(void)
      {
         const static double mul=1.0e-6 / ((double)PentiumCPS);
         static unsigned int top11bits=1, t11;
         unsigned long long cyc;
         unsigned int hi, lo;
         ATL_GetCycleCount(hi, lo);
         t11 = hi & MASK11;
         if (t11 != top11bits)
         {
            assert(top11bits == 1);
            top11bits = t11;
         }
         cyc = hi & (~MASK11);
         cyc = (cyc << 32) | lo;
         return(mul*cyc);
      }
   #endif
#elif defined(ATL_OS_WinNT) /* special code for windows */
   #include <windows.h>
   double ATL_walltime(void)
   {
      static double freqRecip = 0.0;
      LARGE_INTEGER msout;
      unsigned long long myout;
/*
 *    Not thread-safe, but shouldn't cause problems, since all cores
 *    should hopefully get the same answer
 */
      if (freqRecip == 0.0)
      {
         QueryPerformanceFrequency(&msout);
         myout = msout.HighPart;
         myout = (myout<<32) | msout.LowPart;
         freqRecip = 1.0/((double) myout);
      }
      QueryPerformanceCounter(&msout);
      myout = msout.HighPart;
      myout = (myout<<32) | msout.LowPart;
      return(myout*freqRecip);
   }
#elif defined(UseTimes)
   #include <stdlib.h>
   #include <sys/times.h>
   #include <unistd.h>
   double ATL_walltime(void)
   {
      struct tms ts;
      static double ClockTick=0.0;

      if (ClockTick == 0.0) ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
      return( ((double) times(&ts)) * ClockTick);
   }
#elif defined(SUN_HR) /* use sun high resolution timers */
   #include <sys/time.h>
   double ATL_walltime(void)
   {
      return(gethrtime()*1.0e-9);
   }
#elif defined(POSIX_HR) /* use the POSIX HR timers */
   #include <time.h>
   double ATL_walltime(void)
   {
      struct timespec ts;
      double res;

      clock_gettime(CLOCK_REALTIME, &ts);
      res = ts.tv_sec + 1.0e-9 * ts.tv_nsec;
      return(res);
   }
/*
 * Without gcc, I know no standard Windows wall-timer, so use cputime
 */
#elif (!defined(__GNUC__) && (defined(ATL_OS_Win9x) || defined(ATL_OS_WinNT))) \
      || defined(__MINGW32__)
   #include <time.h>
   double ATL_walltime(void)
   {
      static const double CPS = 1.0 / (1.0*CLOCKS_PER_SEC);

      return(clock() * CPS);
   }
#else
   #include <stdlib.h>
   #include <sys/time.h>
   #include <sys/resource.h>
   double ATL_walltime(void)
   {
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return( ((double) tp.tv_sec) + (tp.tv_usec*1.0e-6) );
   }
#endif

