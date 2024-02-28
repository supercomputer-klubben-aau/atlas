/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */


#if defined(ATL_ARCH_TI_C66_BM) && (!defined(ATL_USING_XCC))
   /* The following code is used because ATL_ARCH_TI_C66_BM is defined.    */
   /* Note: timer overhead has been measured; it is 2 cycles. TSCH, TSCL   */
   /* are C66 registers we access directly; TSCL=0 restarts timer. It is   */
   /* C66 specific; it cannot be used when compiling for the host, e.g.    */
   /* ATL_Xwalltime. Also, their time.h file has the wrong CLOCKS_PER_SEC. */
   #include <c6x.h>
   double ATL_cputime(void)
   {
      static int INIT=0;                  /* First time in, reset timer. */
      long long unsigned int now;
      static const double CPS = 1.0 / (1000000000.0);
      double d;

      if (INIT==0)                        /* Reset timer if first run. */
      {
         INIT = 1;                        /* Remember we did it. */
         TSCL = 0;                        /* Ensure hardware time is running. */
      }

      now = _itoll(TSCH, TSCL);           /* Convert timer regs to long long. */
      d = (double) (now);                 /* get as a double. */
      d *= CPS;                           /* Convert to seconds. */
      return(d);                          /* Exit with answer.  */
   } /* END ATL_cputime */
#elif defined(UseClock)
   #include <time.h>
   double ATL_cputime(void)
   {
      clock_t t1;
      static int INIT=0;
      static clock_t t0;
      static const double CPS = 1.0 / (1.0*CLOCKS_PER_SEC);
      double d;

      if (INIT)
      {
         t1 = clock() - t0;
         d = t1 * CPS;
         return(d);
      }
      INIT = 1;
      t0 = clock();
      return(0.0);
   }
#elif defined(POSIX_HR) /* use the POSIX HR timers */
   #include <time.h>
   double ATL_cputime(void)
   {
      struct timespec ts;
      static double t0;
      double res;
      static int INIT = 0;

      if (INIT)
      {
         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
         res = ts.tv_sec + 1.0e-9 * ts.tv_nsec;
         return(res - t0);
      }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&ts);
      t0 = ts.tv_sec + 1.0e-9 * ts.tv_nsec;
      INIT = 1;
      return(0.0);
   }
#elif defined(UseTimes)
   #include <stdlib.h>
   #include <sys/times.h>
   #include <unistd.h>
   double ATL_cputime(void)
   {
      struct tms ts;
      static double ClockTick=0.0;

      if (ClockTick == 0.0) ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
      times(&ts);
      return( ((double) ts.tms_utime) * ClockTick );
   }
#elif defined(SUN_HR) /* use sun high resolution timers */
   #include <sys/time.h>
   double ATL_cputime(void)
   {
      return(gethrvtime()*1.0e-9);
   }
#else
   #include <stdlib.h>
   #include <sys/time.h>
   #include <sys/resource.h>
   double ATL_cputime(void)
   {
      struct rusage ruse;
      getrusage(RUSAGE_SELF, &ruse);
      return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec/1000000.0) );
   }
#endif

