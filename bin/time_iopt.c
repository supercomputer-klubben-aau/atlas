/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ATL_WANT_ILCM
#define INLINE /* don't inline, so we can us funcptr */
#include "atlas_iopt.h"


/*
 * Original lcm written with collaboration between Antoine Petitet &
 * Clint Whaley based around Stein's Algorithm
 */
static int ATL_lcm(const unsigned int M, unsigned const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;
   unsigned long long MN;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else
            {
               MN = M*N;
               return( MN / (1<<gcd) );
            }
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      MN = M*N;
      return( MN / (max<<gcd) );
   }
   return(M);
}
double ATL_walltime(void);
#define cuint const unsigned int
double timeLCM(unsigned int max, int (*lcm)(cuint,cuint))
{
   double t0;
   unsigned int i=1;
   t0 = ATL_walltime();
   do
   {
      unsigned int j=1;
      do
      {
         int good;
         good = lcm(i, j);
         assert(good);     /* use ret so lcm call not dead */
      }
      while(++j != max);
   }
   while(++i != max);
   t0 = ATL_walltime() - t0;
   return(t0);
}
int main(int nargs, char **args)
{
   unsigned int N=512;
   double t0, tg;
   if ( nargs > 1)
      N = atoi(args[1]);
   t0 = timeLCM(N, ATL_lcm);
   tg = timeLCM(N, ATL_iLCM);
   printf("orig=%e, new=%e, speedup=%.2f\n", t0, tg, t0/tg);
   return(0);
}

