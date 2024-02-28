/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ATL_WANT_ILCM
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
int itstLCM(unsigned int max)
{
   unsigned int i=1; /* if it'll work for short, it'll be OK */
   printf("Testing all lcm:\n");
   do
   {
      unsigned int j=1;
      printf("   testing all lcm with M=%u\n", i);
      do
      {
         int good, test;
         good = ATL_lcm(i, j);
         test = ATL_iLCM(i, j);
         if (test != good)
            printf("      LCM=(%u,%u): good=%u, bad=%u\n", i, j, good, test);
         assert(test == good);
      }
      while(++j < max);
      printf("   PASSED all lcm with M=%u\n", i);
   }
   while(++i != max);
   printf("PASSED all lcm!\n");
   return(0);
}
int itstPwr2(void)
{
   unsigned int i=1, pwr2=1, pwrp=0;
   unsigned char log;
   log = ATL_IsPow2(0);
   assert(log);
   do
   {
      unsigned int k;
      log = ATL_IsPow2(i);
      if (i == pwr2) /* really a power of two */
      {
         assert(log);
         ATL_iLeastSetBit(log, i);
         k = pwrp++;
         pwr2 += pwr2;
      }
      else
      {
         assert(!log);
         ATL_iLeastSetBit(log, i);
         for (k=0; (i&(1<<k)) == 0; k++);
      }
      if (k != log)
         printf("i=%x, k=%u, LSB=%u\n", i, k, log);
      assert (k == log);
   }
   while (++i);
   return(0);
}
int main(int nargs, char *args)
{
   assert(itstLCM(1<<16) == 0);
   printf("testing IsPwr2 & LeastSetBit\n");
   assert(itstPwr2() == 0);
   printf("SUCCESS!\n");
   return(0);
}
