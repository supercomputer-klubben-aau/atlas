/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_misc.h"

void Mjoin(PATL,set)(const int N, const TYPE alpha, TYPE *X, const int incX)
/*
 * X = alpha
 */
{
   int i, n;
   if (incX == 1)
   {
      n = N SHIFT;
      i = n >> 5;
      if (i)
      {
         n -= (i << 5);
         do
         {
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
            X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = X[16] = X[17] =
            X[18] = X[19] = X[20] = X[21] = X[22] = X[23] = X[24] = X[25] =
            X[26] = X[27] = X[28] = X[29] = X[30] = X[31] = alpha;
            X += 32;
         }
         while(--i);
      }
      if (n >> 4) /* >= 16 */
      {
         *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = X[8] = X[9] =
         X[10] = X[11] = X[12] = X[13] = X[14] = X[15] = alpha;
         X += 16;
         n -= 16;
      }
      if (n >> 3) /* >= 8 */
      {
         *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = X[7] = alpha;
         X += 8;
         n -= 8;
      }
      switch(n)
      {
         case 1:
            *X = alpha;
            break;
         case 2:
            *X = X[1] = alpha;
            break;
         case 3:
            *X = X[1] = X[2] = alpha;
            break;
         case 4:
            *X = X[1] = X[2] = X[3] = alpha;
            break;
         case 5:
            *X = X[1] = X[2] = X[3] = X[4] = alpha;
            break;
         case 6:
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = alpha;
            break;
         case 7:
            *X = X[1] = X[2] = X[3] = X[4] = X[5] = X[6] = alpha;
            break;
         default:;
      }
   }
   else
   {
      #ifdef TREAL
         for (i=N; i; i--, X += incX) *X = alpha;
      #else
         for (n=incX<<1, i=N; i; i--, X += n) *X = X[1] = alpha;
      #endif
   }
}
