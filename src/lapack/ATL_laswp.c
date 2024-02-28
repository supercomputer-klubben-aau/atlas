/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#include "atlas_lapack.h"

void ATL_laswp(const int N, TYPE *A, const int lda0, const int K1,
               const int K2, const int *piv, const int inci)
{
   #ifdef TCPLX
      const int lda = lda0<<1;
   #else
      #define lda lda0
   #endif
   const int n = K2 - K1;
   int nb = N >> 5;
   const int mr = N - (nb<<5);
   const int incA = lda << 5;
   const int *ipiv;
   int i, ip, i1, i2, KeepOn;
   register int h;
   TYPE *a0, *a1;
   #ifdef TCPLX
      register TYPE r0, r1;
   #else
      register TYPE r;
   #endif

   if (K2 < K1) return;
   if (inci < 0)
   {
      piv -= (K2-1) * inci;
      i1 = K2-1;
      i2 = K1;
   }
   else
   {
      piv += K1*inci;
      i1 = K1;
      i2 = K2-1;
   }

   if (nb)
   {
      do
      {
         ipiv = piv;
         i = i1;
         do
         {
            ip = *ipiv; ipiv += inci;
            if (ip != i)
            {
               a0 = A + (i SHIFT);
               a1 = A + (ip SHIFT);
               for (h=32; h; h--)
               {
                  #ifdef TCPLX
                     r0 = *a0;
                     r1 = a0[1];
                     *a0 = *a1;
                     a0[1] = a1[1];
                     *a1 = r0;
                     a1[1] = r1;
                  #else
                     r = *a0;
                     *a0 = *a1;
                     *a1 = r;
                  #endif
                  a0 += lda;
                  a1 += lda;
               }
            }
            if (inci > 0) KeepOn = (++i <= i2);
            else KeepOn = (--i >= i2);
         }
         while(KeepOn);
         A += incA;
      }
      while(--nb);
   }
   if (mr)
   {
      ipiv = piv;
      i = i1;
      do
      {
         ip = *ipiv; ipiv += inci;
         if (ip != i)
         {
            a0 = A + (i SHIFT);
            a1 = A + (ip SHIFT);
            for (h=mr; h; h--)
            {
               #ifdef TCPLX
                  r0 = *a0;
                  r1 = a0[1];
                  *a0 = *a1;
                  a0[1] = a1[1];
                  *a1 = r0;
                  a1[1] = r1;
               #else
                  r = *a0;
                  *a0 = *a1;
                  *a1 = r;
               #endif
               a0 += lda;
               a1 += lda;
            }
         }
         if (inci > 0) KeepOn = (++i <= i2);
         else KeepOn = (--i >= i2);
      }
      while(KeepOn);
   }
}
