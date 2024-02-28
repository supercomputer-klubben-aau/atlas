/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */

#include "atlas_misc.h"
#include "atlas_level1.h"

void Mjoin(PATL,rotm)(const int N, TYPE *X, const int incX,
                      TYPE *Y, const int incY, const TYPE *P)
{
   int i;
   const TYPE flag = *P;
   TYPE h11, h21, h12, h22, w, z;

   if (N <= 0 || flag == -2.0) return;
   if (flag == ATL_rnone)
   {
      h11=P[1]; h21=P[2]; h12=P[3]; h22=P[4];
      if (incX == 1 && incY == 1)
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = w * h11 + z * h12;
            *Y++ = w * h21 + z * h22;
         }
      }
      else
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = w * h11 + z * h12;
            *Y = w * h21 + z * h22;
            X += incX;
            Y += incY;
         }
      }
   }
   else if (flag == ATL_rzero)
   {
      h21=P[2];
      h12=P[3];
      if (incX == 1 && incY == 1)
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = w + z * h12;
            *Y++ = w * h21 + z;
         }
      }
      else
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = w + z * h12;
            *Y = w * h21 + z;
            X += incX;
            Y += incY;
         }
      }
   }
   else if (flag == ATL_rone)
   {
      h11=P[1];
      h22=P[4];
      if (incX == 1 && incY == 1)
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X++ = w * h11 + z;
            *Y++ = z * h22 - w;
         }
      }
      else
      {
         for (i=N; i; i--)
         {
            w = *X;
            z = *Y;
            *X = w * h11 + z;
            *Y = z * h22 - w;
            X += incX;
            Y += incY;
         }
      }
   }
}
