/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bitvec.h"
void ATL_print1dBV(int col, int N, void *bv)
{
   int i, ierr;
   printf("\n\nBITVEC MAP:\n");
   ierr = ATL_FindFirstSetBitBV(bv, 0);
   if (ierr == -1)
      ierr = N;
   if (col)
   {
      for (i=0; i < N; i++)
      {
         if (i != ierr)
            printf("%d .\n", i%10);
         else
         {
            printf("%d X\n", i%10);
            if (i < N-1)
               ierr = ATL_FindFirstSetBitBV(bv, i+1);
         }
      }
   }
   else
   {
      for (i=0; i < N; i++)
         printf("%d", i%10);
      printf("\n");
      for (i=0; i < N; i++)
      {
         if (i != ierr)
            printf(".");
         else
         {
            printf("X");
            if (i < N-1)
               ierr = ATL_FindFirstSetBitBV(bv, i+1);
         }
      }
      printf("\n");
   }
}
