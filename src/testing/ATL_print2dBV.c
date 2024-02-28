/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bitvec.h"
void ATL_print2dBV(int M, int N, void *bv)
{
   int i, j, k, K=M*N, ierr = -1;

   printf("\n\nBITVEC MAP:\n");
   printf(" ");
   for (j=0; j < N; j++)
      printf("%d",j%10);
   printf("\n");
   ierr = ATL_FindFirstSetBitBV(bv, 0);
   if (ierr == -1)
      ierr = K;
   for (k=i=0; i < M; i++)
   {
      printf("%d", i%10);
      for (j=0; j < N; j++)
      {
         if (k != ierr)
            printf(".");
         else
         {
            printf("X");
            if (ierr+1 < K)
               ierr = ATL_FindFirstSetBitBV(bv, ierr+1);
            else
               ierr = K;
         }
         k++;
      }
      printf("\n");
   }
}

