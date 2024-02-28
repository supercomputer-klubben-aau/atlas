#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int main(int nargs, char **args)
{
   int i, iarr[5];
   if (nargs != 6)
   {
      /* entry in iarr:            0    1    2    3      4 */
      fprintf(stderr, "USAGE: %s <mb> <nb> <mu> <nu> <vlen>\n", args[0]);
      exit(1);
   }
   for (i=0; i < 5; i++)
   {
      int k;
      k = atoi(args[i+1]);
      assert (k > 0);
      if (i == 2 || i == 3)
         assert(k < iarr[i-2]);
      iarr[i] = k;
   }
   i = getsz(iarr[0],iarr[1],iarr[2],iarr[3],iarr[4]);
   printf("BLOCKSIZE=%d\n", i);
   return(0);
}
