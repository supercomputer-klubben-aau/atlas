#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ATL_DEF_RUNTHR
#include "atlas_tprim.h"

ATL_lock_t *GLCK=NULL;
unsigned long long GCNT=0;
unsigned long GN=80000;

void *TestLock(void *vp)
{
   ATL_thread_t *tp=vp;
   unsigned long i;

   for (i=0; i < GN; i++)
   {
      int k;
      k = ATL_lock(GLCK);
      if (k)
         printf("k=%d\n", k);
      assert(!k);
      GCNT++;
      assert(!ATL_unlock(GLCK));
   }
   return(NULL);
}

int main(int nargs, char **args)
{
   unsigned int P = ATL_NTHREADS;
   unsigned long long ans;
   if (nargs == 3 || nargs == 2)
   {
     GN = atol(args[1]);
     if (nargs == 3)
       P = atol(args[2]);
   }
   else if (nargs != 1)
   {
      printf("USAGE: %s [ntests] [nthreads]\n", args[0]);
      return(-1);
   }
   GLCK = malloc(sizeof(ATL_lock_t));
   printf("TESTING LOCK/UNLOCK WITH %lu TESTS AND %u THREADS.\n", GN, P);
   ATL_runThreads(ATL_NTHREADS, TestLock, NULL);
   ans = GN * P;
   if (ans == GCNT)
      printf("PASS.\n");
   else
      printf("FAIL: expected=%lu, got=%lu!\n", (unsigned long)ans,
             (unsigned long)GCNT);
   free(GLCK);
   return(GCNT != ans);
}
