#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_cbc.h"
#include <math.h>  /* for fabs */
/*
 * Use cache-based comm to perform a scalar combine for P threads.
 * This code works on any system with coherent caches (weakly-ordered OK)
 * because the data and boolean sync variables are on the same cache line.
 * We guarantee this by separating each region by ATL_chkgap, which should
 * always be >= cache line size (default value 128).  Then, as long as
 * we can fit all data being combined into the same cache line, coherence
 * will guarantee we have the sync boolean and the data regardless of
 * weakly- or strongly-ordered caches.  At least on the ARM, however,
 * we must memory barrier to prevent OOE from advancing loads above the
 * sync.
 */
/*
 * This routine is used in LU, and it is different from comb_iamax in that
 * rank==0 does post its "we're done" to its sync array before exiting.
 * This means all threads except 0 remain blocked in this routine while 0
 * exits it with the global answer.  In this way we can keep other cores
 * waiting until after core 0 performs pivoting, which avoids doubling the
 * number of syncs required for unblkd LU.
 *
 * It is also different in that we return the global index using a pointer,
 * and the function return is the value to post to the sync array
 */
#ifdef DOFENCE
   #ifdef ATL_dcbc_comb_iamax_nopost0
      #undef ATL_dcbc_comb_iamax_nopost0
      #undef ATL_scbc_comb_iamax_nopost0
      #undef ATL_ccbc_comb_iamax_nopost0
      #undef ATL_zcbc_comb_iamax_nopost0
   #endif
   int Mjoin(PATL,cbc_comb_iamax_nopost0)/* RETURNS: global sync post value */
#else
   int Mjoin(PATL,comb_iamax_nopost0)    /* RETURNS: global sync post value */
#endif
(
   ATL_CUINT P,     /* # of threads in combine */
   ATL_CUINT iam,   /* rank of calling thread in combine */
   int *IDX,        /* INPUT: local index, OUTPUT: global idx */
   TYPE *val,       /* INPUT: local max, OUTPUT: global max */
   void *vchk
)
{
   volatile char *bchk = vchk ? (volatile char*) vchk : ATL_TP_PTR->bchkin;
   volatile char *mybool = bchk + (iam<<ATL_chksh);
   const int mysize = Mmax(ATL_isize, ATL_sizeof);
   volatile TYPE *myval = (volatile TYPE*)(mybool + mysize);
   volatile int *myidx = (volatile int*)(mybool + mysize + (mysize SHIFT));
   const char newv = !(*mybool);
   int idx = *IDX;

   if (iam)
   {
      volatile TYPE *ans = (volatile TYPE*) (bchk + mysize);
      volatile int *hisidx = (volatile int*) (bchk + mysize + (mysize SHIFT));
      *myidx = idx;
      #ifdef TCPLX
         myval[1] = val[1];
      #endif
      *myval = *val;
      #if defined(DOFENCE) && ATL_CBC_WEAK
         #if !ATL_CBC_NOBAR
            ATL_wmembar;
         #else
            ATL_mutex_lock(ATL_TP_PTR->cbcmut);   /* ahhhhhh */
            ATL_mutex_unlock(ATL_TP_PTR->cbcmut); /* the pain, the pain */
         #endif
      #endif
      *mybool = newv;
      while (*bchk != newv);
      *val = *ans;
      #ifdef TCPLX
         val[1] = ans[1];
      #endif
      idx = *hisidx;
      #if defined(DOFENCE) && ATL_CBC_WEAK
         #if !ATL_CBC_NOBAR
            ATL_rmembar;
         #else
            ATL_mutex_lock(ATL_TP_PTR->cbcmut);   /* ahhhhhh */
            ATL_mutex_unlock(ATL_TP_PTR->cbcmut); /* the pain, the pain */
         #endif
      #endif
   }
   else
   {
      int i;
      #if defined(TCPLX)
         TYPE mv = fabs(*val) + fabs(val[1]);
      #else
         TYPE mv = fabs(*val);
      #endif
      for (i=1; i < P; i++)
      {
         ATL_CUINT d = i<<ATL_chksh;
         volatile TYPE *hisval = (volatile TYPE*)(bchk+d+mysize);
         volatile int *hisidx = (volatile int*) (bchk+d+mysize+(mysize SHIFT));
         TYPE hv;
         int hidx;

         while (bchk[d] != newv);  /* wait for his answer to appear */
         hv = fabs(*hisval);
         hidx = *hisidx;
         #ifdef TCPLX
            hv += fabs(hisval[1]);
         #endif
         if (hv > mv || (hv == mv && hidx < idx))
         {
            idx = hidx;
            mv = hv;
            *val = *hisval;
            #ifdef TCPLX
               val[1] = hisval[1];
            #endif
         }
      }
      *myidx = idx;          /* provide global answer for index */
      *myval = *val;         /* and real value */
      #ifdef TCPLX
         myval[1] = val[1];  /* and for imaginary value */
      #endif
   }
   *IDX = idx;
   return(newv);
}
