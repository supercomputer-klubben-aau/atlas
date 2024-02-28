#include "atlas_threads.h"
#include "atlas_misc.h"
#include "atlas_cbc.h"
#if defined(TCPLX) && (defined(COMBMIN) || defined(COMBMAX))
   #include <math.h>  /* for fabs */
#endif
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
#if defined(COMBMIN)
   #define COMB Mjoin(PATL,comb_min)
   #define combvals(v1_, v2_) ((v1_) <= (v2_)) ? (v1_) : (v2_)
#elif defined(COMBMAX)
   #define COMB Mjoin(PATL,comb_max)
   #define combvals(v1_, v2_) ((v1_) >= (v2_)) ? (v1_) : (v2_)
#elif defined(COMBSUM)
   #define COMB Mjoin(PATL,comb_sum)
   #define combvals(v1_, v2_) (v1_) + (v2_)
#else
   #error "Unknown combine!"
#endif
#ifdef TCPLX
   void COMB
   (
      ATL_CUINT P,     /* # of threads in combine */
      ATL_CUINT iam,   /* rank of calling thread in combine */
      TYPE *val,       /* local min */
      void *vchk
   )
#else
   TYPE COMB
   (
      ATL_CUINT P,     /* # of threads in combine */
      ATL_CUINT iam,   /* rank of calling thread in combine */
      TYPE val,        /* local min */
      void *vchk
   )
#endif
{
   volatile char *bchk = vchk ? (volatile char*) vchk : ATL_TP_PTR->bchkin;
   volatile char *mybool = bchk + (iam<<ATL_chksh);
   volatile TYPE *myval = (volatile TYPE*)(mybool + ATL_sizeof);
   const char newv = !(*mybool);

   if (iam)
   {
      volatile TYPE *ans = (volatile TYPE*) (bchk + ATL_sizeof);
      #ifdef TCPLX
         *myval = *val;
         myval[1] = val[1];
      #else
         *myval = val;
      #endif
      *mybool = newv;
      while (*bchk != newv);
      #ifdef TCPLX
         *val = *ans;
         val[1] = ans[1];
      #else
         val = *ans;
      #endif
   }
   else
   {
      int i;
      #if defined(TCPLX) && !defined(COMBSUM)
         TYPE mv = fabs(*val) + fabs(val[1]);
      #endif
      for (i=1; i < P; i++)
      {
         ATL_CUINT d = i<<ATL_chksh;
         volatile TYPE *hisval = (volatile TYPE*)(bchk+d+ATL_sizeof);
         TYPE hv;

         while (bchk[d] != newv);  /* wait for his answer to appear */
         #ifdef TCPLX
            #if defined(COMBSUM)
               *val += *hisval;
               val[1] += hisval[1];
            #else
               hv = fabs(*hisval) + fabs(hisval[1]);
               #if defined(COMBMAX)
               if (hv > mv)
               #else /* COMBMIN */
               if (hv < mv)
               #endif
               {
                  mv = hv;
                  *val = *hisval;
                  val[1] = hisval[1];
               }
            #endif
         #else
            hv = *hisval;
            val = combvals(val, hv);
         #endif
      }
      #ifdef TCPLX
         *myval = *val;      /* provide global answer for real */
         myval[1] = val[1];  /* and for imaginary */
      #else
         *myval = val;       /* provide global answer */
      #endif
      *bchk = newv;          /* signal answer is ready */
   }
   #ifndef TCPLX
      return(val);
   #endif
}
