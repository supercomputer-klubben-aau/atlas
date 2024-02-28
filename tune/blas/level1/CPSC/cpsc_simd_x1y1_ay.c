#include "atlas_misc.h"
#include "atlas_simd.h"
void ATL_UCPSC(const int N, const SCALAR alpha, const TYPE *X, const int incX,
               TYPE *Y, const int incY)
{
   int i, nr, nv;
   ATL_CSZT vlmskb = (ATL_VLENb-1);
/*
 * First, align Y
 */
   for (i=0; i < N && (((size_t)(Y+i))&vlmskb); i++)
      Y[i] = alpha * X[i];
   nr = N-i;                   /* num its remaining */
   nv = nr & (~(ATL_VLEN-1));  /* num of its we can do with vec code */
   if (nv)
   {
      ATL_VTYPE valp;
      nr -= nv;  /* nr now amount of post-vec cleanup */
      nv += i;
      ATL_vbcast(valp, &alpha);
      if (((size_t)(X+i)&vlmskb) == 0)  /* Are both X & Y aligned? */
      {
         for (; i < nv; i += ATL_VLEN)
         {
            ATL_VTYPE x0;
            ATL_vld(x0, X+i);
            ATL_vmul(x0, x0, valp);
            ATL_vst(Y+i, x0);
         }
      }
      else /* Y aligned, X known misaligned */
      {
         for (; i < nv; i += ATL_VLEN)
         {
            ATL_VTYPE x0;
            ATL_vuld(x0, X+i);
            ATL_vmul(x0, x0, valp);
            ATL_vst(Y+i, x0);
         }
      }
   }
   for (i=0; i < N; i++)
      Y[i] = alpha * X[i];
}
