#include "atlas_misc.h"
#include "atlas_simd.h"
#include "atlas_prefetch.h"
#define YDIST 96
void ATL_UAXPBY(const int N, const SCALAR alpha, const TYPE *X, const int incX,
                const SCALAR beta, TYPE *Y, const int incY)
{
   ATL_VTYPE valp, vbet;
   ATL_CSZT vlmskb = (ATL_VLENb-1);
   ATL_SZT i, nr, nv;
   for (i=0; i < N && ((size_t)(Y)&vlmskb); i++, Y++)  /* align Y */
      *Y = beta * *Y + alpha * X[i];
   nr = N - i;
   nv = nr & (~(ATL_VLEN-1));  /* get # of vec iterations remaining */
   if (nv)
   {
      nr -= nv;
      nv += i;
      ATL_vbcast(valp, &alpha);
      ATL_vbcast(vbet, &beta);
      if (((size_t)(X+i)&vlmskb) == 0)  /* Can we use aligned lds for X & Y? */
      {
         for (; i < nv; i += ATL_VLEN, Y += ATL_VLEN)
         {
            ATL_VTYPE x0, y0;
            ATL_vld(y0, Y);
            ATL_vmul(y0, y0, vbet);
            ATL_vld(x0, X+i);
            ATL_vmac(y0, valp, x0);
               ATL_pfl1W(Y+YDIST);
            ATL_vst(Y, y0);
         }
      }
      else  /* Y is aligned, X is not */
      {
         for (; i < nv; i += ATL_VLEN, Y += ATL_VLEN)
         {
            ATL_VTYPE x0, y0;
            ATL_vld(y0, Y);
            ATL_vmul(y0, y0, vbet);
            ATL_vuld(x0, X+i);
            ATL_vmac(y0, valp, x0);
               ATL_pfl1W(Y+YDIST);
            ATL_vst(Y, y0);
         }
      }
   }
   for (; i < N; i++, Y++)
      *Y = *Y * beta + alpha * X[i];
}
