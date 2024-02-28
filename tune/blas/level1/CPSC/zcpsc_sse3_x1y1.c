#include "atlas_misc.h"
#include "atlas_prefetch.h"
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>

void ATL_UCPSC(const int N, const TYPE *alpha, const TYPE *pX, const int incx,
               TYPE *pY, const int incy)
{
   const __m128d vALn={*alpha, alpha[1]};  /* vALn={ia, ra} */
   const __m128d vALs={-alpha[1], *alpha}; /* vALn={ra,-ia} */
   ATL_CSZT N2=N+N;
   ATL_SZT i;

   if ((((size_t)pY)&15L))
   {
      for (i=0; i != N2; i += 2)
      {
         __m128d rX, iX;
         rX = _mm_loaddup_pd(pX+i);        /* rX = {rX, rX} */
         rX = _mm_mul_pd(rX, vALn);        /* rX = {rX*ia, rX*ra} */
         iX = _mm_loaddup_pd(pX+i+1);      /* iX = {iX, iX} */
         iX = _mm_mul_pd(iX, vALs);        /* iX = {iX*ra,-iX*ia} */
         rX = _mm_add_pd(rX, iX);          /* rX={rX*ia+iX*ra, rX*ra--iX*ia} */
         _mm_storeu_pd(pY+i, rX);
      }
   }
   else  /* Y is aligned to 16-bytes */
   {
      for (i=0; i != N2; i += 2)
      {
         __m128d rX, iX;
         rX = _mm_loaddup_pd(pX+i);        /* rX = {rX, rX} */
         rX = _mm_mul_pd(rX, vALn);        /* rX = {rX*ia, rX*ra} */
         iX = _mm_loaddup_pd(pX+i+1);      /* iX = {iX, iX} */
         iX = _mm_mul_pd(iX, vALs);        /* iX = {iX*ra,-iX*ia} */
         rX = _mm_add_pd(rX, iX);          /* rX={rX*ia+iX*ra, rX*ra--iX*ia} */
         _mm_store_pd(pY+i, rX);
      }
   }
}
