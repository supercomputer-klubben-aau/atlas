#include "atlas_misc.h"
#include "atlas_cplxsimd.h"

void ATL_UDOT(const int N, const TYPE *X, const int incx,
              const TYPE *Y, const int incy, SCALAR dot)
{
   ATL_VTYPE vdotR, vdotI, vX, vXr, vY;
   int NV=N, nr, i;

   ATL_vzero(vdotR);
   ATL_vzero(vdotI);
   if (N < ATL_CXVLEN)
   {
      i = N;
      ATL_vcxlduXuYR(vX, X, vY, Y, i);
      goto ONEANDDONE;
   }
/*
 * If neither ptr is alignable, to to unaligned code.  Otherwise, make X
 * be the alignable array, and peel to force X alignment
 */
   if ( ((size_t)X) & (ATL_sizeof-1L) )
   {
      const TYPE *p;
      if ( ((size_t)Y) & (ATL_sizeof-1L) )
         goto NOALIGN;
      p = X;  /* Y can be aligned, and X cannot */
      X = Y;  /* so set X ptr to alignable array */
      Y = p;  /* Y is unalignable in this case */
   }
/*
 * See if either array is already aligned, and let that one be X if so
 */
   i = (int)(((size_t)X) & (ATL_VLENb-1L));
   if (i)
   {
      int j;
      j = (int)(((size_t)Y) & (ATL_VLENb-1L));
      if (!j)  /* if Y is already aligned */
      {        /* accept Y as X/aligned, other vec unaligned */
         const TYPE *p=X;
         X = Y;
         Y = p;
         i = 0;  /* no peel necessary, since Y was aligned! */
      }
      else    /* if neither already aligned, see if either can be aligned */
      {
         int k;
         i = ATL_VLENb-i;  /* # of bytes until aligned ptr */
         j = ATL_VLENb-j;  /* # of bytes until aligned ptr */
         k = ATL_DivBySize(i);
         if (ATL_MulBySize(k) == i)    /* is X vector alignable? */
            i = k;                     /* set i to cplx elts to align X */
         else                          /* X unalignable, how about Y? */
         {
            k = ATL_DivBySize(j);
            if (ATL_MulBySize(k) == j) /* Is Y alignable? */
            {                          /* then make it X */
               const TYPE *p=X;
               X = Y;
               Y = p;
               i = k;
            }
            else                       /* neither vector is alignable! */
               goto NOALIGN;
         }
      }
   }
/*
 * If X is not aligned, but is alignable, peel to force alignment, or if
 * we don't have a full vect it, peel to fully handle the operation
 */
/*
 * If we reach here, X is alignable, i is # of elts to make it aligned
 */
   if (i)                                /* Do I need to peel i elts? */
   {
      ATL_vcxlduXuYR(vX, X, vY, Y, i);
      X += i+i;
      Y += i+i;
   }
   else  /* No align peel, peel full vec iteration to zero dot */
   {
      i = ATL_VLEN>>1;
      ATL_vld(vX, X);
      ATL_vuld(vY, Y);  /* unknown alignment at this point! */
      X += ATL_VLEN;
      Y += ATL_VLEN;
   }
ONEANDDONE:
   ATL_vcxswapRI(vXr, vX);
   ATL_vmul(vdotR, vX, vY);
   ATL_vmul(vdotI, vXr, vY);
   NV = N - i;
   if (!NV)
      goto VEC_REDUCE;
   nr = NV;
   NV >>= ATL_CXVLSH;
   nr -= NV<<ATL_CXVLSH;
   if (NV)  /* have some vector its left */
   {
      if ( ((size_t)Y) & (ATL_VLENb-1L) )  /* if Y not aligned */
          goto NOALIGNY;
      do
      {
         ATL_vld(vX, X);
         ATL_vcxswapRI(vXr, vX);
         ATL_vld(vY, Y);
         ATL_vmac(vdotR, vX, vY);
         ATL_vmac(vdotI, vXr, vY);
         X += ATL_VLEN;
         Y += ATL_VLEN;
      }
      while (--NV);
   }
   if (nr)  /* have a remainder */
   {
      ATL_vcxldXuYR(vX, X, vY, Y, nr);
      ATL_vcxswapRI(vXr, vX);
      ATL_vmac(vdotR, vX, vY);
      ATL_vmac(vdotI, vXr, vY);
   }
VEC_REDUCE:
   ATL_vcxdotcomb(vdotR, vdotI);
   ATL_vcxust1(dot, vdotR);
   return;
NOALIGNY:
   do
   {
      ATL_vld(vX, X);
      ATL_vcxswapRI(vXr, vX);
      ATL_vuld(vY, Y);
      ATL_vmac(vdotR, vX, vY);
      ATL_vmac(vdotI, vXr, vY);
      X += ATL_VLEN;
      Y += ATL_VLEN;
   }
   while (--NV);
   goto CLEANUP;
/*
 * Code where both pointers unalignable
 */
NOALIGN:
   nr = NV & ((ATL_VLEN>>1)-1);
   NV -= nr;
   for (i=0; i < NV; i += (ATL_VLEN>>1))
   {
      ATL_vuld(vX, X);
      ATL_vcxswapRI(vXr, vX);
      ATL_vuld(vY, Y);
      ATL_vmac(vdotR, vX, vY);
      ATL_vmac(vdotI, vXr, vY);
      X += ATL_VLEN;
      Y += ATL_VLEN;
   }
CLEANUP:
   if (nr)
   {
      ATL_vcxlduXuYR(vY, Y, vX, X, nr);
      ATL_vcxswapRI(vXr, vX);
      ATL_vmac(vdotR, vX, vY);
      ATL_vmac(vdotI, vXr, vY);
   }
   goto VEC_REDUCE;
}
