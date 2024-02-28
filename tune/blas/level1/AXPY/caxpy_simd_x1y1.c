#include "atlas_misc.h"
#include "atlas_prefetch.h"
#include "atlas_cplxsimd.h"
#ifdef SCPLX
   #define MY_SIZEOF_REAL ATL_SSIZE
   #define MY_SIZEOF ATL_CSIZE
#else
   #define MY_SIZEOF_REAL ATL_DSIZE
   #define MY_SIZEOF ATL_ZSIZE
#endif
#define YDIST 96

void ATL_UAXPY(const int N, const TYPE *alpha, const TYPE *pX, const int incX,
               TYPE *pY, const int incY)
{
   int NV = N;
   int i, nr;
   ATL_VTYPE rALn, rALs;       /* natural order alpha, scaled/perm alpha */
   ATL_VTYPE rXr, rXi, rYn;

   ATL_vcxuld1(rALn, alpha);
   ATL_vcxPrepAlpha(rALn, rALs);        /* rALs = {ral, -ial}*(VLEN/2) */
                                        /* rALn = {ial,  ral}*(VLEN/2 */
/*
 * If pY has lower alignment than complex length, no peeling can fix Y align
 */
   if ( ((size_t)pY) & (ATL_sizeof-1) )
      goto YBADALIGN;
/*
 * i will be # of bytes past required alignment, peel to align Y
 */
   i = (int)(((size_t)pY) & (ATL_VLENb-1));
   i = ATL_DivBySize(i);  /* make it number of elements */
   i = (i) ? ATL_CXVLEN - i : 0;
   i = (i <= NV) ? i : NV;
/*
 * If Y is not aligned, peel i elts to make it aligned
 */
   if (i)
   {
      ATL_vcxlduXuYR(rYn, pY, rXi, pX, i);
      ATL_vcxsplitRI(rXr, rXi);
      ATL_vmac(rYn, rXi, rALs);
      ATL_vmac(rYn, rXr, rALn);
      ATL_vcxustR(pY, rYn, i);
      NV -= i;
      if (!NV)
         return;
      i += i;
      pX += i;
      pY += i;
   }
   nr = NV;
   NV >>= ATL_CXVLSH;  /* see how many rolled vec iterations we've got */
   nr -= NV<<ATL_CXVLSH;
   if (!NV)
      goto CLEANUP;
/*
 * If after peeling for Y alignment, X is mutually misaligned, go special loop
 */
   i = NV;
   #if ATL_CXSPLDb > MY_SIZEOF_REAL
      if ( ((size_t)pX) & (ATL_CXSPLDb-1) )
         goto XBADALIGN;
   #endif
/*
 * If we reach here, X&Y aligned, and there is at least 1 iter to do
 */
   do
   {
      ATL_vcxsplitRIld(rXr, rXi, pX);      /* rXr={rX,rX..}, rXi={iX,iX...} */
      ATL_vld(rYn, pY);                    /* rYn  = {iY,  rY  } */
      ATL_vmac(rYn, rXi, rALs);            /* rYn += {ral*iX, -ial*iX} */
      ATL_vmac(rYn, rXr, rALn);            /* rYn += {ial*rX, ral*rX} */
      ATL_vst(pY, rYn);
      pY += ATL_VLEN;
      pX += ATL_VLEN;
   }
   while(--i);
/*
 * Cleanup is written to do all ld/st unaligned so all vec loops can share
 * same cleanup.  ldXYR can cause code bloat if called a bunch of different
 * times, and this is an O(1) cost, so not worth it
 */
CLEANUP:
   if (nr)
   {
      #if ATL_CXSPLDb <= MY_SIZEOF_REAL
         ATL_vcxldXuYR(rXi, pX, rYn, pY, nr);
      #else
         ATL_vcxlduXuYR(rXi, pX, rYn, pY, nr);
      #endif
      ATL_vcxsplitRI(rXr, rXi);
      ATL_vmac(rYn, rXi, rALs);
      ATL_vmac(rYn, rXr, rALn);
      ATL_vcxustR(pY, rYn, nr);
   }
   return;
#if ATL_CXSPLDb > MY_SIZEOF_REAL
XBADALIGN:  /* Y is aligned, but X is not */
   do
   {
      ATL_vuld(rXi, pX);
      ATL_vcxsplitRI(rXr, rXi);            /* rXr={rX,rX..}, rXi={iX,iX...} */
      ATL_vld(rYn, pY);                    /* rYn  = {iY,  rY  } */
      ATL_vmac(rYn, rXi, rALs);            /* rYn += {ral*iX, -ial*iX} */
      ATL_vmac(rYn, rXr, rALn);            /* rYn += {ial*rX, ral*rX} */
      ATL_vst(pY, rYn);
      pY += ATL_VLEN;
      pX += ATL_VLEN;
   }
   while(--i);
   goto CLEANUP; /* use goto to avoid code bloat */
#endif

/*
 * Could attempt peel for X alignment, but makes code too nasty
 */
YBADALIGN:
   nr = NV;
   NV >>= ATL_CXVLSH;  /* see how many rolled vec iterations we've got */
   nr -= NV<<ATL_CXVLSH;
   for (i=NV; i; i--)
   {
      ATL_vuld(rYn, pY);                   /* rYn  = {iY,  rY  } */
      #if ATL_CXSPLDb <= MY_SIZEOF_REAL
         ATL_vcxsplitRIld(rXr, rXi, pX);   /* rXr={rX,rX..}, rXi={iX,iX...} */
      #else
         ATL_vuld(rXi, pX);
         ATL_vcxsplitRI(rXr, rXi);         /* rXr={rX,rX..}, rXi={iX,iX...} */
      #endif
                                           /* rALs = {ral, -ial} */
      ATL_vmac(rYn, rXi, rALs);            /* rYn += {ral*iX, -ial*iX} */

                                           /* rALn = {ial,  ral} */
      ATL_vmac(rYn, rXr, rALn);            /* rYn += {ial*rX, ral*rX} */
      ATL_vust(pY, rYn);
      pY += ATL_VLEN;
      pX += ATL_VLEN;
   }
   goto CLEANUP;
}
