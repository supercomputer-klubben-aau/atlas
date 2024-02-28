#include "atlas_bitvec.h"
int ATL_IncorpBV(ATL_BV_t *dst, ATL_BV_t *src, unsigned long pos)
/*
 * Copy entire contents of src into dst, starting at position pos
 * RETURNS: number of bits copied.
 */
{
   int iret=0;
   if (src)
   {
      ATL_BV_t ncpy=(*src), i;
      const ATL_BV_t nr=ncpy&modmskBV, pr=pos&modmskBV;
      ATL_assert(dst);
      ATL_assert(pos+ncpy <= *dst);
      iret = ncpy;
      ncpy >>= shBV;             /* # of full array entries to copy */
      src++;                     /* first src entry to copy */
      dst += 1 + (pos >> shBV);  /* first affected dest entry */
      if (pr)     /* usual case wt pos bits from old, to cancatonate */
      {           /* with 32-pos bits from src */
         unsigned char nnew = bpiBV - pr;
         const ATL_BV_t oldmsk = (1L<<pr)-1;
         ATL_BV_t last = *dst & oldmsk;
         for (i=0; i < ncpy; i++)
         {
            ATL_BV_t s = src[i];
            dst[i] = (s<<pr) | last;
            last = s >> nnew;
         }
         if (nr)  /* partial last block */
         {
            ATL_BV_t d = dst[ncpy], s=src[ncpy], msk=(1L<<nr)-1;
            int nextra = pr+nr - bpiBV;
            d = (d&(~oldmsk))|last;  /* new bits below pr, orig above */
            s = s & msk;             /* kill high bits not copied */
            msk = ~(msk<<pr);        /* nr 0s starting at pr */
            d = d & msk;             /* kill middle bits coming from s */
            d |= s << pr;
            dst[ncpy] = d;
            if (nextra > 0) /* nextra bits spill to next dest entry! */
            {
               int nhandled = nr-nextra;
               s >>= nr-nextra;
               msk = ~((1L<<nextra)-1);
               d = dst[ncpy+1];
               dst[ncpy+1] = (d&msk) | s;
            }
         }
         else /* write pr remaining bits to last dst */
            dst[ncpy] = (dst[ncpy]&(~oldmsk))|last;
      }
      else        /* best case, start on bpiBV boundary */
      {
         for (i=0; i < ncpy; i++)
            dst[i] = src[i];
         if (nr)  /* partial last block */
         {
            ATL_BV_t d = dst[ncpy], s=src[ncpy], msk=allsetBV<<nr;
            d = d & msk;     /* kill low bits overwritten by src */
            s = s & (~msk);  /* kill high bits not copied */
            dst[i] = d | s;
         }
      }
   }
   return(iret);
}
