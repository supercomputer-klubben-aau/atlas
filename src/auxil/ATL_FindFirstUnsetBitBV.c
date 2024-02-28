/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include "atlas_bitvec.h"
long ATL_FindFirstUnsetBitBV(ATL_BV_t *bv, unsigned long bs)
{
   ATL_BV_t ibv;
   if (bv)
   {
      ATL_BV_t n=bv[0], i;
      const ATL_BV_t allset=allsetBV, ln=n&modmskBV, lpos=bs&modmskBV;

      if (bs >= n)
         return(-1);
      n >>= shBV;            /* does not include any partial end block */
      bs >>= shBV;
      bv++;
      if (lpos)
      {
         ibv = bv[bs];
         ibv |= (1L<<lpos)-1L;  /* set skipped bits of first entry */
         if (bs == n)  /* 1st elt also last, maybe partial high bits! */
         {
            if (ln)
               ibv |= ~(allset>>(bpiBV-ln));
            if (ibv != allset)
               goto DONE_SRCH;
            return(-1);
         }
         else if (ibv != allset)
            goto DONE_SRCH;
         bs++;
      }
/*
 *    Look through all full integral elts
 */
      while (bs < n)
      {
         ibv = bv[bs];
         if (ibv != allset)
            goto DONE_SRCH;
         bs++;
      }
/*
 *    Search last, possibly partial, entry
 */
      if (ln)
      {
         ibv = bv[bs];
         ibv |= ~(allsetBV>>(bpiBV-ln));
         if (ibv != allset)
            goto DONE_SRCH;
      }
   }
   return(-1);
DONE_SRCH:  /* if we reach here, at least one bit is set! */
   bs <<= shBV;
   ibv = ~ibv;
   #if bpiBV == 64
   if (ibv & 0xFFFFFFFF)
      goto IN_LOW32;
   ibv >>= 32;
   bs += 32;
   IN_LOW32:
   #endif
   if (ibv & 0xFFFFL)
      goto IN_LOW16;
   ibv >>= 16;
   bs += 16;
   IN_LOW16:
   if (ibv & 0xFFL)
      goto IN_LOW8;
   ibv >>= 8;
   bs += 8;
   IN_LOW8:
   if (ibv & 0xFL)
      goto IN_LOW4;
   ibv >>= 4;
   bs += 4;
   IN_LOW4:
   if (ibv & 0x3L)
      goto IN_LOW2;
   ibv >>= 2;
   bs += 2;
   IN_LOW2:
      if (ibv&1L)
         return(bs);
   return(bs+1);
}
