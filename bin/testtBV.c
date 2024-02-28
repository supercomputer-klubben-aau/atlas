#include "atlas_tbitvec.h"
#include <assert.h>
#define ulong unsigned long

long FindGlobUnset(void *bv)
{
   long nunset=0, nb, i;
   nb = ATL_tGetTotBitsBV(bv);
   for (i=0; i < nb; i++)
      if (!ATL_tIsBitSetBV(bv, i))
         nunset++;
   return(nunset);
}

void setRange(void *bv, ulong b0, ulong mask)
{
   ulong N, i;
   N = ATL_tGetTotBitsBV(bv);
   if (b0 >= N)
      return;
   N -= b0;
   i = (sizeof(long)<<3);
   N = (N > i) ? i : N;
   for (i=0; i < N; i++)
   {
      if (mask & (1L<<i))
         ATL_tSetBitBV(bv, b0+i);
      else
         ATL_tUnsetBitBV(bv, b0+i);
   }
}

void testRange(ulong b0, ulong mask, unsigned int P)
{
   void *bv;
   unsigned long i;
   bv = ATL_tNewBV(b0+bpiBV, P);
   setRange(bv, b0, mask);
   for (i=0; i < bpiBV; i++)
   {
      int ans = (mask>>i)&1;
      assert(ATL_tSetBitBV(bv, b0+i) == ans);
   }
/*   assert(ATL_FindFirstBit() == -1); */
   ATL_tFreeBV(bv);
}

void testSetUnset(ulong len, unsigned int P)
{
   void *bv;
   unsigned long i;
   bv = ATL_tNewBV(len, P);
   for (i=0; i < len; i++)
   {
      unsigned long r;
      r  = ATL_tSetUnsetBitBV(bv, 0);
      if (r != i)
         fprintf(stderr, "   ERR: testSetUns, exp=%ld, got=%ld\n", i, r);
      assert(r == i);
   }
   assert(ATL_tSetUnsetBitBV(bv, 0) == -1);
   i = (len > 2) ? 2 : 0;
   assert(ATL_tUnsetBitBV(bv, i) == 1);
   assert(ATL_tSetUnsetBitBV(bv, 1) == i);
   assert(ATL_tSetUnsetBitBV(bv, 0) == -1);
   if (P > 1)
   {
      unsigned long beg1, end1;
      beg1 = ATL_tGetLocalBoundsBV(bv, 1, &end1);
      if (beg1 < len)
      {
         assert(ATL_tUnsetBitBV(bv, 0) == 1);
         assert(ATL_tUnsetBitBV(bv, beg1) == 1);
         assert(ATL_tSetUnsetBitBV(bv, 1) == beg1);
         assert(ATL_tSetUnsetBitBV(bv, 1) == 0);
         assert(ATL_tSetUnsetBitBV(bv, 0) == -1);
      }
   }
   ATL_tFreeBV(bv);
}

void testGlb2loc(ulong nbits, unsigned int P)
{
   void *gbv;
   ATL_BV_t *lbv;
   long i, pos;

   gbv = ATL_tNewBV(nbits, P);
   lbv = ATL_NewBV(nbits+bpiBV);
   srand48(nbits+(P<<3)+77);
   for (i=0; i < nbits; i += bpiBV)
      setRange(gbv, i, lrand48());

   for (pos=0; pos < bpiBV; pos++)
   {
      long nunset, un=0;
      nunset = ATL_tGlb2locBV(lbv, gbv, pos);
      for (i=0; i < nbits; i++)
      {
         int SET;
         SET = ATL_IsBitSetBV(lbv, pos+i);
         if (!SET)
            un++;
         if (SET != ATL_tIsBitSetBV(gbv, i))
         {
            fprintf(stderr, "Gbl2loc MISMATCH: nbits=%ld, pos=%ld, idx=%ld\n",
                    nbits, pos, i);
            ATL_assert(0);
         }
      }
      if (nunset != un)
         fprintf(stderr, "Gbl2loc NUNSET WRONG: exp=%ld, got=%ld\n",
                 un, nunset);

      ATL_assert(nunset == un);
   }
   ATL_tFreeBV(gbv);
   ATL_FreeBV(lbv);
}

int testSetRange0(void *gbv, ATL_BV_t *lbv, unsigned int nbits,
                  unsigned long pos, unsigned long setmsk)
{
   unsigned long oldmsk, omsk, nerr=0;
   unsigned int onbits=nbits, i;

   ATL_tGlb2locBV(lbv, gbv, 0);  /* local copy of gbv */
   oldmsk = ATL_tSetRangeBV(gbv, &nbits, pos, setmsk);
   if (!nbits && onbits)  /* made no change */
   {
      fprintf(stderr, "   NO BITS AT ALL!!\n");
      return(-1);
   }
   for (omsk=i=0; i < nbits; i++)
   {
      long v;
      if ((setmsk>>i)&1L)
         v = ATL_SetBitBV(lbv, pos+i);
      else
         v = ATL_UnsetBitBV(lbv, pos+i);
      omsk |= v<<i;
   }
   for (i=0; i < nbits; i++)
   {
       int exp, tst;
       exp = ATL_IsBitSetBV(lbv, pos+i);
       tst = ATL_tIsBitSetBV(gbv, pos+i);
       if (exp != tst)
       {
          fprintf(stderr,
                  "   tSetRng ERR: nbits=(%u,%u), pos0=%lu, i=%u, exp=%d\n",
                  nbits, onbits, pos, i, exp);
          nerr++;
       }
   }
   if (omsk != oldmsk)  /* didn't get same return value! */
   {
      fprintf(stderr, "    OLD ERR  exp=%x, got=%x\n", omsk, oldmsk);
      return(-2);       /* return err */
   }
   if (!nerr)
   {
      long *gp = ATL_AlignSafeLS(gbv);
      long tst, exp;
      exp = FindGlobUnset(gbv);
      tst = ATL_tInfoBV(gbv, ATL_TBV_NUNSET);
      if (tst != exp)
      {
         fprintf(stderr, "   NUNSET WRONG: exp=%ld, got=%ld\n", exp, tst);
         return(-3);
      }
   }
   return(nerr);
}

void testSetRange(ulong nbits, unsigned int P)
{
   void *gbv;
   ATL_BV_t *lbv;
   const unsigned long nelt=nbits>>shBV;
   unsigned int off, itst=0;

   srand48(nbits+(P<<3)+77);
   gbv = ATL_tNewBV(nbits, P);
   lbv = ATL_NewBV(nbits);
   for (off=0; off < bpiBV; off++)
   {
      unsigned int len;
      for (len=0; len < bpiBV; len++)
      {
         unsigned long setmsk, elt, pos;
         int err;

         if (len == bpiBV)
            setmsk = allsetBV;
         else
            setmsk = (1L<<len)-1;
         setmsk &= lrand48();
         elt = lrand48()%nelt;
         pos = (elt<<shBV)+off;
         if (pos >= nbits)
            pos -= nbits;
         err = testSetRange0(gbv, lbv, len, pos, setmsk);
         if (err)
         {
            fprintf(stderr, "   FAILED TEST %d: ret=%d, setmsk=%x\n",
                    itst, err, setmsk);
            assert(0);
         }
         itst++;
      }
   }
   ATL_tFreeBV(gbv);
   ATL_FreeBV(lbv);
}

int main(int nargs, char **args)
{
   ATL_BV_t *bv, *bv0, vv;
   unsigned long i;
   unsigned long k;
   int P=ATL_NTHREADS;

   if (nargs > 1)
      P = atoi(args[1]);
/*
 * First, we will sanity-check basic routs
 */
   bv = ATL_tNewBV(4, 7);
   assert(bv);          /* did it return a non-NULL ptr? */
   assert(ATL_tGetTotBitsBV(bv) == 4);  /* 4-bit length? */
   ATL_tFreeBV(bv);
/*
 * Note this code section tests ATL_IsBitSetBV & ATL_SetBitBV
 */
   bv = ATL_tNewBV(33, 4);
   setRange(bv, 0, 0xAC73B2E1);
                            /* 0b1010 1100 0111 0011 1011 0010 1110 0001 */
   assert(ATL_tIsBitSetBV(bv, 0) == 1); /* 0001 */
   assert(ATL_tIsBitSetBV(bv, 1) == 0);
   assert(ATL_tIsBitSetBV(bv, 2) == 0);
   assert(ATL_tIsBitSetBV(bv, 3) == 0);

   assert(ATL_tIsBitSetBV(bv, 4) == 0); /* 1110 */
   assert(ATL_tIsBitSetBV(bv, 5) == 1);
   assert(ATL_tIsBitSetBV(bv, 6) == 1);
   assert(ATL_tIsBitSetBV(bv, 7) == 1);

   assert(ATL_tIsBitSetBV(bv, 8) == 0); /* 0010 */
   assert(ATL_tIsBitSetBV(bv, 9) == 1);
   assert(ATL_tIsBitSetBV(bv,10) == 0);
   assert(ATL_tIsBitSetBV(bv,11) == 0);

   assert(ATL_tIsBitSetBV(bv,12) == 1); /* 1011 */
   assert(ATL_tIsBitSetBV(bv,13) == 1);
   assert(ATL_tIsBitSetBV(bv,14) == 0);
   assert(ATL_tIsBitSetBV(bv,15) == 1);

   assert(ATL_tIsBitSetBV(bv,16) == 1); /* 0011 */
   assert(ATL_tIsBitSetBV(bv,17) == 1);
   assert(ATL_tIsBitSetBV(bv,18) == 0);
   assert(ATL_tIsBitSetBV(bv,19) == 0);

   assert(ATL_tIsBitSetBV(bv,20) == 1); /* 0111 */
   assert(ATL_tIsBitSetBV(bv,21) == 1);
   assert(ATL_tIsBitSetBV(bv,22) == 1);
   assert(ATL_tIsBitSetBV(bv,23) == 0);

   assert(ATL_tIsBitSetBV(bv,24) == 0); /* C=1100 */
   assert(ATL_tIsBitSetBV(bv,25) == 0);
   assert(ATL_tIsBitSetBV(bv,26) == 1);
   assert(ATL_tIsBitSetBV(bv,27) == 1);

   assert(ATL_tIsBitSetBV(bv,28) == 0); /* A=1010 */
   assert(ATL_tIsBitSetBV(bv,29) == 1);
   assert(ATL_tIsBitSetBV(bv,30) == 0);
   assert(ATL_tIsBitSetBV(bv,31) == 1);

   assert(ATL_tIsBitSetBV(bv, 32) == 0);
   assert(ATL_tSetBitBV(bv, 0) == 1); /* 0001 */
   assert(ATL_tSetBitBV(bv, 1) == 0);
   assert(ATL_tSetBitBV(bv, 2) == 0);
   assert(ATL_tSetBitBV(bv, 3) == 0);

   assert(ATL_tSetBitBV(bv, 4) == 0); /* 1110 */
   assert(ATL_tSetBitBV(bv, 5) == 1);
   assert(ATL_tSetBitBV(bv, 6) == 1);
   assert(ATL_tSetBitBV(bv, 7) == 1);

   assert(ATL_tSetBitBV(bv, 8) == 0); /* 0010 */
   assert(ATL_tSetBitBV(bv, 9) == 1);
   assert(ATL_tSetBitBV(bv,10) == 0);
   assert(ATL_tSetBitBV(bv,11) == 0);

   assert(ATL_tSetBitBV(bv,12) == 1); /* 1011 */
   assert(ATL_tSetBitBV(bv,13) == 1);
   assert(ATL_tSetBitBV(bv,14) == 0);
   assert(ATL_tSetBitBV(bv,15) == 1);

   assert(ATL_tSetBitBV(bv,16) == 1); /* 0011 */
   assert(ATL_tSetBitBV(bv,17) == 1);
   assert(ATL_tSetBitBV(bv,18) == 0);
   assert(ATL_tSetBitBV(bv,19) == 0);

   assert(ATL_tSetBitBV(bv,20) == 1); /* 0111 */
   assert(ATL_tSetBitBV(bv,21) == 1);
   assert(ATL_tSetBitBV(bv,22) == 1);
   assert(ATL_tSetBitBV(bv,23) == 0);

   assert(ATL_tSetBitBV(bv,24) == 0); /* C=1100 */
   assert(ATL_tSetBitBV(bv,25) == 0);
   assert(ATL_tSetBitBV(bv,26) == 1);
   assert(ATL_tSetBitBV(bv,27) == 1);

   assert(ATL_tSetBitBV(bv,28) == 0); /* A=1010 */
   assert(ATL_tSetBitBV(bv,29) == 1);
   assert(ATL_tSetBitBV(bv,30) == 0);
   assert(ATL_tSetBitBV(bv,31) == 1);

   assert(ATL_tSetBitBV(bv, 32) == 0);
   ATL_tFreeBV(bv);
   testRange(3, lrand48(), 7);
   testRange(11, lrand48(), 77);
   testSetUnset(P*7,  P);
   testSetUnset(P*7+1,  P);
   if (P > 1)
      testSetUnset(P*7+1,  (P>>1));
   testGlb2loc(777, P);
   testSetRange(777, P);
   testSetRange(111, P);
/*
 * If no assertion failed, print success and return 0 for no error!
 */
   printf("\nSUCCESS!\n\n");
   return(0);
}
