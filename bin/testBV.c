/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "atlas_bitvec.h"

int RangeSet(ATL_BV_t *bv, unsigned int b0, unsigned int b1)
{
   unsigned int i;
   if (b0 > b1 || b1 >= ATL_GetTotBitsBV(bv) || !bv)
      return(0);
   for (i=b0; i <= b1; i++)
      if (!ATL_IsBitSetBV(bv, i))
         return(0);
   return(1);
}
int TestRandRange
(
   unsigned int nbits,
   unsigned int ntest
)
{
   const unsigned int N=(nbits+bpiBV-1)>>shBV;
   unsigned int nerr=0, t;
   ATL_BV_t *bv;
   printf("\nRunning %lu random tests on %u-bit IsBitRangeSetBV:\n",
          ((unsigned long)ntest)*nbits*nbits, nbits);
   srand48(nbits+(ntest<<7));  /* make test repeatable for each N,ntest pair */
   bv = ATL_NewBV(nbits);
/*
 * Check error cases
 */
   assert(ATL_IsBitRangeSetBV(NULL, 0, 0) == 0);
   assert(ATL_IsBitRangeSetBV(bv, 1, 0) == 0);
   assert(ATL_IsBitRangeSetBV(bv, nbits, nbits) == 0);
   assert(ATL_IsBitRangeSetBV(bv, nbits<<1, nbits<<2) == 0);

   for (t=0; t < ntest; t++)
   {
      unsigned int k, b0, b1;
/*
 *    Randomly choose bit patterns for whole vector
 */
      for (k=0; k < N; k++)
         ATL_SetEltBV(bv, k, lrand48());
/*
 *    Now try all possible legal pairs
 */
      for (b0=0; b0 < nbits; b0++)
      {
         for (b1=b0; b1 < nbits; b1++)
         {
            int good, chk;
            chk = ATL_IsBitRangeSetBV(bv, b0, b1);
            good = RangeSet(bv, b0, b1);
            if (chk != good)
            {
               printf("ERROR: nbits=%d, [%d,%d], expected=%d, got=%d\n",
                      nbits, b0, b1, good, chk);
               nerr++;
            }
         }
      }
   }
   ATL_FreeBV(bv);
   if (nerr)
      fprintf(stderr, "NFAILURES=%d\n\n", nerr);
   else
      printf("ALL CASES PASSED\n\n");
   return(nerr);
}
int TestRandPos
(
   unsigned int nbits,   /* total number of bits to have in array */
   unsigned int ntest,   /* # of tests to run */
   int set               /* test FindFirst: 0: Unset : else : Set */
)
{
   ATL_BV_t *bv;
   int i, b;
   unsigned int nelt, t, b0, b1, i0, i1, p0, p1, nerr=0;
   long (*findFirst)(ATL_BV_t *bv, unsigned long bs);
   int (*chgBit)(ATL_BV_t *bv, unsigned long pos);

   printf("\nRunning %u random tests on %u-bit findFirst%sBV:\n", ntest, nbits,
          set ? "Set":"Unset");
   srand(nbits+(ntest<<7));  /* make test repeatable for each N,ntest pair */
   bv = ATL_NewBV(nbits);
   if (set)
   {
      findFirst = ATL_FindFirstSetBitBV;
      chgBit = ATL_SetBitBV;
   }
   else
   {
      findFirst = ATL_FindFirstUnsetBitBV;
      chgBit = ATL_UnsetBitBV;
      ATL_SetAllBitsBV(bv);
   }
   nelt = (nbits + bpiBV-1)/bpiBV;

   assert(findFirst(NULL, 0) == -1);      /* test NULL handling */
   assert(findFirst(bv, nbits) == -1);    /* test out-of-range */
   assert(findFirst(bv, nbits+1) == -1);  /* test out-of-range */
   assert(findFirst(bv, 2*nbits) == -1);  /* test out-of-range */
   assert(findFirst(bv, 0) == -1);        /* test no unset bits found */
   for (t=0; t < ntest; t++)
   {
/*
 *    Randomly generate two bits to unset, put smallest in b0
 */
      b0 = rand() % nbits;
      b1 = rand() % nbits;
      if (b0 > b1)
      {
         unsigned int b=b0;
         b0 = b1;
         b1 = b;
      }
/*
 *    From global bit number, compute integer index (i0/1) and pos (p0/1)
 *    zero out two bits so we can test that least sig is returned.
 */
      i0 = b0 / bpiBV;
      p0 = b0 - i0*bpiBV;
      i1 = b1 / bpiBV;
      p1 = b1 - i1*bpiBV;
      chgBit(bv, b0);
      chgBit(bv, b1);
/*
 *    Test we find first unset bit
 */
      b = findFirst(bv, 0);
      if (b != b0)
      {
         fprintf(stderr, "   FAILED %d: expected=%u, got=%d\n",
                 __LINE__, b0, b);
         if (b == b1)
            fprintf(stderr, "     Are you returning the MOST sig bit?\n");
         nerr++;
      }
/*
 *    Test that skipping unused elements doesn't change answer
 */
      b = findFirst(bv, i0*bpiBV);
      if (b != b0)
      {
      b = findFirst(bv, i0*bpiBV);
         fprintf(stderr, "   FAILED %d: expected=%u, got=%d\n",
                 __LINE__, b0, b);
         fprintf(stderr, "     Is your skp working correctly?\n");
         nerr++;
      }
/*
 *    Test that skipping i0 means we find p1
 */
      if (b1 > b0)
      {
         b = findFirst(bv, b0+1);
         if (b != b1)
         {
         b = findFirst(bv, b0+1);
            fprintf(stderr, "   FAILED %d: expected=%u, got=%d\n",
                    __LINE__, b1, b);
            fprintf(stderr, "     Is your skp working correctly PART 2?\n");
            nerr++;
         }
      }
      if (set)
         ATL_UnsetAllBitsBV(bv);  /* revert to all unset */
      else
         ATL_SetAllBitsBV(bv);    /* revert to all set */
   }
/*
 * Now see if we can find a value in very last element, which will be partial
 * if N%bpiBV != 0
 */
   b0 = nbits-1;
   i0 = nelt-1;
   p0 = b0 - i0*bpiBV;
   if (set)
      ATL_SetBitBV(bv, b0);
   else
      ATL_UnsetBitBV(bv, b0);
   b = findFirst(bv, 0);
   if (b != b0)
   {
      fprintf(stderr, "   FAILED %d: expected=%u, got=%d\n", __LINE__, b0, b);
      if (nbits%bpiBV)
         fprintf(stderr, "     Are you handling partial last int?\n");
      nerr++;
   }
/*
 * Try partial last int with random location
 */
   i1 = nbits%bpiBV;
   if (i1 > 1)
   {
      p0 = rand()%(i1-1);
      b0 = i0*bpiBV + p0;
      if (set)
         ATL_SetBitBV(bv, b0);
      else
         ATL_UnsetBitBV(bv, b0);
      b = findFirst(bv, 0);
      if (b != b0)
      {
         fprintf(stderr, "   FAILED %d: expected=%u, got=%d\n",
                 __LINE__, b0, b);
         fprintf(stderr, "     Are you handling partial last int PART 2?\n");
         nerr++;
      }
   }
/*
 * Put unset bits past end of bv and make sure they are not returned!
 */
   if (i1)
   {
      ATL_BV_t msk=0;
      for (i=0; i < i1; i++)
         msk |= (1L<<i);
      if (set)
         msk = ~msk;
      ATL_SetEltBV(bv, i0, msk);
      b = findFirst(bv, 0);
      if (b != -1)
      {
         fprintf(stderr, "   FAILED %d: expected=%d, got=%d\n",
                 __LINE__, -1, b);
         fprintf(stderr, "     Are you handling partial last int as full?\n");
         nerr++;
      }
   }
/*
 * Need to make test cases for partial first element
 */
   ATL_FreeBV(bv);
   if (nerr)
      fprintf(stderr, "NFAILURES=%d\n\n", nerr);
   else
      printf("ALL CASES PASSED\n\n");
   return(nerr);
}
int TestAllPos  /* returns # of errors found */
(
   int nwords,  /* will alloc bv of size nwords*bpiBV */
   int word,    /* will exhaustively test this word */
   int set      /* test FindFirst: 0: Unset : else : Set */
)
{
   ATL_BV_t *bv, mask;
   int k, i;
   unsigned int nerr=0, skip = word<<shBV;
   long (*findFirst)(ATL_BV_t *bv, unsigned long bs);

   assert(nwords > word);
   bv = ATL_NewBV(nwords<<shBV);
   if (!set)
      ATL_SetAllBitsBV(bv);
   findFirst = (set) ? ATL_FindFirstSetBitBV : ATL_FindFirstUnsetBitBV;
   mask = allsetBV;
   if (set)
      mask = ~mask;
   ATL_SetEltBV(bv, word, mask);
   if (findFirst(bv, skip) != -1)
   {
      fprintf(stderr, "All bits set does not return -1!\n");
      nerr++;
   }
/*
 * Test we can get correct answer with only one unset bit at each loc
 */
   for (i=0; i < bpiBV-1; i++)
   {
      mask = allsetBV & ~(1L<<i);
      if (set)
         mask = ~mask;
      ATL_SetEltBV(bv, word, mask);
      k = findFirst(bv, skip);
      if (k != i+skip)
      {
      k = findFirst(bv, skip);
         fprintf(stderr, "EXPECTED=%d, GOT=%d, skip=%d\n", i+skip, k, skip);
         nerr++;
      }
   }
   if (nerr)            /* don't test ties (confusing) */
      return(nerr);     /* until we can get right ans wt only 1 bit set */
/*
 * Make sure that when we have two locations unset, we always pick least sig
 */
   mask = 0;
   if (set)
      mask = ~mask;
   ATL_SetEltBV(bv, word, mask);
   k = findFirst(bv, skip);
   if (k != skip)
   {
      fprintf(stderr, "TIE: EXPECTED=%d, GOT=%d\n", skip, k);
      nerr++;
   }
/*
 * Test all combinations of two bits set, ensure smallest always picked
 */
   for (i=0; i < bpiBV-1; i++)
   {
      int j;
      mask = allsetBV & ~(1L<<i);
      ATL_SetEltBV(bv, word, mask);
      for (j=i; j < bpiBV; j++)
      {
         ATL_BV_t msk;
         msk = mask & ~(1L<<j);
         if (set)
            msk = ~msk;
         ATL_SetEltBV(bv, word, msk);
         k = findFirst(bv, skip);
         if (k != i+skip)
         {
            fprintf(stderr, "TIE: EXPECTED=%d, GOT=%d\n", i+skip, k);
            nerr++;
         }
      }
   }
   ATL_FreeBV(bv);
/*
 * Test that the AllBits funcs don't mess with remainder values
 */
   for (i=0; i < bpiBV; i++)
   {
      bv = ATL_NewBV((2<<shBV)+i);
      assert(ATL_GetEltBV(bv, 0) == 0);
      assert(ATL_GetEltBV(bv, 1) == 0);
      if (i)
         assert(ATL_GetEltBV(bv, 2) == 0);
      ATL_SetAllBitsBV(bv);
      assert(ATL_GetEltBV(bv, 0) == allsetBV);
      assert(ATL_GetEltBV(bv, 1) == allsetBV);
      if (i)  /* should have i bits set, rest unset */
      {
         ATL_BV_t v;
         v = ATL_GetEltBV(bv, 2);
         for (k=0; k < i; k++)
            assert(v&(1L<<k));
         for (k=i; k < bpiBV; k++)
            assert(!(v&(1L<<k)));
      }
/*
 *    Set all bits (even past end)
 */
      ATL_SetEltBV(bv, 0, allsetBV);
      ATL_SetEltBV(bv, 1, allsetBV);
      if (i)
         ATL_SetEltBV(bv, 2, allsetBV);
      ATL_UnsetAllBitsBV(bv);
      assert(ATL_GetEltBV(bv, 0) == 0);
      assert(ATL_GetEltBV(bv, 1) == 0);
      if (i)  /* should have i bits unset, rest set */
      {
         ATL_BV_t v;
         v = ATL_GetEltBV(bv, 2);
         for (k=0; k < i; k++)
            assert(!(v&(1L<<k)));
         for (k=i; k < bpiBV; k++)
            assert(v&(1L<<k));
      }
/*
 *    Set all bits (even past end)
 */
      ATL_SetEltBV(bv, 0, allsetBV);
      ATL_SetEltBV(bv, 1, allsetBV);
      if (i)
         ATL_SetEltBV(bv, 2, allsetBV);
      ATL_ReverseAllBitsBV(bv);
      assert(ATL_GetEltBV(bv, 0) == 0);
      assert(ATL_GetEltBV(bv, 1) == 0);
      if (i)  /* should have i bits unset, rest set */
      {
         ATL_BV_t v;
         v = ATL_GetEltBV(bv, 2);
         for (k=0; k < i; k++)
            assert(!(v&(1L<<k)));
         for (k=i; k < bpiBV; k++)
            assert(v&(1L<<k));
      }
      ATL_ReverseAllBitsBV(bv);  /* all bits set again */
      assert(ATL_GetEltBV(bv, 0) == allsetBV);
      assert(ATL_GetEltBV(bv, 1) == allsetBV);
      if (i)
         assert(ATL_GetEltBV(bv, 2) == allsetBV);
      ATL_FreeBV(bv);
   }
   if (nerr)
      fprintf(stderr, "FindFirst%sBit: NFAILURES=%d\n",
              set ? "Set":"Unset", nerr);
   else
      printf("FindFirst%sBit: ALL CASES PASSED\n", set ? "Set":"Unset");
   return(nerr);
}

int testIncorp(long seed)
{
   int sl;
   int nerr=0;
   ATL_BV_t *src, *dst;
/*
 * Sanity test with 0 bounds (1 bounds tested in main loop)
 */
   srand48(seed);
   src = ATL_NewBV(3);
   dst = ATL_NewBV(5);         /*  00000 (right least sig) */
   ATL_ReverseAllBitsBV(src);  /*    111 */
   ATL_IncorpBV(dst, src, 1);  /* 001110 */
   nerr = (ATL_IsBitSetBV(dst, 0) != 0);
   nerr |= (ATL_IsBitSetBV(dst, 1) != 1);
   nerr |= (ATL_IsBitSetBV(dst, 2) != 1);
   nerr |= (ATL_IsBitSetBV(dst, 3) != 1);
   nerr |= (ATL_IsBitSetBV(dst, 4) != 0);
   nerr |= (ATL_IsBitSetBV(dst, 5) != 0);
   if (nerr)
      fprintf(stderr, "INCORP MISMATCH: exp=%x, got=%x\n", 0xE, dst[1]);
   assert(!nerr);
   ATL_FreeBV(src);
   ATL_FreeBV(dst);
   for (sl=0; sl < 4*bpiBV; sl++)  /* source length */
   {
      const int nselt=(sl+bpiBV-1)>>shBV;
      ATL_BV_t pos;
      src = ATL_NewBV(sl);
      for (pos=0; pos < nselt; pos++)
         ATL_SetEltBV(src, pos, lrand48()); /* set src to random bit pattern */
      for (pos=0; pos <= 2*bpiBV; pos++)
      {
         ATL_BV_t i;
         const int dl = sl+pos+7;
         dst = ATL_NewBV(dl);         /* all 0s */
         ATL_ReverseAllBitsBV(dst);  /* all 1s */
         ATL_IncorpBV(dst, src, pos);
         for (i=0; i < sl; i++)
         {
            if (ATL_IsBitSetBV(src, i) != ATL_IsBitSetBV(dst, i+pos))
            {
               fprintf(stderr, "INCORP MISMATCH SLEN=%ld, POS=%ld: bit=%ld\n",
                       sl, pos, i);
               nerr++;
               ATL_assert(!nerr);
            }
         }
         for (i=0; i < pos; i++)
         {
            if (!ATL_IsBitSetBV(dst, i))
            {
               fprintf(stderr, "POS=%ld, SLEN=%ld: LOWER BOUND %ld ZEROED!\n",
                       pos, sl, i);
               nerr++;
               ATL_assert(!nerr);
            }
         }
         for (i=pos+sl; i < dl; i++)
         {
            if (!ATL_IsBitSetBV(dst, i))
            {
               fprintf(stderr, "POS=%ld, SLEN=%ld: UPPER BOUND %ld ZEROED!\n",
                       pos, sl, i);
               nerr++;
               ATL_assert(!nerr);
            }
         }
         ATL_FreeBV(dst);
      }
      ATL_FreeBV(src);
   }
   return(nerr);
}
int main(int nargs, char **args)
{
   ATL_BV_t *bv, *bv0, vv;
   unsigned int i;
   unsigned int k;
/*
 * First, we will sanity-check basic routs
 */
   bv = ATL_NewBV(4);
   assert(bv);          /* did it return a non-NULL ptr? */
   assert(ATL_GetTotBitsBV(bv) == 4);  /* 4-bit length? */
/*
 * Basic set/GetElt test
 */
   assert(ATL_SetEltBV(NULL, 0, 0) == 0);
   assert(ATL_SetEltBV(bv, 1, 0) == 0);
   assert(ATL_SetEltBV(bv, 2, 0) == 0);
   k = bpiBV+4;
   bv = ATL_ExpandBV(bv, k);    /* expand bv to require 2 ints to store */
   assert(bv);
   assert(ATL_GetTotBitsBV(bv) == k);  /* k-bit length? */
/*
 * -1 is just a bitpattern, that is all bits set, is why this should work
 * even though argument is unsigned!
 */
   ATL_SetEltBV(bv, 0, ~((ATL_BV_t)(0)));
   assert(ATL_GetEltBV(bv, 0) == allsetBV);
   assert(ATL_GetEltBV(bv, 1) == 0);
   bv0 = bv;
/*
 * Check that "expanding" vec to smaller size has no affect
 */
   bv = ATL_ExpandBV(bv, 4);
   assert(bv == bv0);
   assert(ATL_GetTotBitsBV(bv) == k);

   ATL_SetAllBitsBV(NULL);
   ATL_SetAllBitsBV(bv);
   assert(ATL_GetEltBV(bv, 0) == allsetBV);
   assert(ATL_GetEltBV(bv, 1) == 0xF);

   ATL_UnsetAllBitsBV(NULL);
   ATL_SetEltBV(bv, 1, -1);              /* partial vec to all bits set */
   ATL_UnsetAllBitsBV(bv);
   assert(ATL_GetEltBV(bv, 0) == 0);     /* all bits unset for full int is 0 */
   vv = ATL_GetEltBV(bv, 1);
   if (vv != ~0xFL)                      /* only low 4 bits should change */
   {
      printf("expected=%lx, got=%lx!\n", vv, ~0xFL);
      assert(vv == ~(0xFL)); /* only low 4 bits changed */
   }
   ATL_FreeBV(bv);                       /* try to free bitvec */
/*
 * Let's see if your ATL_GetTotBitsBV works!
 */
   bv = ATL_NewBV(5);
   i = ATL_GetTotBitsBV(bv);
   if (i != 5)
   {
      fprintf(stderr, "ATL_GetTotBitsBV: expected=5, got=%d!\n", i);
      assert(i == 5);
   }
   bv = ATL_ExpandBV(bv, 31);
   i = ATL_GetTotBitsBV(bv);
   if (i != 31)
   {
      fprintf(stderr, "ATL_GetTotBitsBV: expected=31, got=%d!\n", i);
      assert(i == 31);
   }
   bv = ATL_ExpandBV(bv, bpiBV+1);
   assert(ATL_GetTotBitsBV(bv) == bpiBV+1);
   ATL_FreeBV(bv);

   bv = ATL_NewBV(8);
   assert(ATL_GetEltBV(bv, 0) == 0);
   ATL_SetAllBitsBV(bv);
   assert(ATL_GetEltBV(bv, 0) == 0xFF);
   ATL_FreeBV(bv);
/*
 * Note this code section tests ATL_IsBitSetBV & ATL_SetBitBV
 */
   bv = ATL_NewBV(33);
   assert(ATL_SetEltBV(bv, 0, 0xAC73B2E1) == 0);
                            /* 0b1010 1100 0111 0011 1011 0010 1110 0001 */
   assert(ATL_IsBitSetBV(bv, 0) == 1); /* 0001 */
   assert(ATL_IsBitSetBV(bv, 1) == 0);
   assert(ATL_IsBitSetBV(bv, 2) == 0);
   assert(ATL_IsBitSetBV(bv, 3) == 0);

   assert(ATL_IsBitSetBV(bv, 4) == 0); /* 1110 */
   assert(ATL_IsBitSetBV(bv, 5) == 1);
   assert(ATL_IsBitSetBV(bv, 6) == 1);
   assert(ATL_IsBitSetBV(bv, 7) == 1);

   assert(ATL_IsBitSetBV(bv, 8) == 0); /* 0010 */
   assert(ATL_IsBitSetBV(bv, 9) == 1);
   assert(ATL_IsBitSetBV(bv,10) == 0);
   assert(ATL_IsBitSetBV(bv,11) == 0);

   assert(ATL_IsBitSetBV(bv,12) == 1); /* 1011 */
   assert(ATL_IsBitSetBV(bv,13) == 1);
   assert(ATL_IsBitSetBV(bv,14) == 0);
   assert(ATL_IsBitSetBV(bv,15) == 1);

   assert(ATL_IsBitSetBV(bv,16) == 1); /* 0011 */
   assert(ATL_IsBitSetBV(bv,17) == 1);
   assert(ATL_IsBitSetBV(bv,18) == 0);
   assert(ATL_IsBitSetBV(bv,19) == 0);

   assert(ATL_IsBitSetBV(bv,20) == 1); /* 0111 */
   assert(ATL_IsBitSetBV(bv,21) == 1);
   assert(ATL_IsBitSetBV(bv,22) == 1);
   assert(ATL_IsBitSetBV(bv,23) == 0);

   assert(ATL_IsBitSetBV(bv,24) == 0); /* C=1100 */
   assert(ATL_IsBitSetBV(bv,25) == 0);
   assert(ATL_IsBitSetBV(bv,26) == 1);
   assert(ATL_IsBitSetBV(bv,27) == 1);

   assert(ATL_IsBitSetBV(bv,28) == 0); /* A=1010 */
   assert(ATL_IsBitSetBV(bv,29) == 1);
   assert(ATL_IsBitSetBV(bv,30) == 0);
   assert(ATL_IsBitSetBV(bv,31) == 1);

   assert(ATL_IsBitSetBV(bv, 32) == 0);
   assert(ATL_SetBitBV(bv, 0) == 1); /* 0001 */
   assert(ATL_SetBitBV(bv, 1) == 0);
   assert(ATL_SetBitBV(bv, 2) == 0);
   assert(ATL_SetBitBV(bv, 3) == 0);

   assert(ATL_SetBitBV(bv, 4) == 0); /* 1110 */
   assert(ATL_SetBitBV(bv, 5) == 1);
   assert(ATL_SetBitBV(bv, 6) == 1);
   assert(ATL_SetBitBV(bv, 7) == 1);

   assert(ATL_SetBitBV(bv, 8) == 0); /* 0010 */
   assert(ATL_SetBitBV(bv, 9) == 1);
   assert(ATL_SetBitBV(bv,10) == 0);
   assert(ATL_SetBitBV(bv,11) == 0);

   assert(ATL_SetBitBV(bv,12) == 1); /* 1011 */
   assert(ATL_SetBitBV(bv,13) == 1);
   assert(ATL_SetBitBV(bv,14) == 0);
   assert(ATL_SetBitBV(bv,15) == 1);

   assert(ATL_SetBitBV(bv,16) == 1); /* 0011 */
   assert(ATL_SetBitBV(bv,17) == 1);
   assert(ATL_SetBitBV(bv,18) == 0);
   assert(ATL_SetBitBV(bv,19) == 0);

   assert(ATL_SetBitBV(bv,20) == 1); /* 0111 */
   assert(ATL_SetBitBV(bv,21) == 1);
   assert(ATL_SetBitBV(bv,22) == 1);
   assert(ATL_SetBitBV(bv,23) == 0);

   assert(ATL_SetBitBV(bv,24) == 0); /* C=1100 */
   assert(ATL_SetBitBV(bv,25) == 0);
   assert(ATL_SetBitBV(bv,26) == 1);
   assert(ATL_SetBitBV(bv,27) == 1);

   assert(ATL_SetBitBV(bv,28) == 0); /* A=1010 */
   assert(ATL_SetBitBV(bv,29) == 1);
   assert(ATL_SetBitBV(bv,30) == 0);
   assert(ATL_SetBitBV(bv,31) == 1);

   assert(ATL_SetBitBV(bv, 32) == 0);
   #if bpiBV == 64
      assert(ATL_GetEltBV(bv, 0) == 0x1FFFFFFFF);
   #else
      assert(ATL_GetEltBV(bv, 0) == 0xFFFFFFFF);
   #endif
   #if bpiBV > 32
   bv = ATL_ExpandBV(bv, bpiBV+1);
   ATL_SetEltBV(bv, 0, 0xF854B36D00000000L);
   assert(ATL_GetEltBV(bv, 0) ==  0xF854B36D00000000L);
   for (i=0; i < 32; i++)
   {
      if (ATL_IsBitSetBV(bv, i+32))
         assert(ATL_SetBitBV(bv, i) == 0);
   }
   assert(ATL_GetEltBV(bv, 0) ==  0xF854B36DF854B36DL);
   assert((ATL_GetEltBV(bv, 0)>>32) == 0xF854B36D);
   #endif
   ATL_FreeBV(bv);
/*
 * Check error exits of accessor functions
 */
   bv = ATL_NewBV(65);
   assert(ATL_UnsetBitBV(NULL, 100000) == -1);
   assert(ATL_UnsetBitBV(bv, 100000) == -1);
   assert(ATL_UnsetBitBV(bv, 65) == -1);
   assert(ATL_SetBitBV(NULL, 100000) == -1);
   assert(ATL_SetBitBV(bv, 100000) == -1);
   assert(ATL_SetBitBV(bv, 65) == -1);
   assert(ATL_GetEltBV(NULL, 100000) == 0);
   assert(ATL_GetEltBV(bv, 100000) == 0);
   assert(ATL_GetEltBV(bv, 65) == 0);
   assert(ATL_IsBitSetBV(NULL, 100000) == 0);
   assert(ATL_IsBitSetBV(bv, 100000) == 0);
   assert(ATL_IsBitSetBV(bv, 65) == 0);
   ATL_FreeBV(bv);
   testIncorp(0xfc135);
   assert(TestAllPos(6,  5, 1) == 0);
   assert(TestAllPos(6,  3, 1) == 0);
   assert(TestAllPos(6,  0, 1) == 0);
   assert(TestAllPos(6,  5, 0) == 0);
   assert(TestAllPos(6,  3, 0) == 0);
   assert(TestAllPos(6,  0, 0) == 0);

   assert(TestRandPos(7777, 8000, 1) == 0);
   assert(TestRandPos(6400, 8000, 1) == 0);
   assert(TestRandPos(5555, 8000, 0) == 0);
   assert(TestRandPos(3200, 8000, 0) == 0);
   assert(TestRandRange(7777, 2) == 0);
/*
 * If no assertion failed, print success and return 0 for no error!
 */
   printf("\nSUCCESS!\n\n");
   return(0);
}
