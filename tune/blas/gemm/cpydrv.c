#include "atlas_mmgen.h"
#define uint unsigned int

ATL_cpnode_t **SplitList(char pre, ATL_cpnode_t *cb, int *NLIST, uint **FLAGS)
/*
 * Splits cb according to GetCopyMatches, returns number of non-NULL lists
 * in *NLIST, and the non-NULL flag values in FLGS.
 * RETURNS: *NLIST-length array of lists
 */
{
   int ik; /* loop over kernels, which can be incoherent */
   int nlist=0, nalloc;
   int *flags;
   ATL_cpnode_t **lists;

   nalloc = 8;
   lists = malloc(nalloc*sizeof(ATL_cpnode_t *));
   flags = malloc(nalloc*sizeof(int));
   assert(lists && flags);

   CopyApplyBlasRules(cb);
   for (ik=CPF_SYRK-1; ik <= CPF_TRSM; ik++)
   {
      int ic; /* CBLK unset or set */
      for (ic=0; ic < 2; ic++)
      {
         const int B0=(ic)?CPF_BE1:0, BN=(ic)?CPF_BE0:1;
         int id; /* TOBLK unset or set */
         for (id=0; id < 2; id++)
         {
            int ia;
            for (ia=CPF_AL1; ia <= CPF_ALX; ia++)
            {
               int ib;

               for (ib=B0; ib <= BN; ib++)
               {
                  const int HN = (ic || pre == 'd' || pre == 's') ? 0 : 1;
                  int ih;
                  for (ih=0; ih <= HN; ih++)
                  {
                     ATL_cpnode_t *bp;
                     int flag;
                     flag = (id<<CPF_TOBLK)|(ic<<CPF_CBLK)|(1<<ia);
                     flag |= (ik >= CPF_SYRK) ? (1<<ik) : 0;
                     if (ic)
                        flag |= 1<<ib;
                     else
                        flag |= (ib<<CPF_TRANS)|(ih<<CPF_CONJ);
                     bp = GetCopyMatches(cb, flag);
                     if (bp)
                     {
                        if (nlist == nalloc)
                        {
                           nalloc = nlist + 4;
                           lists = realloc(lists,nalloc*sizeof(ATL_cpnode_t *));
                           flags = realloc(flags,nalloc*sizeof(int));
                        }
                        lists[nlist] = bp;
                        flags[nlist++] = flag;
                     }
                  }
               }
            }
         }
      }
   }
   *NLIST = nlist;
   *FLAGS = flags;
   return(lists);
}

unsigned int CountRegion(int N, ATL_cpnode_t **cbs)
{
   const uint mask=(CPF_ALLKERN|(1<<CPF_CBLK)|(1<<CPF_TOBLK)),
              flag=(cbs[0]->flag)&mask;
   uint i;

   for (i=1; i < N; i++)
      if (flag != ((cbs[i]->flag)&mask))
         return(i);
   return(N);
}

void DoShareBV(const uint N,  ATL_cpnode_t **cbs)
/*
 * cbs is N-length array of cpnodes that must cohere (use same idx for same cp)
 * This routine creates a bitmask in ->bv.  Having bit 0 set means that this
 * cpy func appears in the first list, and so on.
 */
{
   int i;

   for (i=0; i < N; i++)
   {
      ATL_cpnode_t *p0;
      for (p0=cbs[i]; p0; p0 = p0->next)
      {
         int j;

         p0->bv |= (1<<i);
         for (j=i+1; j < N; j++)
         {
            ATL_cpnode_t *p1;
            p1 = FindCopy_cohere(cbs[j], p0);
            if (p1)
            {
               p0->bv |= (1<<j);
               p1->bv |= (1<<i);
            }
         }
      }
   }
}

uint *CountDictShared(const uint N, ATL_cpnode_t **cbs)
/*
 * Counts how many dictating entries their are in each of the N cp queus of cbs
 * there.  A dictating entry is any one where its not the two extreme cases:
 * shared by everyone (so any file can dictate w/o problem) or shared by
 * nobody except this queue (private, so we can use any index we like).
 */
{
   uint *cnts, i;
   const uint mall=(1<<N)-1;
   cnts = malloc(N*sizeof(uint));
   assert(cnts);
   for (i=0; i < N; i++)
   {
      const uint mprv=(1<<i);
      uint bv, cnt=0;
      ATL_cpnode_t *p;
      for (p=cbs[i]; p; p = p->next)
      {
         const uint flg=p->bv;
         if (flg != mall && flg != mprv)
            cnt++;
      }
      cnts[i] = cnt;
   }
   return(cnts);
}


int NumSetBits(const uint N, uint flg)
/*
 * RETURNS: how many bits set in 1st N bits of flg
 */
{
   uint nset=0, i;
   for (i=0; i < N; i++)
      if (flg&(1<<i))
         nset++;
   return(nset);
}

ATL_cpnode_t *findMaxSetBits(const uint N, ATL_cpnode_t *cb)
/*
 * RETURNS: node in cb that has the most bits set in bv
 */
{
   if (cb)
   {
      ATL_cpnode_t *cp, *mp=cb;;

      if (cb->next)
      {
         uint nset;
         nset = NumSetBits(N, cb->bv);
         for (cp=cb->next; cp; cp = cp->next)
         {
            uint ns = NumSetBits(N, cp->bv);
            if (ns > nset)
            {
               mp = cp;
               nset = ns;
            }
         }
      }
      return(mp);
   }
   return(NULL);
}

void moveCopiesInMask
(
   const uint N,         /* number of files being created */
   uint *cnts,           /* count of nodes in cbs */
   ATL_cpnode_t **cbs,   /* unsorted copies that need to cohere */
   ATL_cpnode_t **srtb,  /* queue of sorted-by-idx that already cohere */
   ATL_cpnode_t *mp,     /* cohere pattern to move */
   const int msk         /* bitvec indicating kern to move cbs->srtb */
)
{
   uint i;
   for (i=0; i < N; i++)
   {
      if (msk & (1<<i))
      {
         ATL_cpnode_t *cp;
         cp = FindCopy_cohere(cbs[i], mp);
         assert(cp);
         cbs[i] = RemoveCPNodeFromQ(cbs[i], cp);
         cnts[i]--;
         srtb[i] = ATL_JoinCPQs(srtb[i], cp);
      }
   }
}

uint coherePatt2Arr     /* RETURNS: number of unique coherence patterns */
(
   ATL_cpnode_t *cb,    /* list of all copies in file */
   uint *bp             /* uint array long enough to hold all unique bitpatts */
)
{
   ATL_cpnode_t *cp;
   int n=0;
   for (cp=cb; cp; cp = cp->next)
   {
      uint patt = cp->bv;
      uint i;

      for (i=0; i < n && bp[i] != patt; i++); /* see if pattern already there */
      if (i == n)
         bp[n++] = patt;
   }
   return(n);
}

uint maxSpanIdx
(
   uint B,   /* max number of bits in pattern */
   uint N,   /* number of patterns to search in bp */
   uint *bp  /* N-length cnt:patterns */
)
{
   if (N)
   {
      uint vmax = (*bp)>>B, imax=0, i;
      for (i=1; i < N; i++)
      {
         uint v = bp[i] >> N;
         if (v > vmax)
         {
            vmax = v;
            imax = i;
         }
      }
      return(imax);
   }
   return(0);
}

uint sortPattBySpan
(
   uint B,   /* max number of bits in pattern */
   uint N,   /* length of bit pattern array */
   uint *bp  /* 1st entry length, rest bit patterns */
)
/*
 * Simple greatest-to-least selection sort on number of bits set in patterns
 */
{
   uint i;
   for (i=0; i < N-1; i++)
   {
      uint j;
      j = maxSpanIdx(B, N-i, bp+i);
      if (j != i)
      {
         uint t;
         t = bp[i];
         bp[i] = bp[j];
         bp[j] = t;
      }
   }
}

uint addSetCnt2Patt
(
   uint B,   /* max number of bits in pattern */
   uint N,   /* length of bit pattern array */
   uint *bp  /* 1st entry length, rest bit patterns */
)
/*
 * Adds length of bit pattern to pattern encoding above bit N
 */
{
   const uint mall=(1<<B)-1;
   uint i;
   for (i=0; i < N; i++)
   {
      uint pat=(bp[i] & mall);
      if (pat)
         pat |= NumSetBits(B, pat)<<B;
      bp[i] = pat;
   }
}

uint selDrawFile
(
   int N,      /* total number of files */
   uint *cnts, /* N-length array of ints with remaining kern count */
   uint **bps  /* N-length array of arrays of bitpatterns */
)
/*
 * Tries to select the best file to drive the coherence search.  At present,
 * we return the file with the least number of kerns that must be cohered,
 * with idea that once we move this to 0, we can remove the number of files
 * that must be cohere.  When two files have same kern count, we break ties
 * by using the one with the most number of bitpatterns, as that will maximize
 * our flexibility in finding a spanning solution
 */
{
   uint idx=0, minv=0, i;
   for (i=0; i < N; i++)
   {
      const uint cnt = cnts[i];
      if (cnt)
      {
         if (!minv || cnt < minv)
         {
            idx = i;
            minv = cnt;
         }
         else if (cnt == minv) /* tie on count! */
         {
            if (bps[i] > bps[idx])
            {
               idx = i;
               minv = cnt;
            }
         }
      }
   }
   return(idx);
}

uint findMaxSpanR(uint N, uint mcoh, uint msk, uint mneed, uint *cnts,
                  uint **bps, uint mprv, uint *seq)
/*
 * Searches through all files in mneed, and returns the sequence of the mask
 * that spans mneed to the greatest possible extent
 * RETURNS: or of all sequence masks
 */
{
   const uint mall=(1<<N)-1;
   uint i, maxset=0, maxM=0;

   if (!(mneed & (mneed-1)))     /* if only 1 bit set can't fix here,*/
   {                             /* since all patterns have at least 2 bits */
      *seq = 0;                  /* set, so finish sequence */
      return(msk);               /* and return all bits we managed to fix */
   }
   for (i=0; i < N; i++)
   {
      if (cnts[i] && (mneed&(1<<i)))
      {
         uint *bp=bps[i]+1;
         const uint NP = bp[-1];
         uint j, ns;
/*
 *       Look through all bitpatts in file, and see find one that maximizes span
 */
         for (j=0; j < NP; j++)
         {
            const uint m=bp[j]&mall;
            uint NEW=0;

            if (m == mneed)
            {
               *seq = m;
               seq[1] = 0;
               return(msk|m);
            }
            else if ((m|mneed) != mneed) /* if m not a subset of mneed */
               continue;                 /* can't use it */
            else if (!(m&mneed))         /* if m contains no bits we need */
               continue;                 /* can't use it */
            ns = NumSetBits(N, m);
/*
 *          If two sequences have same number of unsatisfied dep, do tiebreak
 *          based on how many holes we can plug using kerns that appear in only
 *          one file!
 */
            if (ns > maxset)
               NEW = 1;
            else if (ns == maxset)
            {
               uint m0, m1;
               m0 = (maxM|mprv)&mneed;
               m1 = (m|mprv)&mneed;
               NEW = (NumSetBits(N, m0) < NumSetBits(N, m1));
            }
            if (NEW)
            {
               maxset = ns;
               *seq = maxM = m;
            }
         }
/*
 *    If I reach here, I was unable to find a pattern that completely covers
 *    mneed, so take the biggest span found, and then see if I can complete
 *    from there.
 */
         mneed ^= maxM;
         msk |= maxM;
         if (maxM)  /* only continue if making progress! */
            msk = findMaxSpanR(N, mcoh, msk, mneed, cnts, bps, mprv, seq+1);
         else
            *seq = 0;
      }
   }
/*
 * If we reach here, no sequence got everything we needed, so return largest
 * span we found, whose series pattern indices are encoded in seq
 */
   seq[1] = 0;
   return(msk);
}

uint findMaxSpan1(uint N, uint mcoh, uint fi, uint idx, uint *cnts, uint **bps,
                  uint mprv, uint *seq)
/*
 * searches thru all files for maximum span using a greedy algorithm on span
 */
{
   const uint mall=(1<<N)-1;
   uint msk, mneed;
   msk = bps[fi][idx+1] & mall;
   mneed = mcoh ^ msk;  /* remove parts of mcoh we've found so far */
   *seq = msk;
   if (mneed)
   {
      uint k;
/*
 *    Look for files where bit is set in mcoh, but not msk
 */
      for (k=0; k < N; k++)
      {
         if (!(mneed&(1<<k)))
            continue;
         msk = findMaxSpanR(N, mcoh, msk, mneed, cnts, bps, mprv, seq+1);
      }
   }
   else
      seq[1] = 0;
   return(msk);
}

uint findMaxSpan(uint N, uint mcoh, uint id, uint *cnts, uint **bps, uint mprv,
                 uint *seq)
/*
 * Maximize the number of matching set bits with mcoh using any bit pattern
 * in file id.  Return immediately if we can match all set bits.
 */
{
   uint *sq0, *sq;
   uint ncoh, n, msk, maxset=0, maxmsk=0, i;
/*
 * sq0 will be used to find possible sequences, with best copied to seq
 */
   sq0 = malloc(N*sizeof(uint));
   assert(sq0);
   *sq0 = 0;          /* sq0 starts out empty */
   sq = bps[id] + 1;
   n = sq[-1];
   assert(n > 0);
   for (i=0; i < n; i++)
   {
      uint ns, NEW=0;
      msk = findMaxSpan1(N, mcoh, id, i, cnts, bps, mprv, sq0);
      ns = NumSetBits(N, msk);
      if (ns > maxset)
         NEW = 1;
      else if (ns == maxset)
      {
         uint m0 = (maxmsk|mprv)&mcoh, m1 = (msk|mprv)&mcoh;
         NEW = (NumSetBits(N, m0) < NumSetBits(N, m1));
      }
      if (NEW)
      {
         uint j;
         for (j=0; (seq[j] = sq0[j]); j++); /* copy seq including NULL end */
         maxmsk = msk;
         maxset = NumSetBits(N, msk);
      }
      if (msk == mcoh) /* fully spanned, no need to keep searching */
         break;
   }
   free(sq0);
   return(maxmsk);
}

void setCopyOrder(const uint N,  ATL_cpnode_t **cbs)
/*
 * Given N copy files that need to cohere (have the same routine at same idx)
 * Sort into an order that minimizes the number of NULLs that must be stored.
 */
{
   ATL_cpnode_t **newb, **prvb, *cb, *cp;
   const uint mall=(1<<N)-1; /* mask fo all N bits set */
   uint mcoh=mall, mprv=mall;  /* assume all files must cohere at beginning */
   uint i, mincnt=0, imin=0, n=N;
   uint *cnts, *seq;
   uint **bps; /* coherence bit patterns, count 1st entry */
/*
 * Create bitvec showing what all files each node must show up in
 */
   DoShareBV(N,  cbs);
   newb = malloc((N+N)*(sizeof(ATL_cpnode_t*)+sizeof(uint)));
   assert(newb);
   prvb = newb + N;
   cnts = (uint*) (prvb + N);
   seq = cnts + N;   /* sequence of patterns used to span */
   newb[0] = cb = CPmoveMaskBVNewQ(mall, mall, cbs);
/*
 * If we've got copy kerns appearing in all files, that will be our first
 * idx entries and they are guaranteed to cohere, so just let 1st queue
 * dictate order to all others
 */
   if (cb)
   {
      for (i=1; i < N; i++)
      {
         ATL_cpnode_t *cp;
         cp = CPmoveMaskBVNewQ(mall, mall, cbs+i);
         newb[i] = CopyDictateOrder(cb, cp);
      }
   }
/*
 * If no copy appears everywhere, start with empty sorted queues
 */
   else
      for (i=1; i < N; i++)
         newb[i] = NULL;
/*
 * Any copy appearing in only one file can't cause a coherence problem,
 * any remaining after plugging in NULLs can just be dumped at end of file
 */
   for (i=0; i < N; i++)
      prvb[i] = CPmoveMaskBVNewQ(mall, 1<<i, cbs+i);
/*
 * The shortest queue of shared files will be finished first, in order
 * to allow us to reduce N, and ensure the short file doesn't need a large
 * number of NULLs to keep matching larger files
 */
   for (mincnt=imin=i=0; i < N; i++)
   {
      uint cnt;
      cnts[i] = cnt = ATL_CountNumberOfCPNodes(cbs[i]);
      if (cnt == 0)
      {
         if (prvb[i])
            newb[i] = ATL_JoinCPQs(newb[i], prvb[i]);
         prvb[i] = NULL;
         n--;
         mcoh ^= 1<<i;
      }
      else if (!mincnt || cnt < mincnt)
      {
         mincnt = cnt;
         imin = i;
      }
   }
/*
 * Construct list of all coherence bit patterns that appear in each active file
 */
   bps = malloc(sizeof(uint*)*N);
   assert(bps);
   for (i=0; i < N; i++)
   {
      if (cnts[i] && mcoh & (1<<i))
      {
         uint np, j;
         uint *bp;
/*
 *       Going to store npatt in 1st entry, can use bitpattern 0 for binary
 */
         np = (1<<n);  /* number of diff bit patterns based on binary */
         j = cnts[i] + 1; /* can't have more than copies! */
         np = Mmin(j, np);
         bp = malloc(sizeof(uint)*np);
         assert(bp);
         bps[i] = bp++;
         j = coherePatt2Arr(cbs[i], bp);
         assert(j <= np);
         bp[-1] = j;
         if (i == 5 && bps[i])
         {
            uint k;
            fprintf(stderr, "bps[%d]={(%u,%x)", j, bp[0]>>N, mall&bp[0]);
            for (k=1; k < j; k++)
               fprintf(stderr, ",(%u,%x)\n", bp[k]>>N, mall&bp[k]);
            fprintf(stderr, "}\n");
         }
         addSetCnt2Patt(N, j, bp);
         sortPattBySpan(N, j, bp);  /* sort greatest-to-least span */
      }
      else
         bps[i] = NULL;
   }
/*
 * While we've still got files that must cohere
 */
   while (mcoh)
   {
      uint ns, id;  /* index of draw file */
      uint b, bv, nspan, mprv, mneed;
      if (!(mcoh & (mcoh-1)))
      {
         uint j;
         fprintf(stderr, "mcoh=%u,%x\n", mcoh, mcoh);
         fprintf(stderr, "ncnt={%u", cnts[0]);
         for (j=1; j < N; j++)
            fprintf(stderr, ",%u", cnts[j]);
         fprintf(stderr, "}\n");
         fprintf(stderr, "bps={%d", bps[0]?bps[0][0]:-1);
         for (j=1; j < N; j++)
            fprintf(stderr, ",%d", bps[j]?bps[j][0]:-1);
         fprintf(stderr, "}\n");
         for (j=0; !(mcoh&(1<<j)); j++);
         fprintf(stderr, "bps[%u][1]=%u,%x\n", j, bps[j][1], bps[j][1]);
      }
      assert((mcoh & (mcoh-1))); /* impossible to cohere wt no other files! */
      id = selDrawFile(N, cnts, bps);  /* find file to draw from */
      for (mprv=i=0; i < N; i++)       /* mprv bitvec showing if we have */
         if (prvb[i])                  /* private kerns that can be used */
            mprv |= 1<<i;              /* to plug NULLs; used as tiebreak */
      bv = findMaxSpan(N, mcoh, id, cnts, bps, mprv, seq);
/*
 *    Loop over all bit patterns, and remove the kerns involved in them
 */
      for (i=0; seq[i]; i++)
      {
         const uint bp0=seq[i];
         uint bp=bp0, j=0;

/*
 *       Now, make changes to all files metadata & kern Q involved in bp
 */
         do
         {
            ATL_cpnode_t *nxtp;
            while (!(bp&(1<<j)))  /* find next file set in bit pattern bp */
               j++;
            cp = CPfindMaskedBV(cbs[j], mall, bp0);
            if (cp == NULL)
            {
               uint np = bps[j] ? bps[j][0] : 0;
               uint ii;
               fprintf(stderr, "cbs[%u]=%p, flg=%x\n", j, cbs[j],
                       ((cbs[j]->bv)&mall));
               fprintf(stderr, "N=%u, j=%u, bp=%x,%x\n", N, j, bp0, bp);
               fprintf(stderr, "bps[%u]={%u ", j, np);
               for (ii=1; ii <= np; ii++)
                  fprintf(stderr, ",%x", bps[j][ii]);
               fprintf(stderr, "}\n");
               assert(cp);
            }
            nxtp = cp->next;
            cbs[j] = RemoveCPNodeFromQ(cbs[j], cp);
            newb[j] = ATL_JoinCPQs(newb[j], cp);
            if (--cnts[j] == 0)  /* 1 less node in file, if 0 */
            {
               mcoh ^= 1<<j;     /* no longer have to cohere wt others*/
               if (prvb[j])
                  newb[j] = ATL_JoinCPQs(newb[j], prvb[j]);
               n--;
               assert(bps[j][0] == 1);
               bps[j][0] = 0;
            }
/*
 *          If this was the last entry of that bit pattern in this file,
 *          we must reduce the number bit patterns!
 */
            else if (!CPfindMaskedBV(nxtp, mall, bp0))
            {
               uint *ip = bps[j];
               uint np = *ip;
               uint k;
               assert(np > 1);
               *ip = --np;
               ip++;
               for (k=0; (ip[k]&mall) != bp0; k++);
               assert(k <= np);
               for (; k < np; k++)
                  ip[k] = ip[k+1];
            }
            bp ^= 1<<j;
            j++;
         }
         while (bp);
      }
/*
 *    If we still need kernels for bits not set in the span (bv), try filling
 *    from prvb before adding NULL
 */
      mneed = mcoh & (~bv);
      if (mneed)
      {
         for (i=0; i < N; i++)
         {
            if (mneed&(1<<i))
            {
               cp = prvb[i];
               if (cp)
               {
                  prvb[i] = cp->next;
                  cp->next = NULL;
               }
               else
                  cp = GetCPNode();
               newb[i] = ATL_JoinCPQs(newb[i], cp);
            }
         }
      }
   }
/*
 * Free coherence bit pattern arrays
 */
   for (i=0; i < N; i++)
      if (bps[i])
         free(bps[i]);
   free(bps);
/*
 * Put sorted list back into input array
 */
   for (i=0; i < N; i++)
      cbs[i] = newb[i];
   free(newb);
}

uint *CountMoveShared(const uint N,  ATL_cpnode_t **cbs)
{
   uint *cnts;
   ATL_cpnode_t *cp=cbs[0];
   uint i, ioff, nxtoff;

   cnts = malloc(N*sizeof(uint));
   assert(cnts);
   ioff = GetOffset(&cp->flag, cp);
   nxtoff = GetOffset(&cp->next, cp);
   for (i=0; i < N; i++)
   {
      uint j;
/*
 *    Mark all files duplicated between i and remaining lists
 */
      for (j=i+1; j < N; j++)
         CopyMarkDup_cohere(cbs[i], cbs[j]);
/*
 *    Move all duplicated entries to beginning of list
 */
      cbs[i] = CopySortMarkedFirst(cbs[i]);
      cnts[i] = CountListMaskALL(cbs[i], (1<<CPF_TMP), nxtoff, ioff);
   }
   return(cnts);
}

void SortByCount(uint N, uint *nshar,  ATL_cpnode_t **cbs, uint *flgs)
{
   uint i;
   for (i=0; i < N-1; i++)
   {
      uint imin=i, min=nshar[i], j;
      for (j=i+1; j < N; j++)
      {
         if (nshar[j] < min)
         {
            imin = j;
            min = nshar[j];
         }
      }
      if (imin != i)
      {
         void *vp = cbs[i];
         j = imin;
         cbs[i] = cbs[j];
         cbs[j] = vp;
         nshar[j] = nshar[i];
         nshar[i] = min;
         min = flgs[i];
         flgs[i] = flgs[j];
         flgs[j] = min;
      }
   }
}

void DictateOrders(const uint N, ATL_cpnode_t **cbs)
{
   uint i;
   for (i=0; i < N-1; i++)
   {
      uint j;
      for (j=i+1; j < N; j++)
         cbs[j] = CopyDictateOrder(cbs[i], cbs[j]);
   }
}

void CohereLists(const uint N, ATL_cpnode_t **cbs, uint *flgs)
{
   int n, i=0;
   do
   {
      unsigned int *nshar;

      n = CountRegion(N-i, cbs);
      #if 0
         nshar = CountMoveShared(n, cbs);
         SortByCount(n, nshar, cbs, flgs);
         free(nshar);
         DictateOrders(n, cbs);
      #else
         setCopyOrder(n,  cbs);
      #endif

      cbs += n;
      flgs += n;
      i += n;
   }
   while (i < N);
}
void CreateMasterIdx(char pre)
{
   ATL_cpnode_t *cb, **cbs, *cp;
   FILE *fp;
   unsigned int *flgs;
   unsigned int N, i;
   const unsigned int msk=~((1<<CPF_MVEC)|(1<<CPF_NVEC)|(1<<CPF_ASM)|
      (1<<CPF_TMP));
   char fn[32]={'r','e','s','/',pre,'m','a','s','t','e','r','.','C','P','I',0};

   cb = ReadCPFileWithPath(pre, "res", "cpyPERF.CPS");
   assert(cb);
   for (cp=cb; cp; cp = cp->next)
   {
      if (cp->rout)
         free(cp->rout);
      cp->rout = GetCopyName(cp, 0);
   }
   cbs = SplitList(pre, cb, &N, &flgs);
   KillAllCPNodes(cb);
   CohereLists(N, cbs, flgs);
/*
 * Write out master list of used copy kernels
 */
   fp = fopen(fn, "w");
   assert(fp);
   fprintf(fp, "%u\n", N);
   for (i=0; i < N; i++)
   {
      char *sp;
      unsigned int flag;
      cb = cbs[i];
      flag = cb->flag & msk;
      sp = CopyFlag2Str(flag);
      fprintf(fp, "%x '%s'\n", flag, sp);
      sprintf(fn+5, "cpyPERF_%s.CPS", sp);
      WriteCPFile(fn, cb);
   }
   fclose(fp);

   for (i=0; i < N; i++)
      KillAllCPNodes(cbs[i]);
   free(cbs);
   free(flgs);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   fprintf(stderr,
"This program is a driver for searching & generating all copy kernels.\n"
"Output file of cumulative views: res/<pre>full.CPS, used for xcpygen.\n"
"Delete this file if you get kerns you don't need in a prior search.\n"
"You can either invoke this driver routine multiple times, or repeat\n"
"the -V flag.\n");

   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr,"USAGE: %s [flags:\n", name);
   fprintf(stderr,"   -p [d,s,c,z]\n");
   fprintf(stderr, "   -o /path : specify directory to generate into\n");
   fprintf(stderr,
"   -V [C,A]a=1,N,X Cb=0,1,N,X [A,C]d=I,F K=0,1 S=A,C,M,U <name> mmview.sum:\n"
          );
   fprintf(stderr,
"      If args to left of = are omitted, these vals don't appear in the view.\n"
"      Args to right of = may be omitted, which means all values are set.\n"
"      <name> and <mmview> must appear in that order.\n"
"      A/Cd=[I,F]: copy Into or From col-major (I,F means both directions).\n"
"      A/Ca/b give the list of needed alpha/beta for that matrix copy.\n"
"      K : change how K-cleanup is done, default 1:\n"
"          0: Kerns self-clean, so no padding if MOD(kbB,ku)==0, else ku-pad\n"
"          1: Kernels use other kernels for K-cleanup.\n"
"      S=A,C,M,U: don't store [A,C] copyI, matmulI, or unrollings, resp.\n"
"      Suppressed copies will still be added to master list, just not in view\n"
"      <name> is a unique string that will be used in all header files.\n"
"      mmview.sum: list of amm kerns demanding the copies.\n");
   exit(ierr ? ierr : -1);
}

ATL_view_t *GetFlags
   (int nargs, char **args, char *PRE, char **PTH)
{
   ATL_view_t *vb=NULL, *vp;
   char *pth=NULL;
   int i, minSz=24*24;
   char pre='d';

   for (i=1; i < nargs; i++)
   {
      char *nam;
      int k, flag;
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         pre = tolower(args[i][0]);
         assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
         break;
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         pth = DupString(args[i]);
         break;
      case 'V':    /* -V C/Aa=1,N,X Cb=0,1,N,X C/Ad=I,F K=[0,1] S=A,C,M,U */
         if (++i >= nargs)/*    <name> mmview.sum */
            PrintUsage(args[0], i-1, NULL);
         flag = 0;
         for (k=0; k < 7; k++)
         {
            char mt, wt;
            int st;
            mt = args[i][0];
            if (mt == 'S')
            {
               char *cp = args[i] + 1;
               assert(*cp == '=');
               do
               {
                  cp++;
                  switch(*cp)
                  {
                  case 'A':
                     flag |= 1<<CPV_NOACP;
                     break;
                  case 'C':
                     flag |= 1<<CPV_NOCCP;
                     break;
                  case 'M':
                     flag |= 1<<CPV_NOMMI;
                     break;
                  case 'U':
                     flag |= 1<<CPV_NOUNR;
                     break;
                  default:
                     PrintUsage(args[0], i-1, NULL);
                  }
                  cp++;
               }
               while(*cp == ',');
            }
            else if (mt == 'K')
            {
               char *cp = args[i] + 1;
               assert(*cp == '=');
               if (cp[1] == '0')
                  flag |= 1<<CPV_SELFK;
               else if (cp[1] == '1')
                  flag &= ~(1<<CPV_SELFK);
               else
                  PrintUsage(args[0], i-1, NULL);
            }
            else if (mt != 'C' && mt != 'A')
               break;
            else
            {
               wt = args[i][1];
               if (wt != 'a' && wt != 'b' && wt != 'd')
                  break;
               if (args[i][2] != '=')
                  break;
               if (wt == 'b')
               {
                  assert(mt == 'C');
                  flag |= CPV_ScalStr2bits(args[i]+3, CPV_BE1C);
               }
               else if (wt == 'a')
                  flag |= CPV_ScalStr2bits(args[i]+3,
                          (mt=='C')?CPV_AL1C:CPV_AL1A);
               else if (wt == 'd')
                  flag |= CPV_DirStr2bits(args[i]+3,
                          (mt=='C')?CPV_C2BLK:CPV_A2BLK);
               else
                  PrintUsage(args[0], i-1, NULL);
            }
            if (++i >= nargs)
               PrintUsage(args[0], i-1, NULL);
         }
         nam = DupString(args[i]);
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         vp = ATL_NewView(flag, nam, DupString(args[i]));
         vp->next = vb;
         vb = vp;
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   assert(vb);
   if (!pth)
      pth = DupString("tmp");
   *PTH = pth;
   *PRE = pre;
   return(vb);
}


void DoGen(char pre, char *path)
{
   char *ln;
   ln = NewMergedString("./xcpygen -i res/XcpyPERF.CPS -o ", path);
   ln[17] = pre;
   assert(!Sys2File(ln, NULL));
   free(ln);
}

void DoGenH(char pre, char *path, ATL_view_t *vb)
{
   ATL_view_t *vp;
   int plen;

   plen=strlen(path);
   for (vp=vb; vp; vp = vp->next)
   {
      char *ln, *sp, *va;
      int len, d;

      va = View2Args(vp);
      len = 20 + strlen(va) + plen;
      ln = malloc(len*sizeof(char));
      assert(ln);
      d = sprintf(ln, "./xcphgen -p %c %s -o %s", pre, va, path);
      assert(d < len);
      assert(!Sys2File(ln, NULL));
      free(va);
      free(ln);
   }
}

/*
 * xcpydrv -Vs -o cpylst.CPS | xcpysearch -o cpyPERF.CPS | xcpygen/xcphgen
 */
int main(int nargs, char **args)
{
   ATL_view_t *vb=NULL, *vp;  /* view base & ptr */
   ATL_cpnode_t *cb=NULL, *op;   /* cpsearch output files */
   char *path, *ln;
   char pre;

   vb = GetFlags(nargs, args, &pre, &path);
   CreateMasterIdx(pre);
   DoGen(pre, path);
   DoGenH(pre, path, vb);
   KillAllViews(vb);
   free(path);
   return(0);
}
