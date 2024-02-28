#ifndef ATLAS_MMPARSE_H
   #define ATLAS_MMPARSE_H

#include "atlas_genparse.h"
#include "atlas_enum.h"
/*
 * BLAS kernel types that we've potentially got generators for
 */
#define ATL_KGEMM       0
#define ATL_KSYRK       1
#define ATL_KSYMM       2 /* not yet supported */
#define ATL_KTRMM       3 /* not yet supported */
#define ATL_KTRSM       4 /* only supported for non-generated case */
#define ATL_KGECPFA     5 /* TA trans setting, TB row/col matrix access */
#define ATL_KGECP2A     6 /* TA trans setting, TB row/col matrix access */
#define ATL_KGECPFC     7
#define ATL_KGECP2C     8
#define ATL_KSKCPFC     9 /* syrk C copy */
#define ATL_KSKCP2C    10 /* syrk C copy */
/*
 * Structure bit position in flag
 */
#define MMF_SINGLE      0  /* 1: single precision, else double */
#define MMF_COMPLEX     1  /* 1: complex type, else real */
#define MMF_LDCTOP      2  /* 1: load C before K-loop, 0: ld after */
#define MMF_KRUNTIME    3  /* 1: K dim is run-time variable? */
#define MMF_X87         4  /* 1: requires the Intel x87 unit */
#define MMF_AOUTER      5  /* 1: MNK loop order, 0: NMK loop order */
#define MMF_PFACOLS     6  /* 1: prefetch next mu cols of A */
#define MMF_PFABLK      7  /* 1: prefetch next KBxNB block of A */
#define MMF_PFBCOLS     8  /* 1: prefetch next nu cols of B */
#define MMF_PFCELTS     9  /* 1: pf elts of C at top of loop, load at bottom */
#define MMF_L14NB      10  /* 1: need to fit all 3 matrices+nextA in L1 */
#define MMF_MVA        11  /* 1: A expected to change between calls */
#define MMF_MVB        12  /* 1: B expected to change between calls */
#define MMF_MVC        13  /* 1: C expected to change between calls */
#define MMF_CPC        14  /* 1: copy out to col of C like rank-K */
#define MMF_KVEC       15  /* 1: vectorize on K-dim. 0: do not vec on K dim */
#define MMF_KUISKB     16  /* 1: only works for fully unrolled, constant K */
#define MMF_NOBCAST    17  /* 1: use ld/splat to get B, rather than bcast */
#define MMF_FKO        18  /* 1: compile using FKO, else [d,s]MC */
/*
 * Right now, NOBCAST overrides BREG1.  Could eventually make it so
 * BREG1+NOBCAST means we do splat wt 0 dep distance, but not supported now
 */
#define MMF_BREG1      19  /* 1: use only 1 reg for B, else use NU */
#define MMF_KCLN       20  /* 1: cleaner for K dim is required */
#define MMF_ALLTRANS   21  /* 1: do C^T = B^T * A^T, not C=A*B */
/*
 * Some definitions for triangular kernels
 */
#define MMF_RIGHT      22  /* 1: Right, else Left */
#define MMF_UPPER      23  /* 1: Upper, else Lower */
#define MMF_MAXBIT     23
#define MMF_ALLPF ( (1<<MMF_PFACOLS)|(1<<MMF_PFABLK)|(1<<MMF_PFBCOLS) \
                    |(1<<MMF_PFCELTS) )
#define MMF_MVSET  ( (1<<MMF_MVA) | (1<<MMF_MVB) | (1<<MMF_MVC) | (1<<MMF_CPC) )
#define MMF_MVDEF  ( (1<<MMF_MVA) | (1<<MMF_MVB) )
#define MMF_DIFFMSK (~((1<<MMF_KCLN)|MMF_MVSET|(1<<MMF_RIGHT)|(1<<MMF_UPPER)))
#define MMF_STKBTS (1<<MMF_KCLN)
#define MMF_DEFAULT ( (1<<MMF_AOUTER) | MMF_MVDEF )
#ifndef  FLAG_IS_SET
   #define FLAG_IS_SET(field_, bit_) ( ((field_) & (1<<(bit_))) != 0 )
#endif
#define ATL_MMF_MVGET(field_) (((field_) >> MMF_MVA) & 0xF)
#define ATL_MMF_MVPUT(field_, v_) \
   (field_) = ( ((field_) & ~MMF_MVSET) | (((v_) & 0xF) << MMF_MVA) )

typedef struct MMStg ATL_mmstg_t;
struct MMstg
{
   int extelts;
   short ID, flag;
   char *szfunc;
};

#define uchar unsigned char
#define uint  unsigned int
typedef struct MMNode ATL_mmnode_t;
struct MMNode
{
   ATL_mmnode_t *next;
   double mflop[8];             /* 1st entry perf using mbB, nbB, kbB */
   char *rout, *auth, *comp, *cflags;
   char *str;                   /* tmp string used in generation */
   char *genstr;                /* system(genstr) will generate gened kernel */
   char *exflags;               /* extra flags to pass test/time call */
   char *moves;                 /* -DMove[A,B,C] to use, NULL default */
   int ID, mu, nu, ku;          /* ID, and unrolling on each loop */
   int kbmin, kbmax;            /* min/max KB this kernel can handle */
   int SSE;                     /* 0: no SSE, 1: SSE1 req, 2: SSE2 req, etc */
   int vlen;                    /* vector length, 0 or 1 if scalar code */
   int pref;                    /* pref strategy used */
   int pfLS;                    /* pref line size */
   int ivar;                    /* int param used for various purposes */
   int mbB, nbB, kbB;           /* best blocking dims found by search */
   int szA, szB, szC;           /* size in elts for block; 0: use normal */
   int szExtra;                 /* extra elts needed at end of alloc */
   int stgA, stgB, stgC;        /* storage ID, 0: access major */
   enum ATLAS_TRANS TA, TB;
   int asmbits;   /* bitfield indicating which assembly(ies) is required */
   uint flag;
   uchar blask;   /* 0:ammm, 1:syrk, 2:symm, 3:trmm */
};
#undef uchar
#undef uint

#ifndef ATL_DEF_MMFLAG
   #define ATL_DEF_MMFLAG @up@(rt)F_DEFAULT
#endif
/* procedure 1 */
static ATL_mmnode_t *GetMMNode(void)
{
   ATL_mmnode_t *p;
   p = calloc(1, sizeof(ATL_mmnode_t));
   assert(p);
   p->TA = AtlasTrans; p->TB = AtlasNoTrans;
   p->flag = MMF_DEFAULT;
   return(p);
}

/* procedure 2 */
static void CopyMMNode(ATL_mmnode_t *p, ATL_mmnode_t *dup)
{
   ATL_mmnode_t *nxt=p->next;
   if (p->str)
      free(p->str);
   if (p->exflags)
      free(p->exflags);
   if (p->genstr)
      free(p->genstr);
   if (p->rout)
      free(p->rout);
   if (p->cflags)
      free(p->cflags);
   if (p->comp)
      free(p->comp);
   if (p->auth)
      free(p->auth);
   memcpy(p, dup, sizeof(ATL_mmnode_t));
   if (dup->str)
      p->str = DupString(dup->str);
   if (dup->exflags)
      p->exflags = DupString(dup->exflags);
   if (dup->genstr)
      p->genstr = DupString(dup->genstr);
   if (dup->rout)
      p->rout = DupString(dup->rout);
   if (dup->cflags)
      p->cflags = DupString(dup->cflags);
   if (dup->comp)
      p->comp = DupString(dup->comp);
   if (dup->auth)
      p->auth = DupString(dup->auth);
   if (dup->moves)
      p->moves = DupString(dup->moves);
   p->next = nxt;
}
/* procedure 3 */
static ATL_mmnode_t *CloneMMNode(ATL_mmnode_t *dup)
{
   ATL_mmnode_t *p;
   p = malloc(sizeof(ATL_mmnode_t));
   assert(p);
   memcpy(p, dup, sizeof(ATL_mmnode_t));
   if (dup->str)
      p->str = DupString(dup->str);
   if (dup->exflags)
      p->exflags = DupString(dup->exflags);
   if (dup->genstr)
      p->genstr = DupString(dup->genstr);
   if (dup->rout)
      p->rout = DupString(dup->rout);
   if (dup->cflags)
      p->cflags = DupString(dup->cflags);
   if (dup->comp)
      p->comp = DupString(dup->comp);
   if (dup->auth)
      p->auth = DupString(dup->auth);
   if (dup->moves)
      p->moves = DupString(dup->moves);
   p->next = NULL;
   return(p);
}

/* procedure 4: clones a queue of MM structs */
static ATL_mmnode_t *CloneMMQueue(ATL_mmnode_t *dupb)
{
   ATL_mmnode_t *p, *pd, *nb;
   if (!dupb)
      return(NULL);
   p = nb = CloneMMNode(dupb);
   for (pd=dupb->next; pd; pd = pd->next)
   {
      p->next = CloneMMNode(pd);
      p = p->next;
   }
   return(nb);
}

/* procedure 5: clones a queue of strided MM structs */
static ATL_mmnode_t *CloneStridedMMQueue
(
   ATL_mmnode_t *dupb,   /* queue of nodes to clone */
   int stride               /* increment between nodes to take */
)
/*
 * Creates a queue of cloned nodes from dupb; move stride each time
 * (stride must be >= 1); i.e. skip stride-1 structs in original queue
 */
{
   ATL_mmnode_t *p, *pd, *nb;
   int i;

   if (!dupb)
      return(NULL);
   if (stride == 1)
      return(CloneMMQueue(dupb));
   assert(stride > 1);
   p = nb = CloneMMNode(dupb);
   pd = nb;
   while(pd)
   {
      for (i=0; i < stride && pd; i++, pd = pd->next);
      if (pd)
      {
         p->next = CloneMMNode(pd);
         p = p->next;
      }
      else
         p->next = NULL;
   }
   return(nb);
}

/* procedure 6 */
static ATL_mmnode_t *KillMMNode(ATL_mmnode_t *die)
{
   ATL_mmnode_t *p=NULL;
   if (die)
   {
      p = die->next;
      if (die->rout)
         free(die->rout);
      if (die->auth)
         free(die->auth);
      if (die->comp)
         free(die->comp);
      if (die->cflags)
         free(die->cflags);
      if (die->str)
         free(die->str);
      if (die->genstr)
         free(die->genstr);
      if (die->exflags)
         free(die->exflags);
      if (die->moves)
         free(die->moves);
      free(die);
   }
   return(p);
}

/* procedure 7: safely remove nukeme from Q, reseting all links */
static ATL_mmnode_t *RemoveMMNodeFromQ
(
   ATL_mmnode_t *Q,     /* queue of nodes */
   ATL_mmnode_t *nukeme /* node to remove from queue */
)
/*
 * Removes nukeme from Q, sets nukeme->next=NULL, and returns updated Q
 */
{
   ATL_mmnode_t *p, *prev;

   if (!nukeme)
      return(Q);
   assert(Q);
   if (Q == nukeme)
   {
      Q = Q->next;
      nukeme->next = NULL;
      return(Q);
   }
   prev = Q;
   for (p=Q->next; p && p != nukeme; p = p->next)
      prev = p;
   assert(p);
   prev->next = nukeme->next;
   nukeme->next = NULL;
   return(Q);
}

/* procedure 8: swap two nodes in Q */
static ATL_mmnode_t *SwapMMNodeInQ
(
   ATL_mmnode_t *Q,     /* queue of nodes */
   ATL_mmnode_t *p0,    /* first node to be swap */
   ATL_mmnode_t *p1     /* second node to be swap */
)
/*
 * Puts p0 in queue at place p1 was, and p1 in place p0 was.
 * RETURNS: possibly changed Q
 */
{
   ATL_mmnode_t *p, *prv0=NULL, *prv1=NULL, *nxt;
   unsigned int fnd=0;

   if (!Q || !p0 || !p1 || p1 == p0) /* need different pair to swap */
      return(Q);
   for (p=Q; fnd != 3 && p; p = p->next)
   {
      if ((fnd&1) == 0) /* still looking for p0 */
      {
         if (p == p0)
            fnd |= 1;
         else
         {
            prv0 = p;
         }
      }
      if ((fnd&2) == 0) /* still looking for p1 */
      {
         if (p == p1)
            fnd |= 2;
         else
         {
            prv1 = p;
         }
      }
   }
   if (fnd != 3)
   {
      fprintf(stderr, "\nFND=%u!!\n\n", fnd);
      assert(fnd == 3);
   }
   if (prv0 == p1)  /* p1->p0-> */
   {
      if (prv1)
         prv1->next = p0;
      p1->next = p0->next;
      p0->next = p1;
   }
   else if (prv1 == p0) /* p0->p1-> */
   {
      if (prv0)
         prv0->next = p1;
      p0->next = p1->next;
      p1->next = p0;
   }
   else /* neither predecessor of other */
   {
      nxt = p0->next;
      p0->next = p1->next;
      p1->next = nxt;
      if (prv0)
         prv0->next = p1;
      if (prv1)
         prv1->next = p0;
   }
   if (Q == p1)
      Q = p0;
   else if (Q == p0)
      Q = p1;
   return(Q);
}

/* procedure 9 */
static ATL_mmnode_t *KillMMNodeFromQ
(
   ATL_mmnode_t *Q,     /* queue of nodes */
   ATL_mmnode_t *nukeme /* node to remove from queue */
)
{
   Q = RemoveMMNodeFromQ(Q, nukeme);
   KillMMNode(nukeme);
   return(Q);
}

/* procedure 10: frees all ->rout entries in queue b */
static void KillAllMMRouts(ATL_mmnode_t *b)
{
   while (b)
   {
      if (b->rout)
      {
         free(b->rout);
         b->rout = NULL;
      }
      b = b->next;
   }
}

/* procedure 10 */
void *KillAllMMStrings(ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;

   for (mp=mb; mp; mp = mp->next)
   {
      if (mp->moves)
         free(mp->moves);
      mp->moves = NULL;
      if (mp->exflags)
         free(mp->exflags);
      mp->exflags = NULL;
      if (mp->genstr)
         free(mp->genstr);
      mp->genstr = NULL;
      if (mp->str)
         free(mp->str);
      mp->str = NULL;
      if (mp->cflags)
         free(mp->cflags);
      mp->cflags = NULL;
      if (mp->comp)
         free(mp->comp);
      mp->comp = NULL;
      if (mp->auth)
         free(mp->auth);
      mp->auth = NULL;
      if (mp->rout)
         free(mp->rout);
      mp->rout = NULL;
   }
}

/* procedure 11: frees all ->genstr entries in queue b */
static void KillAllMMGenstr(ATL_mmnode_t *b)
{
   while (b)
   {
      if (b->genstr)
      {
         free(b->genstr);
         b->genstr = NULL;
      }
      b = b->next;
   }
}

/* procedure 12 */
static void KillAllMMNodes(ATL_mmnode_t *die)
{
   while (die)
      die = KillMMNode(die);
}

/* procedure 13 */
static void ATL_SubGoodGccInMMNodes
(
   ATL_mmnode_t *bp   /* queue to make sub in */
)
/*
 *  Gets GOODGCC (from Make.inc), and substitutes it for all comp == "gcc"
 *  in the queue.  This gets us mandatory flags like -pg,-m64,etc.
 */
{
   ATL_mmnode_t *kp;  /* queue to make sub in */
   char *gcc;
   gcc = GetGoodGcc();
   for (kp=bp; kp; kp = kp->next)
   {
      if (kp->comp && !strcmp(kp->comp, "gcc"))
      {
         free(kp->comp);
	 kp->comp = DupString(gcc);
      }
   }
}

/* procedure 14 */
static void ATL_UnsubGoodGccInMMNodes
(
   ATL_mmnode_t *bp   /* queue to make reverse sub in */
)
/*
 *  Gets GOODGCC (from Make.inc); Any comp string matching that is switched
 *  back to "gcc".  This is usually necessary so that output files don't
 *  use an old GOODGCC that lacks something like -pg.
 */
{
   ATL_mmnode_t *kp;  /* queue to make sub in */
   char *gcc;
   gcc = GetGoodGcc();
   for (kp=bp; kp; kp = kp->next)
   {
      if (kp->comp && !strcmp(kp->comp, gcc))
      {
         free(kp->comp);
	 kp->comp = DupString("gcc");
      }
   }
}

/* procedure 15 */
static void ResubGoodGccInMMNodes
(
   ATL_mmnode_t *bp   /* queue to make sub in */
)
/*
 * Takes gcc compiler that use GOODGCC, and replaces them with "gcc"
 * to help portability
 */
{
   ATL_mmnode_t *kp;  /* queue to make sub in */
   char *gcc;
   gcc = GetGoodGcc();
   for (kp=bp; kp; kp = kp->next)
   {
      if (kp->comp && !strcmp(kp->comp, gcc))
      {
         free(kp->comp);
	 kp->comp = DupString("gcc");
      }
   }
}

/* procedure 16 */
static int ATL_CountNumberOfMMNodes
(
    ATL_mmnode_t *bp   /* queue to count */
)
{
   int i;
   for (i=0; bp; i++, bp = bp->next);
   return(i);
}

/* procedure 17 */
static ATL_mmnode_t *ATL_LastMMNode(ATL_mmnode_t *bp)
/*
 * RETURNS: pointer to last node in queue
 */
{
   ATL_mmnode_t *p;
   if (!bp)
      return(NULL);
   for (p=bp; p->next; p = p->next);
   return(p);
}

/* procedure 18: adds q1 to end of q0 */
static ATL_mmnode_t *ATL_JoinMMQs(ATL_mmnode_t *q0, ATL_mmnode_t *q1)
{
   ATL_mmnode_t *mp;
   if (!q1)
      return(q0);
   if (!q0)
      return(q1);
   mp = ATL_LastMMNode(q0);
   mp->next = q1;
   return(q0);
}

/* procedure 19: finds node with max mflop[imf]  */
static ATL_mmnode_t *FindMaxMflopMMQ
(
   ATL_mmnode_t *bp,   /* queue to be searched */
   int imf
)
/*
 * RETURNS: ptr to structure containing max value in mflop[imf]
 */
{
   ATL_mmnode_t *mp, *p;
   double mfm;
   if (!bp)
      return(NULL);
   mp = bp;
   mfm = mp->mflop[imf];
   for (p=bp->next; p; p = p->next)
   {
      const double mf=p->mflop[0];
      if (mf > mfm)
      {
         mfm = mf;
         mp = p;
      }
   }
   return(mp);
}
/* procedure 20: finds max integer at ip0 in struct */
static ATL_mmnode_t *FindMaxIntInMMQ
(
   ATL_mmnode_t *bp,   /* queue to be searched */
   void *ip0           /* ptr to integer withinin node bp */
)
/*
 * RETURNS: ptr to structure containing max int value at byte offset
 *          offset in struct
 */
{
   ATL_mmnode_t *mp=NULL, *p;
   int *ip;
   int val;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   if (!bp)
      return(NULL);

   mp = bp;
   ip = (int*)(((char*)bp) + offset);
   val = *ip;
   for (p=bp->next; p; p = p->next)
   {
      ip = (int*)(((char*)p) + offset);
      if (*ip > val)
      {
         mp = p;
         val = *ip;
      }
   }
   return(mp);
}
/* procedure 21: finds min integer at ip0 in struct */
static ATL_mmnode_t *FindMinIntInMMQ
(
   ATL_mmnode_t *bp,   /* queue to be searched */
   void *ip0           /* ptr to integer withinin node bp */
)
/*
 * RETURNS: ptr to structure containing min int value at byte offset
 *          offset in struct
 */
{
   ATL_mmnode_t *mp=NULL, *p;
   int *ip;
   int val;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   if (!bp)
      return(NULL);

   mp = bp;
   ip = (int*)(((char*)bp) + offset);
   val = *ip;
   for (p=bp->next; p; p = p->next)
   {
      ip = (int*)(((char*)p) + offset);
      if (*ip < val)
      {
         mp = p;
         val = *ip;
      }
   }
   return(mp);
}

/* procedure 22: finds first integer equal to val at ip0 in struct */
static ATL_mmnode_t *FindIntValInMMQ
(
   ATL_mmnode_t *bp,   /* queue to be searched */
   void *ip0,          /* ptr to integer withinin node bp */
   int val             /* value being searched for */
)
/*
 * RETURNS: ptr to first structure containing value val at byte offset
 *          offset in struct, or NULL if no such value found
 */
{
   ATL_mmnode_t *mp=NULL, *p;
   int *ip;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   if (!bp)
      return(NULL);

   for (p=bp; p; p = p->next)
   {
      ip = (int*)(((char*)p) + offset);
      if (*ip == val)
         return(p);
   }
   return(NULL);
}

/* procedure 23: sorts Q from least-to-greatest on int val at ip0 in struc */
static ATL_mmnode_t *SortMMQByIntVal
(
   ATL_mmnode_t *bp,   /* queue to be sorted */
   void *ip0           /* ptr to integer withinin node bp to sort on*/
)
/*
 * RETURNS: possibly new queue base, sorted from least-to-greatest on int at ip0
 */
{
   ATL_mmnode_t *sb=NULL, *p;
   int *ip;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   if (!bp)
      return(NULL);

   while(bp)
   {
      ip = (int*)(((char*)bp) + offset);
      p = FindMaxIntInMMQ(bp, ip);
      bp = RemoveMMNodeFromQ(bp, p);
      p->next = sb;
      sb = p;
   }
   return(sb);
}

/* procedure 24: reverses order in Q */
static ATL_mmnode_t *ReverseMMQ(ATL_mmnode_t *bp)
/*
 * RETURNS: new base ptr of reversed queue
 */
{
   ATL_mmnode_t *nb=NULL, *p;
   while(bp)
   {
      p = bp;
      bp = bp->next;
      p->next = nb;
      nb = p;
   }
   return(nb);
}

/* procedure 25: places all nodes wt int value val at ip0 in new queue */
static ATL_mmnode_t *YankMMNodesByIntVal
(
   ATL_mmnode_t **bp0,  /* queue to be searched */
   void *ip0,          /* ptr to integer withinin node *bp */
   int val             /* value to be yanked out of original Q */
)
/*
 * Finds all nodes that have the integeral value val stored in position
 * ip0-bp0 in nodes.  These nodes are removed from bp0, and placed in
 * their own queue, which is returned.  bp0 is modified in the process.
 * RETURNS: ptr to queue of nodes wt integer value val
 */
{
   ATL_mmnode_t *bp=(*bp0), *p, *valb=NULL, *vp;
   int *ip;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   while(bp)
   {
      p = FindIntValInMMQ(bp, (((char*)bp)+offset), val);  /* find node */
      if (!p)       /* if there are no more in bp, we are done */
         break;
      bp = RemoveMMNodeFromQ(bp, p);   /* remove it from original queue */
/*
 *    Add node at front of new value-only queue
 */
      if (valb)
      {
         vp->next = p;
         vp = p;
      }
      else
         vp = valb = p;
   }
   *bp0 = bp;
   return(valb);
}

/* procedure 26 */
static ATL_mmnode_t *ATL_SortMMNodesByMflop
(
   int imf,            /* which mflop entry to sort on */
   ATL_mmnode_t *bp    /* queue to be sorted */
)
/*
 * kills original queue, and returns a greatest-to-least sorted queue
 * on p->mflop[imf].  Does it with O(N^2) alg, but if this is a bottleneck,
 * we never get here because timing takes an eternity.
 */
{
   ATL_mmnode_t *p, *prev, *sb=NULL;   /* ptr, prev, sorted base */
   ATL_mmnode_t *minp, *minprev;
   double mf;

/*
 * Sort from greatest-to-least by always adding smallest entry in old
 * list to head of greatest-to-least list
 */
   while (bp)
   {
/*
 *    Find slowest remaining kernel
 */
      mf = bp->mflop[imf];
      for (minp=prev=bp, p=bp->next; p; p = p->next)
      {
         if (p->mflop[imf] < mf)
         {
            minp = p;
            mf = p->mflop[imf];
            minprev = prev;
         }
         prev = p;
      }
/*
 *    Remove it from unsorted queue, and add as new head of sorted
 */
      if (minp == bp)
      {
         bp = bp->next;
         minp->next = sb;
      }
      else   /* in the middle of unsorted queue */
      {
         minprev->next = minp->next;
         minp->next = sb;
      }
      sb = minp;
   }
   return(sb);
}

/* procedure 27 */
static ATL_mmnode_t *ParseMMLine(char *ln)
/*
 * Given a line from a mm index file (with multiple lines pasted together
 * into one line (ln), return a structure describing that line.
 */
{
   ATL_mmnode_t *p;
   char *sp;
   int itmp;
   char ch;

   p = GetMMNode();

   sp = strstr(ln, "VLEN=");
   if (sp)
      p->vlen = atoi(sp+4+1);
   else
      p->vlen = 1;
   sp = strstr(ln, "BLASK=");
   if (sp)
      p->blask = atoi(sp+5+1);
   else
      p->blask = 0;
   sp = strstr(ln, "STGC=");
   if (sp)
      p->stgC = atoi(sp+4+1);
   else
      p->stgC = 0;
   sp = strstr(ln, "STGB=");
   if (sp)
      p->stgB = atoi(sp+4+1);
   else
      p->stgB = 0;
   sp = strstr(ln, "STGA=");
   if (sp)
      p->stgA = atoi(sp+4+1);
   else
      p->stgA = 0;
   sp = strstr(ln, "SZEXTRA=");
   if (sp)
      p->szExtra = atoi(sp+7+1);
   else
      p->szExtra = 0;
   sp = strstr(ln, "PFLS=");
   if (sp)
      p->pfLS = atoi(sp+4+1);
   else
      p->pfLS = 0;
   sp = strstr(ln, "IVAR=");
   if (sp)
      p->ivar = atoi(sp+4+1);
   else
      p->ivar = 0;
   sp = strstr(ln, "PREF=");
   if (sp)
      p->pref = atoi(sp+4+1);
   else
      p->pref = 0;
   sp = strstr(ln, "KBMAX=");
   if (sp)
      p->kbmax = atoi(sp+5+1);
   else
      p->kbmax = 0;
   sp = strstr(ln, "KBMIN=");
   if (sp)
      p->kbmin = atoi(sp+5+1);
   else
      p->kbmin = 0;
   sp = strstr(ln, "KU=");
   if (sp)
      p->ku = atoi(sp+2+1);
   else
      p->ku = 0;
   sp = strstr(ln, "NU=");
   if (sp)
      p->nu = atoi(sp+2+1);
   else
      p->nu = 0;
   sp = strstr(ln, "MU=");
   if (sp)
      p->mu = atoi(sp+2+1);
   else
      p->mu = 0;
   sp = strstr(ln, "MB=");
   if (sp)
      p->mbB = atoi(sp+2+1);
   else
      p->mbB = 0;
   sp = strstr(ln, "NB=");
   if (sp)
      p->nbB = atoi(sp+2+1);
   else
      p->nbB = 0;
   sp = strstr(ln, "KB=");
   if (sp)
      p->kbB = atoi(sp+2+1);
   else
      p->kbB = 0;
   sp = strstr(ln, "SSE=");
   if (sp)
      p->SSE = atoi(sp+3+1);
   else
      p->SSE = 0;

   sp = strstr(ln, "ID=");
   if (sp)
      p->ID = atoi(sp+2+1);
   else
      p->ID = 0;

   sp = strstr(ln, "OPMV=");
   if (sp)
   {
      int imv;
      imv = atoi(sp+5);
      ATL_MMF_MVPUT(p->flag, imv);
   }

   sp = strstr(ln, "UPPER=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<MMF_UPPER);
      else
         p->flag &= ~(1<<MMF_UPPER);
   }
   sp = strstr(ln, "RIGHT=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<MMF_RIGHT);
      else
         p->flag &= ~(1<<MMF_RIGHT);
   }
   sp = strstr(ln, "COMPLEX=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_COMPLEX);
      else
         p->flag &= ~(1<<MMF_COMPLEX);
   }
   sp = strstr(ln, "FKO=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<MMF_FKO);
      else
         p->flag &= ~(1<<MMF_FKO);
   }
   sp = strstr(ln, "KCLN=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<MMF_KCLN);
      else
         p->flag &= ~(1<<MMF_KCLN);
   }
   sp = strstr(ln, "BREG1=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<MMF_BREG1);
      else
         p->flag &= ~(1<<MMF_BREG1);
   }
   sp = strstr(ln, "NOBCAST=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_NOBCAST);
      else
         p->flag &= ~(1<<MMF_NOBCAST);
   }
   sp = strstr(ln, "KUISKB=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_KUISKB);
      else
         p->flag &= ~(1<<MMF_KUISKB);
   }
   sp = strstr(ln, "KVEC=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<MMF_KVEC);
      else
         p->flag &= ~(1<<MMF_KVEC);
   }
   sp = strstr(ln, "L14NB=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<MMF_L14NB);
      else
         p->flag &= ~(1<<MMF_L14NB);
   }
   sp = strstr(ln, "PFCELTS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFCELTS);
      else
         p->flag &= ~(1<<MMF_PFCELTS);
   }
   sp = strstr(ln, "PFBCOLS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFBCOLS);
      else
         p->flag &= ~(1<<MMF_PFBCOLS);
   }
   sp = strstr(ln, "PFABLK=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_PFABLK);
      else
         p->flag &= ~(1<<MMF_PFABLK);
   }
   sp = strstr(ln, "PFACOLS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFACOLS);
      else
         p->flag &= ~(1<<MMF_PFACOLS);
   }
   sp = strstr(ln, "AOUTER=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_AOUTER);
      else
         p->flag &= ~(1<<MMF_AOUTER);
   }
   sp = strstr(ln, "KRUNTIME=");
   if (sp)
   {
      if (atoi(sp+8+1))
         p->flag |= (1<<MMF_KRUNTIME);
      else
         p->flag &= ~(1<<MMF_KRUNTIME);
   }
   sp = strstr(ln, "LDCTOP=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_LDCTOP);
      else
         p->flag &= ~(1<<MMF_LDCTOP);
   }
   sp = strstr(ln, "X87=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<MMF_X87);
      else
         p->flag &= ~(1<<MMF_X87);
   }

   sp = strstr(ln, "MFLOP=");
   if (sp)
      GetDoubleArr(sp+6, 8, p->mflop);

   sp = strstr(ln, "ASM=");
   if (sp)
      p->asmbits = asmNames2bitfield(sp+4);


   sp = strstr(ln, "TB='");
   if (sp)
   {
      ch = tolower(sp[4]);
      if (ch == 'n')
         p->TB = AtlasNoTrans;
      else if (ch == 'c')
         p->TB = AtlasConjTrans;
      else if (ch == 't')
         p->TB = AtlasTrans;
      else
         assert(0);
   }
   sp = strstr(ln, "TA='");
   if (sp)
   {
      ch = tolower(sp[4]);
      if (ch == 'n')
         p->TA = AtlasNoTrans;
      else if (ch == 'c')
         p->TA = AtlasConjTrans;
      else if (ch == 't')
         p->TA = AtlasTrans;
      else
         assert(0);
   }

   sp = strstr(ln, "CFLAGS='");
   if (sp)
      p->cflags = GetSingleQuoteString(sp+6+1);
   else
      p->cflags = NULL;

   sp = strstr(ln, "COMP='");
   if (sp)
      p->comp = GetSingleQuoteString(sp+4+1);
   else
      p->comp = NULL;

   sp = strstr(ln, "AUTH='");
   if (sp)
      p->auth = GetSingleQuoteString(sp+4+1);
   else
      p->auth = NULL;

   sp = strstr(ln, "ROUT='");
   if (sp)
      p->rout = GetSingleQuoteString(sp+4+1);
   else
      p->rout = NULL;

   return(p);
}

/* procedure 18 */
static void PrintMMLine(FILE *fpout, ATL_mmnode_t *np)
{
   int i, j, k;
   char ta, tb;

   if (!np)
      return;
   if (!np->rout)
      np->ID = 0;
   if (np->TA == AtlasConjTrans) ta = 'C';
   else if (np->TA == AtlasTrans) ta = 'T';
   else ta = 'N';
   if (np->TB == AtlasConjTrans) tb = 'C';
   else if (np->TB == AtlasTrans) tb = 'T';
   else tb = 'N';
   fprintf(fpout, "ID=%d ROUT='%s' AUTH='%s' TA='%c' TB='%c' \\\n",
           np->ID, np->rout ? np->rout : "generated",
           np->auth ? np->auth : "R. Clint Whaley", ta, tb);

   fprintf(fpout, "   ");
   i = 3;
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   i += fprintf(fpout, "OPMV=%d ", ATL_MMF_MVGET(np->flag));
   if (np->blask)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "BLASK=%d ", np->blask);
   }
   if (np->vlen)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "VLEN=%d ", np->vlen);
   }
   if (np->pfLS)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PFLS=%d ", np->pfLS);
   }
   if (np->pref)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PREF=%d ", np->pref);
   }
   if (np->kbmax)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KBMAX=%d ", np->kbmax);
   }
   if (np->kbmin)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KBMIN=%d ", np->kbmin);
   }
   if (np->ivar)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "IVAR=%d ", np->ivar);
   }
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   i += fprintf(fpout, "KU=%d ", np->ku);
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   i += fprintf(fpout, "NU=%d ", np->nu);
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   i += fprintf(fpout, "MU=%d ", np->mu);

   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   if (np->mbB != 0)
      i += fprintf(fpout, "MB=%d ", np->mbB);
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   if (np->nbB != 0)
      i += fprintf(fpout, "NB=%d ", np->nbB);
   if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
   if (np->kbB != 0)
      i += fprintf(fpout, "KB=%d ", np->kbB);
   if (FLAG_IS_SET(np->flag, MMF_UPPER))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "UPPER=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_RIGHT))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "RIGHT=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_COMPLEX))
   {
      if (i+7 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "COMPLEX=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_BREG1))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "BREG1=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_KCLN))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KCLN=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_FKO))
   {
      if (i+3 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "FKO=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_NOBCAST))
   {
      if (i+7 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "NOBCAST=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_KUISKB))
   {
      if (i+6 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KUISKB=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_KVEC))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KVEC=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_L14NB))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "L14NB=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_PFCELTS))
   {
      if (i+7 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PFCELTS=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_PFBCOLS))
   {
      if (i+7 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PFBCOLS=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_PFABLK))
   {
      if (i+6 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PFABLK=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_PFACOLS))
   {
      if (i+7 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "PFACOLS=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_AOUTER))
   {
      if (i+6 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "AOUTER=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_KRUNTIME))
   {
      if (i+8 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KRUNTIME=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_LDCTOP))
   {
      if (i+6 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "LDCTOP=1 ");
   }
   if (FLAG_IS_SET(np->flag, MMF_X87))
   {
      if (i+3 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "X87=1 ");
   }

   k = 0;  /* no need to write mflop */
   for (j=0; j < 8; j++)
   {
      if (np->mflop[j] != 0.0)
      {
         k = 1;
         break;
      }
   }
   if (k)
   {
      if (i > 3) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "MFLOP=%le", np->mflop[0]);
      for (j=8-1; j && np->mflop[j] == 0.0; j--);
      for (k=1; k <= j; k++)
         i += fprintf(fpout, ",%le", np->mflop[k]);
   }
   if (np->asmbits)
   {
      if (i > 40) { fprintf(fpout, " \\\n   "); i = 3; }
      for (j=0; !(np->asmbits & (1<<j)); j++);
      assert(j < NASMD);
      i += fprintf(fpout, "  ASM=%s", ASMNAM[j]);
      for (j++; j < NASMD; j++)
         if (np->asmbits & (1<<i))
            i += fprintf(fpout, ",%s", ASMNAM[j]);
   }
   if (np->cflags)
   {
      if (i+6+strlen(np->cflags) > 70)
         { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "  CFLAGS='%s'", np->cflags);
   }
   if (np->comp)
   {
      if (i+4+strlen(np->comp) > 70)
         { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "  COMP='%s'", np->comp);
   }
   if (i)
      fprintf(fpout, "\n");
}

/* procedure 28 */
static ATL_mmnode_t *KillMMNodesByFlag(int flag, int msk, ATL_mmnode_t *bp)
{
   ATL_mmnode_t *p=bp;
   while(p)
   {
      ATL_mmnode_t *next = p->next;
      if ((p->flag & msk) != (flag & msk))
        bp = KillMMNodeFromQ(bp, p);
      p = next;
   }
   return(bp);
}

/* procedure 29 */
static void PrintMMNodes(FILE *fpout, ATL_mmnode_t *bp)
{
   while (bp)
   {
      PrintMMLine(fpout, bp);
      bp = bp->next;
   }
}

/* procedure 30 */
static void WriteMMFile(char *file, ATL_mmnode_t *nq)
{
   FILE *fpout;

   if (!file || !strcmp(file, "stdout"))
      fpout = stdout;
   else if (!strcmp(file, "stderr"))
      fpout = stderr;
   else
   {
      fpout = fopen(file, "w");
      assert(fpout);
   }
   PrintMMNodes(fpout, nq);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}

/* procedure 31 */
static void WriteMMFileWithPath
   (char pre, char *path, char *file, ATL_mmnode_t *nq)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   WriteMMFile(ln, nq);
}

/* procedure 32 */
static ATL_mmnode_t *ReadMMFile(char *file)
/*
 * Reads in a standard ATLAS parsable MM index file, and returns a
 * list of all the kernels defined there.
 */
{
   ATL_mmnode_t *nq=NULL, *p;
   FILE *fpin;
   char *ln, *sp;
   int i, j, KeepOn, len;

   if (!file || !strcmp(file, "stdin"))
      fpin = stdin;
   else
      fpin = fopen(file, "r");
   if (!fpin)
      return(NULL);
   nq = p = GetMMNode();
   while (ln = GetJoinedLines(fpin))
   {
      if (ln[0] != '#' && ln[0] != '\0')
      {
         p->next = ParseMMLine(ln);
         p = p->next;
      }
   }
   fclose(fpin);
   return(KillMMNode(nq));
}

/* procedure 33 */
static ATL_mmnode_t *ReadMMFileWithPath
   (char pre, char *path, char *file)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   return(ReadMMFile(ln));
}


/* procedure 34 */
static ATL_mmnode_t *DelRepeatedMMKernels(ATL_mmnode_t *bp)
/*
 * Deletes any repeated IDs
 */
{
   ATL_mmnode_t *prev, *p, *np;
   int ID;

   for (p=bp; p; p = p->next)
   {
      ID = p->ID;
      prev = p;
      do
      {
         for (np=p->next; np && np->ID != ID; np = np->next)
            prev = np;
         if (np)  /* found duplicate */
            prev->next = KillMMNode(np);
      }
      while (np);
   }
   return(bp);
}

/* procedure 35 */
static ATL_mmnode_t *DelBadArchMMKernels(ATL_mmnode_t *bp)
/*
 * Weeds out kernels that require SSE/assembly that we haven't got
 */
{
   int asmb=0, die;
   ATL_mmnode_t *p, *prev;
   #ifdef ATL_GAS_ARM64
      asmb |= (1<<8);
   #endif
   #ifdef ATL_GAS_ARM
      asmb |= (1<<7);
   #endif
   #ifdef ATL_GAS_MIPS
      asmb |= (1<<6);
   #endif
   #ifdef ATL_GAS_PARISC
      asmb |= (1<<5);
   #endif
   #ifdef ATL_GAS_PPC
      asmb |= (1<<4);
   #endif
   #ifdef ATL_GAS_SPARC
      asmb |= (1<<3);
   #endif
   #ifdef ATL_GAS_x8664
      asmb |= (1<<2);
   #endif
   #ifdef ATL_GAS_x8632
      asmb |= (1<<1);
   #endif

   prev = p = bp;
   while (p)
   {
      die = (p->asmbits) ? !(asmb & p->asmbits) : 0;
      #ifndef ATL_SSE3
         if (p->SSE)
         {
            die |= (p->SSE >= 3);
            #ifndef ATL_SSE2
               die |= (p->SSE >= 2);
            #endif
            #ifndef ATL_SSE1
               die |= (p->SSE >= 1);
            #endif
         }
      #endif
      if (die)
      {
         if (p == bp)
            bp = p = KillMMNode(p);
         else
            prev->next = p = KillMMNode(p);
      }
      else
      {
         prev = p;
         p = p->next;
      }
   }
   return(bp);
}


/* procedure 36 */
static int MMKernsCompat(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels have compatible storage for A,B & C, 0 otherwise
 */
{
   if (p0->blask != p1->blask)
      return(0);
   if (p0->nu != p1->nu)
      return(0);
   if (p0->mu != p1->mu)
      return(0);
   if (p0->stgC != p1->stgC)
      return(0);
   if (p0->stgB != p1->stgB)
      return(0);
   if (p0->stgA != p1->stgA)
      return(0);
   if ((p0->flag | p1->flag)&(1<<MMF_KVEC)) /* KVEC routs have additional req */
   {
      if (p0->flag&(1<<MMF_KVEC) != p1->flag&(1<<MMF_KVEC))
         return(0);
      if (p0->vlen != p1->vlen)
         return(0);
   }
   return(1);
}

/* procedure 37 */
static int MMKernCleansK(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if p1 can provide K-cleanup for p0, else 0
 * NOTE: upper 3 bits of p->flag must indicate whether KU=4,3,2 are legeal
 *       KU for cleanup (most sig to least).  If any kbB that kern p0 is used
 *       for is not a multiple of a given KU, then it is ruled out, and so has
 *       this bit unset.
 * To do K-cleanup in our scheme, MAX(ku, kbmin) <= 4, kbmax == 0, K must
 * by a runtime variable, and it must have the same storage for all matrices.
 * We could avoid requiring same  A/B storage, but this would require using
 * different copy routs for KB0 blk, which would stress i-cache even worse,
 * and might force us to bring in more copy funcs overall.  Since mu/nu
 * present in C, this relaxation wouldn't help much anyway, so blow it off.
 * If we wanted to support this, could occasionally use an Mvec kernel to
 * clean a KVEC, or vice versa.
 */
{
   if (!p0 || !p1)
      return(0);
   if (p1->kbmax != 0)
      return(0);
   if (!FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
   if (FLAG_IS_SET(p0->flag, MMF_KVEC)) /* K-vectorized ops on vlen muls */
   {
      if (p1->kbmin && p1->kbmin != p0->vlen)
         return(0);
      if (p1->ku != p0->vlen) /* remember, vlen's must be same p0/p1 */
         return(0);
   }
   else /* M/N vectorized */
   {
      const int KU=p1->ku;
      if (p1->kbmin > KU)  /* kernel must handle all K cases, incl KU! */
         return(0);
      if (KU > 1)
      {
         if (KU > 4)
            return(0);
         if (!(p0->flag & (1<<(31+KU-4))))
            return(0);
      }
   }
   return(MMKernsCompat(p0, p1));
}

/* procedure 38 */
ATL_mmnode_t *MMGetGenCaseKClean(ATL_mmnode_t *mp)
/*
 * RETURNS: generator-supported case for K-cleanup for p0
 * genstr & filename left unset so user can specialize further
 */
{
   ATL_mmnode_t *gp;
   if (mp->stgA || mp->stgB || mp->stgC)  /* can't generate for formats */
      return(NULL);                       /* gen doesn't understand */
   gp = GetMMNode();
   gp->vlen = mp->vlen;
   gp->blask = mp->blask;
   gp->kbB = mp->kbB;
   gp->nbB = mp->nbB;
   gp->mbB = mp->mbB;
   gp->nu = mp->nu;
   gp->mu = mp->mu;
   gp->flag = mp->flag & (MMF_MVSET | (1<<MMF_KVEC));
   gp->flag |= 1<<MMF_KRUNTIME;
/*
 * May have problem if generator doesn't support VLEN!
 */
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      gp->ku = mp->vlen;
   else
      gp->ku = 1;
   return(gp);
}

/* procedure 39 */
ATL_mmnode_t *MMCompatKernPresent(ATL_mmnode_t *mb, ATL_mmnode_t *mt)
/*
 * RETURNS: ptr to mt-storage-compatible node in mb, else NULL
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
      if (MMKernsCompat(mt, mp))
         return(mp);
   return(NULL);
}

/* procedure 40 */
static int MMKernsSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same except for blocking, 0 otherwise
 */
{
   const int flg0=p0->flag&MMF_DIFFMSK, flg1=p1->flag&MMF_DIFFMSK;
/*
 * TRSM and GEMM kernels are same thing, but other blask differences mean no
 */
   if (p0->blask != p1->blask)
   {
      unsigned int OK;
      OK = (p0->blask == ATL_KGEMM && p1->blask == ATL_KTRSM) |
           (p1->blask == ATL_KGEMM && p0->blask == ATL_KTRSM);

      if (!OK)
         return(0);
   }
/*
 * Two generated kernels are the same if mu,nu,ku,VLEN,flag are the same.
 * However, if KUISKB, generated kernels with differing KBs are not same!
 * NOTE: pref & pfLS do not affect output, and so are not checked!
 * NOTE: any extension of generator functionality should extend this check!
 */
   if (p0->ID == 0 && p1->ID == 0)
   {
      if (FLAG_IS_SET(p0->flag, MMF_KUISKB) && p0->kbB != p1->kbB)
         return(0);
      return(p0->mu == p1->mu && p0->nu == p1->nu && p0->ku == p1->ku &&
             p0->vlen == p1->vlen && flg0 == flg1);
   }
/*
 * If both are user kernels, then they may be repeats.  For user kernels,
 * they are the same if both ID and flag match, else they are not.
 */
   else if (p0->ID > 0 && p1->ID > 0)
      return(p0->ID == p1->ID && flg0 == flg1);
   return(0);  /* Can't be the same if above criteria fails */
}
/* procedure 41 */
static int MMKernsPerfSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * ignoring blocking params, are kernels the same?  This version does not
 * ignore setting that effect only performance (eg, prefetch settings).
 */
{
   if (p0->ID == 0 && p1->ID == 0 &&
       (p0->pref != p1->pref || p0->pfLS != p1->pfLS))
      return(0);
   return(MMKernsSame(p0, p1));
}

/* procedure 42 */
static int MMKernCompsSame(ATL_mmnode_t *p0, ATL_mmnode_t *p1)
/*
 * RETURNS: 1 if kernels are the same including KB, else 0
 */
{
/*
 * Kernels are not the same if one has compile-time K and other runtime
 */
   if (FLAG_IS_SET(p0->flag, MMF_KRUNTIME) !=
       FLAG_IS_SET(p1->flag, MMF_KRUNTIME))
      return(0);
/*
 * Kernels not same if both compile-time K with differing KB
 */
   if (!FLAG_IS_SET(p0->flag, MMF_KRUNTIME) && p0->kbB != p1->kbB)
      return(0);
/*
 * Generated kernels aren't same if prefetch is different
 */
   if ((p0->ID | p1->ID) == 0 && (p0->pref != p1->pref || p0->pfLS != p1->pfLS))
      return(0);
   return(MMKernsSame(p0, p1));
}

/* procedure 43 */
static ATL_mmnode_t *MMKernIsPresent(ATL_mmnode_t *mmb, ATL_mmnode_t *mmp)
/*
 * RETURNS: 1 if kernel compilation matching mmp is in list mmb, 0 otherwise
 */
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
      if (MMKernsSame(mmp, mp))
         return(mp);
   return(NULL);
}

/* procedure 44 */
ATL_mmnode_t *MMGetKCleanQ(ATL_mmnode_t *db, ATL_mmnode_t *mb)
/*
 * RETURNS: new queue with cloned nodes from mb that can clean any kernel
 *          present in db
 */
{
   ATL_mmnode_t *dp, *kb=NULL;;
   for (dp=db; dp; dp = dp->next)
   {
      ATL_mmnode_t *mp;
      mp = mb;
      while ( (mp = MMCompatKernPresent(mp, dp)) )
      {
         if (!MMKernIsPresent(kb, mp))
         {
            ATL_mmnode_t *np;
            np = CloneMMNode(mp);
            np->next = kb;
            kb = np;
         }
         mp = mp->next;
      }
   }
   return(kb);
}


/* procedure 45 */
static ATL_mmnode_t *MMKernCompIsPresent(ATL_mmnode_t *mmb, ATL_mmnode_t *mmp)
/*
 * RETURNS: ptr to kernel compilation matching mmp is in list mmb, else NULL
 */
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
      if (MMKernCompsSame(mmp, mp))
         return(mp);
   return(NULL);
}

/* procedure 46 */
static ATL_mmnode_t *MMPerfKernIsPresent(ATL_mmnode_t *mmb, ATL_mmnode_t *mmp)
/*
 * RETURNS: 1 if kernel compilation matching mmp is in list mmb, 0 otherwise
 */
{
   ATL_mmnode_t *mp;
   for (mp=mmb; mp; mp = mp->next)
      if (MMKernsPerfSame(mmp, mp))
         return(mp);
   return(NULL);
}

/* procedure 47 */
static ATL_mmnode_t *AddUniquePerfMMKernsToList
   (ATL_mmnode_t *mmb, ATL_mmnode_t *newb)
/*
 * RETURNS: mmb with any kernels present in newb that aren't in mmb added to it
 * (1) changes mmb, does not change newb
 * (2) new nodes are added to BEGINNING of mmb (leaving mmb unsorted)!
 */
{
   ATL_mmnode_t *mp;

   for (mp=newb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      p = MMPerfKernIsPresent(mmb, mp);
      if (p)
         p->flag |= mp->flag & MMF_STKBTS;
      else
      {
         p = CloneMMNode(mp);
         p->next = mmb;
         mmb = p;
      }
   }
   return(mmb);
}
/* procedure 48 */
static ATL_mmnode_t *AddUniqueMMKernsToList
   (ATL_mmnode_t *mmb, ATL_mmnode_t *newb)
/*
 * RETURNS: mmb with any kernels present in newb that aren't in mmb added to it
 * (1) changes mmb, does not change newb
 * (2) new nodes are added to BEGINNING of mmb (leaving mmb unsorted)!
 */
{
   ATL_mmnode_t *mp;

   for (mp=newb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      p = MMKernIsPresent(mmb, mp);
      if (p)
         p->flag |= mp->flag & MMF_STKBTS;
      else
      {
         p = CloneMMNode(mp);
         p->next = mmb;
         mmb = p;
      }
   }
   return(mmb);
}

/* procedure 49 */
static ATL_mmnode_t *AddUniqueMMKernCompList
   (ATL_mmnode_t *mmb, ATL_mmnode_t *newb)
/*
 * RETURNS: mmb with any kerncomps present in newb that aren't in mmb added
 *          to it
 * (1) changes mmb, does not change newb
 * (2) new nodes are added to BEGINNING of mmb (leaving mmb unsorted)!
 */
{
   ATL_mmnode_t *mp;
   for (mp=newb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      p = MMKernCompIsPresent(mmb, mp);
      if (p)
         p->flag |= mp->flag & MMF_STKBTS;
      else
      {
         p = CloneMMNode(mp);
         p->next = mmb;
         mmb = p;
      }
   }
   return(mmb);
}
/* procedure 50 */
int MMNamelen(char pre, char *nm, ATL_mmnode_t *mp, int kb)
{
   int len=5;                         /* ATL_<pre> */

   len = strlen(nm);
   if (len == 4)  /* may be non-gemm kernel (eg., syrk, trmm) */
   {              /* only one of these, so it gets fixed name */
      if (nm[0] == 's' && nm[1] == 'y' && nm[2] == 'r' && nm[3] == 'k')
         return(12); /* ATL_<pre>amsyrkK */
   }
   len += 5;  /* ATL_<pre> */

   len += NumDecDigits(mp->ID);       /* ATL_<pre><nm><ID> */
   len += 1 + NumDecDigits(kb);       /* ATL_<pre><nm><ID>_<kb> */
   len += 1 + NumDecDigits(mp->vlen); /* +[m,k]<vlen>, V=[m,k]*/
   len += 1 + NumHexDigits(mp->pref); /* +p<pref> */
   len += 1 + NumDecDigits(mp->pfLS); /* +x<LS> */
   len += 1 + NumHexDigits(mp->flag); /* +_<flg> */
   len += 1 + NumDecDigits(mp->mu);   /* +_mu */
   len += 1 + NumDecDigits(mp->nu);   /* +_nu */
   len += 1 + NumDecDigits(mp->ku);   /* +_ku */
   return(len);     /* string terminator/beta/file ext not included! */
} /* ATL_<pre><nm><ID>_<kb>[m,n]<vlen>p<pf>x<LS>_<flg>_MUxNUxKU */

/*
 * A kernel name is how a kernel is called from code, and is set by
 * redefining ATL_USERMM, so we see that it can be unrelated to filename.
 * However, during install we will make all filenames the same as the
 * kernel name, with the following exceptions:
 * 1. _bX is done through recompilation, so is replaced by file ext (.c, etc.)
 * 2. KB will not match the kernel name in the (a,b) cases below:
 *    a. For ID != 0, KB is always 0 (same source file used regardless of
 *       compile flags, etc).
 *    b. For kerns that that need compile-time K, but can simply be recompiled
 *       for each differing KB, we will encode KB=1 (never a valid KB).
 *    c. For kerns taking a runtime KB, KB is 0
 *    d. For kerns that only supply an exact KB, KB is set to mp->kbB.
 *
 * During installation, we can use mp->str to temporarily store the filename,
 * while we use mp->rout to store the kernel (call) name, w/o beta suffix.
 */
/* procedure 51 */
int SprintMMName(char *name, char pre, char *nm, ATL_mmnode_t *mp, int kb)
/*
 * kb >= 0 : print kb in name, else don't
 */
{
   int i, ku=mp->ku;
   unsigned int flag = mp->flag & MMF_DIFFMSK;

   if (mp->blask == 1)
   {
      if (mp->mu == mp->nu && !FLAG_IS_SET(mp->flag, MMF_RIGHT))
         strcpy(name, "ATL_XsqsyrkK");
      else
         strcpy(name, "ATL_XumsyrkK");
      name[4] = pre;
      return(12);
   }
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
      ku = mp->kbB;
   /* ATL_<pre><nm><ID>_kb[m,n]<vlen>p<pf>x<LS>_<flg>_MUxNUxKU */
   i = sprintf(name, "ATL_%c%s%d_%d%c%dp%xx%d_%x_%dx%dx%d", pre, nm, mp->ID, kb,
               FLAG_IS_SET(mp->flag, MMF_KVEC) ? 'k':'m', mp->vlen,
               mp->pref, mp->pfLS, flag, mp->mu, mp->nu, mp->ku);
   return(i);
}

/* procedure 52 */
static char *GetMMKernName(char pre, char *nm, ATL_mmnode_t *mp)
/*
 * Get a string of form ATL_<pre><nm>suff, where suff encodes all info
 * required to differentiate between amm kernels.  Name must be suffixed
 * later with _b[n,1,0] to be fully qualified, and there is room left in
 * string to add that.
 */
{
   int i, h, kb=0;
   char *name;

   if (FLAG_IS_SET(mp->flag,MMF_KUISKB) || !FLAG_IS_SET(mp->flag,MMF_KRUNTIME))
      kb = mp->kbB;
   i = MMNamelen(pre, nm, mp, kb);
   name = malloc(i+3+1);  /* 3 for beta name, +1 is for string terminator */
   assert(name);
   h = SprintMMName(name, pre, nm, mp, kb);
   assert(h <= i);
   return(name);
}

/* procedure 53 */
char *GetMMLabelName(char pre, ATL_mmnode_t *mp)
/*
 * This func used to get a decent filename during search for printing.  A given search
 * doesn't have to worry about malloc/dealloc, can just call wt NULL when done
 */
{
   static char *nm=NULL;
   static int L=0;
   int nL;
   if (!mp)
   {
      free(nm);
      nm = NULL;
      L = 0;
      return(NULL);
   }
   if (mp->ID && mp->rout)
      return(mp->rout);
   nL = MMNamelen(pre, "amm", mp, mp->kbB);
   if (nL > L)
   {
      free(nm);
      nm = GetMMKernName(pre, "amm", mp);
   }
   else
   {
      nL = SprintMMName(nm, pre, "amm", mp, mp->kbB);
      assert(nL <= L);
   }
   return(nm);
}
/* procedure 54 */
static char *GetMMFilename(char pre, char *nm, ATL_mmnode_t *mp)
/*
 * Get a string of form ATL_<pre><nm>suff, where suff encodes all info
 * required to differentiate between amm kernels
 */
{
   int i, h, kb, flg, pf, pfLS;
   char *name;

   if (mp->ID)
#ifdef ATL_GENERATE  /* generators get real name */
      kb = 0;
#else  /* tuners get default name */
      return(DupString("ATL_tmp.c"));
#endif
   else if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
      kb = mp->kbB;
   else if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      kb = 0;
   else
      kb = 1;
   flg = mp->flag;
   pf = mp->pref;
   pfLS = mp->pfLS;
   if (mp->ID)       /* pref & flag don't matter for non-genned filenames */
      mp->pref = mp->pfLS = mp->flag = 0;
   else  /* remove from flag all non-filename vals */
      mp->flag &= ~(MMF_MVSET|(1<<MMF_KCLN)|(MMF_ALLTRANS));
   i = MMNamelen(pre, nm, mp, kb);
   name = malloc(i+2+1);           /* 2 file ext, 1 str term */
   assert(name);
   h = SprintMMName(name, pre, nm, mp, kb);
   mp->flag = flg;  /* restore flag/pref to correct vals */
   mp->pref= pf;
   mp->pfLS = pfLS;
   assert(h <= i);
   name[h] = '.';
   if (mp->ID == 0)
      name[h+1] = 'c';
   else if (FLAG_IS_SET(mp->flag, MMF_FKO))
      name[h+1] = 'B';
   else if (mp->rout) /* user supplied kern, get ext from orig name */
   {
      int k;
      k = strlen(mp->rout)-1;
      assert(mp->rout[k-1] == '.');
      name[h+1] = mp->rout[k];
   }
   else /* user-supplied w/o name, assume C */
      name[h+1] = 'c';
   name[h+2] = '\0';
   return(name);
}

/* procedure 55 */
static int MMSplitByFlagAny
(
   const int mask,      /* if any flags in mask are present, kern->KFLG */
   ATL_mmnode_t **MMB,  /* IN/OUT: original base ptr */
   ATL_mmnode_t **KFLG  /* IN/OUT: base ptr for all runtime-K kerns  */
)
/*
 * Splits original list MMB wt kerns wt any bit in mask set moved to KFLG,
 * while rest emain in MMB.  This possibly changes both lists, though no node
 * is lost.
 * RETURNS: number of nodes removed from MMB and added to KFLG
 */
{
   ATL_mmnode_t *mb=(*MMB), *kb=(*KFLG), *mp;
   int nk=0, nxtoff, bvoff;

   nxtoff = GetOffset(&mb->next, mb);
   bvoff =  GetOffset(&mb->flag, mb);
   while ( (mp = FindNodeWithMaskOR(mb, nxtoff, bvoff, mask)) )
   {
      mb = RemoveNodeFromList(mb, mp, nxtoff);
      mp->next = kb;
      kb = mp;
      nk++;
   }
   *MMB = mb;
   *KFLG = kb;
   return(nk);
}

/* procedure 56 */
static void MMApplyMoves2Flags
(
   ATL_mmnode_t *mmb,  /* kernel to set MMF_MV[A,B,C] flag bits */
   int mvBits          /* last 3 bits: MOVE_[CBA] */
)
{
   const unsigned int mvMSK = ~MMF_MVSET, mvSET = (mvBits&7)<<MMF_MVA;
   ATL_mmnode_t *mmp;
   for (mmp=mmb; mmp; mmp = mmp->next)
      mmp->flag = ((mmp->flag) & mvMSK) | mvSET;
}

/* procedure 57 */
void MMIvar2str(ATL_mmnode_t *mb)
/*
 * Sets str to point to the mmnode of the K cleaner; str is NULL if the routine
 * does not need K cleanup.  If it is self-cleaning, gets its own node address.
 * For any non-self-cleaning node needing cleanup, ivar is set to -1.
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
   {
      if (FLAG_IS_SET(mp->flag, MMF_KCLN))
      {
         if (mp->str)
            free(mp->str);
         if (mp->ivar == 0)
            mp->str = (char*) mp;
         else
         {
            const int cnt = mp->ivar - 1;
            int i;
            ATL_mmnode_t *kp;
            for (i=0, kp=mb; i < cnt && kp; kp = kp->next, i++);
            assert(kp);
            mp->str = (char*) kp;
            mp->ivar = -1;
         }
      }
      else
         mp->str = NULL;
   }
}
/* procedure 58 */
void MMstr2Ivar(ATL_mmnode_t *mb)
/*
 * Sets ivar according to p->str: if p->str=NULL, or p = p->str, ivar=0, else
 * ivar equals the index in mb of the K cleanup routine.  Note ivar actually
 * stores idx+1, so 0 means no externel K clean required.
 */
{
   ATL_mmnode_t *mp;

   for (mp=mb; mp; mp = mp->next)
   {
      mp->ivar = 0;
      if (mp->str && mp != (void*)mp->str)
      {
         int cnt;
         ATL_mmnode_t *p, *fnd=(void*)mp->str;
         for (cnt=0, p=mb; p && p != fnd; cnt++, p = p->next);
         assert(p);
         mp->ivar = cnt + 1;
      }
      mp->str = NULL;
   }
}
/* procedure 59 */
ATL_mmnode_t *MMFindKCleanStr(ATL_mmnode_t *mb)
/*
 * RETURNS: first node with str that is not NULL or equal to node address.
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
   {
      if (mp->str && mp->str != (void*)mp)
         return(mp);
   }
   return(NULL);
}

/* procedure 60 */
ATL_mmnode_t *MMFindStrMatch(ATL_mmnode_t *mb, void *st)
/*
 * RETURNS: first node with str equal to st
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
      if (mp->str == st)
         return(mp);
   return(NULL);
}

/* procedure 61 */
static void MMWinnowByOrder(ATL_mmnode_t *mb, int imf, double tol)
/*
 * Deletes any node from mb that isn't faster than the preceeding node
 */
{
   ATL_mmnode_t *mp;
   if (!mb)
      return;
   if (tol == 0.0)
      tol = 1.0;
   mp = mb;
   while (mp->next)
   {
      if (mp->mflop[imf] >= mp->next->mflop[imf]*tol)
         mp->next = KillMMNode(mp->next);
      else
         mp = mp->next;
   }
}

/* procedure 62 */
static ATL_mmnode_t *MMWinnowByKU(ATL_mmnode_t *mb, int minKB, double tol)
/*
 * Successively does following:
 *  (1) Move highest performing kernel from mmb to retained list
 *  (2) Get rid of all kernels that don't provide speedup over retained KUs
 * ASSUMES: timings in mflop[0] correct for this range
 */
{
   ATL_mmnode_t *rb=NULL;
   while (mb)
   {
      ATL_mmnode_t *mp;
      int KEEP=1;                  /* assume this kernel useful */
      mp = FindMaxMflopMMQ(mb, 0);
      if (rb)                      /* are there kerns that might do same job? */
      {
         ATL_mmnode_t *kp;
         const unsigned int ku=mb->ku;
/*
 *       See if this guy can handle a KU case better than existing kernels
 */
         for (kp=rb; kp; kp = kp->next)
         {
            const unsigned int rku=kp->ku;
            if (rku == ku)
            {
               KEEP = mp->mflop[0]*tol >= kp->mflop[0];
               if (!KEEP)
                  break;
            }
            else
            {
               const unsigned int gap = (rku > ku) ? rku-ku : ku-rku;
               if (gap <= 4) /* retained kern handles all cases */
               {
                  const unsigned int kb = ((minKB+ku-1)/ku) * ku;
                  double pen=(kb-3.0);
                  pen /= kb;
                  KEEP = mp->mflop[0]*tol >= pen*kp->mflop[0];
                  if (!KEEP)
                     break;
               }
            }
         }
      }
      mb = RemoveMMNodeFromQ(mb, mp);
      if (KEEP)
      {
         mp->next = rb;
         rb = mp;
      }
      else
         KillMMNode(mp);
   }
   return(rb);
}
/* procedure 63 */
void SetExtBFkoRout(char *rout)
{
   int len;
   len = strlen(rout);
   if (rout[len-1] != 'B' && rout[len-2] == '.' )
      rout[len-1] = 'B';
}
#endif  /* end atlas_mmparse.h guard */
