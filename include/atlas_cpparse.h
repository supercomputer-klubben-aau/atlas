#ifndef ATLAS_CPPARSE_H
   #define ATLAS_CPPARSE_H

#include "atlas_genparse.h"
#include "atlas_enum.h"
#ifndef uint
   #define uint unsigned int
#endif

#define ATL_CPMAXGEN  6  /* maximum generated code */
#define ATL_CPMAXACCMAJ  4  /* max access-major copy */
#define ATL_CPMVECAB  1
#define ATL_CPMVECC   2
#define ATL_CPKVECAB  3
#define ATL_CPKVECC   4
#define ATL_CPBLKMJAB 5
#define ATL_CPBLKMJC  6

#define CPF_TOBLK    0  /* set: copy to blk format;  unset: opposite direc */
#define CPF_CBLK     1  /* set: C format;  unset:A/B */
#define CPF_ABLK     2  /* only used in old ammgen-style */
#define CPF_MVEC     3  /* 1: vectorize on M dim */
#define CPF_NVEC     4  /* 1: vectorize on N dim */
#define CPF_TRANS    5  /* 1: transpose, else NoTrans */
#define CPF_CONJ     6  /* 1: Conjugate, else NoConj */
#define CPF_SINGLE   7  /* 1: single, else double */
#define CPF_REAL     8  /* 1: real, else complex */
#define CPF_BE1      9
#define CPF_BEN     10
#define CPF_BEX     11
#define CPF_BE0     12
#define CPF_AL1     13
#define CPF_ALN     14
#define CPF_ALX     15
#define CPF_SYRK    16
#define CPF_SYMM    17
#define CPF_TRMM    18
#define CPF_TRSM    19
#define CPF_ASM     20 /* only set during gen */
#define CPF_TMP     21 /* always unset after use */
#define CPF_MAXB    21 /* last used bit in flag */

#define CPF_DEFAULT   0  /* default is no bits set */
#define CPF_ALLALP ((1<<CPF_AL1)|(1<<CPF_ALN)|(1<<CPF_ALX))
#define CPF_ALLBET ((1<<CPF_BE0)|(1<<CPF_BE1)|(1<<CPF_BEN)|(1<<CPF_BEX))
#define CPF_ALLKERN ((1<<CPF_SYRK)|(1<<CPF_SYMM)|(1<<CPF_TRMM)|(1<<CPF_TRSM))

typedef struct CPNode ATL_cpnode_t;
struct CPNode
{
   ATL_cpnode_t *next;
   char *auth, *comp, *cflags;
   char *rout;     /* callable routine name w/o alpha/beta suffix */
   char *genstr;   /* string that will generate, or NULL for hand-tuned */
   char *exflags;
   char *str;
   double mflop[4];
   int ID;         /* ID unique in <pre>copy.idx */
   int STGID;      /* storage ID predefined vals: [1-6] */
   int mu, nu;     /* unrolling for C-M&N, A-K&M, B-K&N */
   int mb, nb;     /* M/N dims to use for timing */
   int vlen;       /* > 1: vector length used in copy generation */
   int kvec;       /* > 1: vec with size kvec along K dim */
   int rtlen;      /* # of chars/bytes in rout */
   uint flag;      /* main bitvec, using CPF_ macros */
   uint asmbits;   /* bitfield indicating which assembly(ies) is required */
   uint bv;        /* temp bitvec, used to determine coherence */
};
/* procedure 1 */
static ATL_cpnode_t *GetCPNode(void)
{
   ATL_cpnode_t *p;
   p = calloc(1, sizeof(ATL_cpnode_t));
   assert(p);
   p->flag = CPF_DEFAULT;
   return(p);
}

/* procedure 2 */
static void CopyCPNode(ATL_cpnode_t *p, ATL_cpnode_t *dup)
{
   ATL_cpnode_t *nxt=p->next;
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
   memcpy(p, dup, sizeof(ATL_cpnode_t));
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
   p->next = nxt;
}
/* procedure 3 */
static ATL_cpnode_t *CloneCPNode(ATL_cpnode_t *dup)
{
   ATL_cpnode_t *p;
   p = malloc(sizeof(ATL_cpnode_t));
   assert(p);
   memcpy(p, dup, sizeof(ATL_cpnode_t));
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
   p->next = NULL;
   return(p);
}

/* procedure 4: clones a queue of CP structs */
static ATL_cpnode_t *CloneCPQueue(ATL_cpnode_t *dupb)
{
   ATL_cpnode_t *p, *pd, *nb;
   if (!dupb)
      return(NULL);
   p = nb = CloneCPNode(dupb);
   for (pd=dupb->next; pd; pd = pd->next)
   {
      p->next = CloneCPNode(pd);
      p = p->next;
   }
   return(nb);
}

/* procedure 5: clones a queue of strided CP structs */
static ATL_cpnode_t *CloneStridedCPQueue
(
   ATL_cpnode_t *dupb,   /* queue of nodes to clone */
   int stride               /* increment between nodes to take */
)
/*
 * Creates a queue of cloned nodes from dupb; move stride each time
 * (stride must be >= 1); i.e. skip stride-1 structs in original queue
 */
{
   ATL_cpnode_t *p, *pd, *nb;
   int i;

   if (!dupb)
      return(NULL);
   if (stride == 1)
      return(CloneCPQueue(dupb));
   assert(stride > 1);
   p = nb = CloneCPNode(dupb);
   pd = nb;
   while(pd)
   {
      for (i=0; i < stride && pd; i++, pd = pd->next);
      if (pd)
      {
         p->next = CloneCPNode(pd);
         p = p->next;
      }
      else
         p->next = NULL;
   }
   return(nb);
}

/* procedure 6 */
static ATL_cpnode_t *KillCPNode(ATL_cpnode_t *die)
{
   ATL_cpnode_t *p=NULL;
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
      free(die);
   }
   return(p);
}

/* procedure 7: safely remove nukeme from Q, reseting all links */
static ATL_cpnode_t *RemoveCPNodeFromQ
(
   ATL_cpnode_t *Q,     /* queue of nodes */
   ATL_cpnode_t *nukeme /* node to remove from queue */
)
/*
 * Removes nukeme from Q, sets nukeme->next=NULL, and returns updated Q
 */
{
   ATL_cpnode_t *p, *prev;

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
static ATL_cpnode_t *SwapCPNodeInQ
(
   ATL_cpnode_t *Q,     /* queue of nodes */
   ATL_cpnode_t *p0,    /* first node to be swap */
   ATL_cpnode_t *p1     /* second node to be swap */
)
/*
 * Puts p0 in queue at place p1 was, and p1 in place p0 was.
 * RETURNS: possibly changed Q
 */
{
   ATL_cpnode_t *p, *prv0=NULL, *prv1=NULL, *nxt;
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
static ATL_cpnode_t *KillCPNodeFromQ
(
   ATL_cpnode_t *Q,     /* queue of nodes */
   ATL_cpnode_t *nukeme /* node to remove from queue */
)
{
   Q = RemoveCPNodeFromQ(Q, nukeme);
   KillCPNode(nukeme);
   return(Q);
}

/* procedure 10: frees all ->rout entries in queue b */
static void KillAllCPRouts(ATL_cpnode_t *b)
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

/* procedure 11: frees all ->genstr entries in queue b */
static void KillAllCPGenstr(ATL_cpnode_t *b)
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
static void KillAllCPNodes(ATL_cpnode_t *die)
{
   while (die)
      die = KillCPNode(die);
}

/* procedure 13 */
static void ATL_SubGoodGccInCPNodes
(
   ATL_cpnode_t *bp   /* queue to make sub in */
)
/*
 *  Gets GOODGCC (from Make.inc), and substitutes it for all comp == "gcc"
 *  in the queue.  This gets us mandatory flags like -pg,-m64,etc.
 */
{
   ATL_cpnode_t *kp;  /* queue to make sub in */
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
static void ATL_UnsubGoodGccInCPNodes
(
   ATL_cpnode_t *bp   /* queue to make reverse sub in */
)
/*
 *  Gets GOODGCC (from Make.inc); Any comp string matching that is switched
 *  back to "gcc".  This is usually necessary so that output files don't
 *  use an old GOODGCC that lacks something like -pg.
 */
{
   ATL_cpnode_t *kp;  /* queue to make sub in */
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
static void ResubGoodGccInCPNodes
(
   ATL_cpnode_t *bp   /* queue to make sub in */
)
/*
 * Takes gcc compiler that use GOODGCC, and replaces them with "gcc"
 * to help portability
 */
{
   ATL_cpnode_t *kp;  /* queue to make sub in */
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
static int ATL_CountNumberOfCPNodes
(
    ATL_cpnode_t *bp   /* queue to count */
)
{
   int i;
   for (i=0; bp; i++, bp = bp->next);
   return(i);
}

/* procedure 17 */
static ATL_cpnode_t *ATL_LastCPNode(ATL_cpnode_t *bp)
/*
 * RETURNS: pointer to last node in queue
 */
{
   ATL_cpnode_t *p;
   if (!bp)
      return(NULL);
   for (p=bp; p->next; p = p->next);
   return(p);
}

/* procedure 18: adds q1 to end of q0 */
static ATL_cpnode_t *ATL_JoinCPQs(ATL_cpnode_t *q0, ATL_cpnode_t *q1)
{
   ATL_cpnode_t *mp;
   if (!q1)
      return(q0);
   if (!q0)
      return(q1);
   mp = ATL_LastCPNode(q0);
   mp->next = q1;
   return(q0);
}

/* procedure 19: finds node with max mflop[imf]  */
static ATL_cpnode_t *FindMaxMflopCPQ
(
   ATL_cpnode_t *bp,   /* queue to be searched */
   int imf
)
/*
 * RETURNS: ptr to structure containing max value in mflop[imf]
 */
{
   ATL_cpnode_t *mp, *p;
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
static ATL_cpnode_t *FindMaxIntInCPQ
(
   ATL_cpnode_t *bp,   /* queue to be searched */
   void *ip0           /* ptr to integer withinin node bp */
)
/*
 * RETURNS: ptr to structure containing max int value at byte offset
 *          offset in struct
 */
{
   ATL_cpnode_t *mp=NULL, *p;
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
static ATL_cpnode_t *FindMinIntInCPQ
(
   ATL_cpnode_t *bp,   /* queue to be searched */
   void *ip0           /* ptr to integer withinin node bp */
)
/*
 * RETURNS: ptr to structure containing min int value at byte offset
 *          offset in struct
 */
{
   ATL_cpnode_t *mp=NULL, *p;
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
static ATL_cpnode_t *FindIntValInCPQ
(
   ATL_cpnode_t *bp,   /* queue to be searched */
   void *ip0,          /* ptr to integer withinin node bp */
   int val             /* value being searched for */
)
/*
 * RETURNS: ptr to first structure containing value val at byte offset
 *          offset in struct, or NULL if no such value found
 */
{
   ATL_cpnode_t *mp=NULL, *p;
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
static ATL_cpnode_t *SortCPQByIntVal
(
   ATL_cpnode_t *bp,   /* queue to be sorted */
   void *ip0           /* ptr to integer withinin node bp to sort on*/
)
/*
 * RETURNS: possibly new queue base, sorted from least-to-greatest on int at ip0
 */
{
   ATL_cpnode_t *sb=NULL, *p;
   int *ip;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   if (!bp)
      return(NULL);

   while(bp)
   {
      ip = (int*)(((char*)bp) + offset);
      p = FindMaxIntInCPQ(bp, ip);
      bp = RemoveCPNodeFromQ(bp, p);
      p->next = sb;
      sb = p;
   }
   return(sb);
}

/* procedure 24: reverses order in Q */
static ATL_cpnode_t *ReverseCPQ(ATL_cpnode_t *bp)
/*
 * RETURNS: new base ptr of reversed queue
 */
{
   ATL_cpnode_t *nb=NULL, *p;
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
static ATL_cpnode_t *YankCPNodesByIntVal
(
   ATL_cpnode_t **bp0,  /* queue to be searched */
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
   ATL_cpnode_t *bp=(*bp0), *p, *valb=NULL, *vp;
   int *ip;
   const int offset = (int)((char*)((char*) ip0) - ((char*)bp));

   while(bp)
   {
      p = FindIntValInCPQ(bp, (((char*)bp)+offset), val);  /* find node */
      if (!p)       /* if there are no more in bp, we are done */
         break;
      bp = RemoveCPNodeFromQ(bp, p);   /* remove it from original queue */
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
static ATL_cpnode_t *ATL_SortCPNodesByMflop
(
   int imf,            /* which mflop entry to sort on */
   ATL_cpnode_t *bp    /* queue to be sorted */
)
/*
 * kills original queue, and returns a greatest-to-least sorted queue
 * on p->mflop[imf].  Does it with O(N^2) alg, but if this is a bottleneck,
 * we never get here because timing takes an eternity.
 */
{
   ATL_cpnode_t *p, *prev, *sb=NULL;   /* ptr, prev, sorted base */
   ATL_cpnode_t *minp, *minprev;
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
static ATL_cpnode_t *ParseCPLine(char *ln)
/*
 * Given a line from a cp index file (with multiple lines pasted together
 * into one line (ln), return a structure describing that line.
 */
{
   ATL_cpnode_t *p;
   char *sp;
   int itmp;
   char ch;

   p = GetCPNode();

   sp = strstr(ln, "NU=");
   if (sp)
      p->nu = atoi(sp+2+1);
   else
      p->nu = 1;
   sp = strstr(ln, "MU=");
   if (sp)
      p->mu = atoi(sp+2+1);
   else
      p->mu = 1;
   sp = strstr(ln, "NB=");
   if (sp)
      p->nb = atoi(sp+2+1);
   else
      p->nb = 0;
   sp = strstr(ln, "MB=");
   if (sp)
      p->mb = atoi(sp+2+1);
   else
      p->mb = 0;
   sp = strstr(ln, "STGID=");
   if (sp)
      p->STGID = atoi(sp+5+1);
   else
      p->STGID = 0;
   sp = strstr(ln, "ID=");
   if (sp)
      p->ID = atoi(sp+2+1);
   else
      p->ID = 0;
   sp = strstr(ln, "KVEC=");
   if (sp)
      p->kvec = atoi(sp+4+1);
   else
      p->kvec = 0;
   sp = strstr(ln, "VLEN=");
   if (sp)
      p->vlen = atoi(sp+4+1);
   else
      p->vlen = 0;
   sp = strstr(ln, "BETA=");
   if (sp)
   {
      int cnt=0;
      sp += 4 + 1;
      GET_BETA:
      {
         switch(*sp)
         {
         case '0':
            p->flag |= (1L<<CPF_BE0);
            break;
         case '1':
            p->flag |= (1L<<CPF_BE1);
            break;
         case 'N':
            p->flag |= (1L<<CPF_BEN);
            break;
         case 'X':
            p->flag |= (1L<<CPF_BEX);
            break;
         default:
            assert(0);
         }
      }
      if (*(++sp) == ',')
      {
         sp++;
         goto GET_BETA;
      }
   }
   sp = strstr(ln, "ALPHA=");
   if (sp)
   {
      int cnt=0;
      sp += 5 + 1;
      GET_ALPHA:
      {
         switch(*sp)
         {
         case '1':
            p->flag |= (1L<<CPF_AL1);
            break;
         case 'N':
            p->flag |= (1L<<CPF_ALN);
            break;
         case 'X':
            p->flag |= (1L<<CPF_ALX);
            break;
         default:
            assert(0);
         }
      }
      if (*(++sp) == ',')
      {
         sp++;
         goto GET_ALPHA;
      }
   }
   sp = strstr(ln, "TRMM=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_TRMM);
      else
         p->flag &= ~(1<<CPF_TRMM);
   }
   sp = strstr(ln, "TRSM=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_TRSM);
      else
         p->flag &= ~(1<<CPF_TRSM);
   }
   sp = strstr(ln, "SYMM=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_SYMM);
      else
         p->flag &= ~(1<<CPF_SYMM);
   }
   sp = strstr(ln, "SYRK=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_SYRK);
      else
         p->flag &= ~(1<<CPF_SYRK);
   }
   sp = strstr(ln, "REAL=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_REAL);
      else
         p->flag &= ~(1<<CPF_REAL);
   }
   sp = strstr(ln, "SINGLE=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<CPF_SINGLE);
      else
         p->flag &= ~(1<<CPF_SINGLE);
   }
   sp = strstr(ln, "CONJ=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_CONJ);
      else
         p->flag &= ~(1<<CPF_CONJ);
   }
   sp = strstr(ln, "TRANS=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<CPF_TRANS);
      else
         p->flag &= ~(1<<CPF_TRANS);
   }
   sp = strstr(ln, "NVEC=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_NVEC);
      else
         p->flag &= ~(1<<CPF_NVEC);
   }
   sp = strstr(ln, "MVEC=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_MVEC);
      else
         p->flag &= ~(1<<CPF_MVEC);
   }
   sp = strstr(ln, "ABLK=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_ABLK);
      else
         p->flag &= ~(1<<CPF_ABLK);
   }
   sp = strstr(ln, "CBLK=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<CPF_CBLK);
      else
         p->flag &= ~(1<<CPF_CBLK);
   }
   sp = strstr(ln, "TOBLK=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<CPF_TOBLK);
      else
         p->flag &= ~(1<<CPF_TOBLK);
   }

   sp = strstr(ln, "MFLOP=");
   if (sp)
      GetDoubleArr(sp+6, 4, p->mflop);

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
static void PrintCPLine(FILE *fpout, ATL_cpnode_t *np)
{
   int i, j, k;
   char ta, tb;

   if (!np)
      return;
   if (!np->rout)
      np->ID = 0;
   if (np->rout)
      fprintf(fpout, "ID=%d STGID=%d ROUT='%s' AUTH='%s' \\\n   ", np->ID,
              np->STGID, np->rout, np->auth?np->auth:"R. Clint Whaley");
   else
      fprintf(fpout, "ID=%d STGID=%d AUTH='%s' \\\n   ",
              np->ID, np->STGID, np->auth ? np->auth : "R. Clint Whaley");
   i = 3;
   if (np->flag & CPF_ALLBET)
   {
      int b, f = (np->flag >> CPF_BE1) & 15;
      int ch[4] = {'1', 'N', 'X', '0'};
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      for (k=0; k < 4; k++, f >>= 1)
         if (f&1)
            break;
      if (f&1)
      {
         i += fprintf(fpout, "BETA=%c", ch[k]);
         for (f >>= 1,k++; k < 4; k++, f >>= 1)
            if (f&1)
               i += fprintf(fpout, ",%c", ch[k]);
         i += fprintf(fpout, " ");
      }
   }
   if (np->flag & CPF_ALLALP)
   {
      int b, f = (np->flag >> CPF_AL1) & 15;
      int ch[4] = {'1', 'N', 'X', '0'};
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      for (k=0; k < 3; k++, f >>= 1)
         if (f&1)
            break;
      if (f&1)
      {
         i += fprintf(fpout, "ALPHA=%c", ch[k]);
         for (f >>= 1,k++; k < 3; k++, f >>= 1)
            if (f&1)
               i += fprintf(fpout, ",%c", ch[k]);
         i += fprintf(fpout, " ");
      }
   }
   if (np->nb)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "NB=%d ", np->nb);
   }
   if (np->mb)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "MB=%d ", np->mb);
   }
   if (np->kvec)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "KVEC=%d ", np->kvec);
   }
   if (np->vlen)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "VLEN=%d ", np->vlen);
   }
   if (np->nu)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "NU=%d ", np->nu);
   }
   if (np->mu)
   {
      if (i > 70) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "MU=%d ", np->mu);
   }
/*
 * Special code for printing NULLs in copy list
 */
   if ((np->vlen|np->mu|np->nu|np->flag|np->kvec|np->mb|np->nb)==0)
      i += fprintf(fpout, "MU=0 NU=0 ");
   if (FLAG_IS_SET(np->flag, CPF_REAL))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "REAL=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_SINGLE))
   {
      if (i+6 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "SINGLE=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_CONJ))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "CONJ=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_TRANS))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "TRANS=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_NVEC))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "NVEC=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_MVEC))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "MVEC=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_ABLK))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "ABLK=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_CBLK))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "CBLK=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_TOBLK))
   {
      if (i+5 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "TOBLK=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_TRSM))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "TRSM=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_TRMM))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "TRMM=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_SYMM))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "SYMM=1 ");
   }
   if (FLAG_IS_SET(np->flag, CPF_SYRK))
   {
      if (i+4 > 72) { fprintf(fpout, " \\\n   "); i = 3; }
      i += fprintf(fpout, "SYRK=1 ");
   }

   k = 0;  /* no need to write mflop */
   for (j=0; j < 4; j++)
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
      for (j=4-1; j && np->mflop[j] == 0.0; j--);
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
static ATL_cpnode_t *KillCPNodesByFlag(int flag, int msk, ATL_cpnode_t *bp)
{
   ATL_cpnode_t *p=bp;
   while(p)
   {
      ATL_cpnode_t *next = p->next;
      if ((p->flag & msk) != (flag & msk))
        bp = KillCPNodeFromQ(bp, p);
      p = next;
   }
   return(bp);
}

/* procedure 29 */
static void PrintCPNodes(FILE *fpout, ATL_cpnode_t *bp)
{
   while (bp)
   {
      PrintCPLine(fpout, bp);
      bp = bp->next;
   }
}

/* procedure 30 */
static void WriteCPFile(char *file, ATL_cpnode_t *nq)
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
   PrintCPNodes(fpout, nq);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}

/* procedure 31 */
static void WriteCPFileWithPath
   (char pre, char *path, char *file, ATL_cpnode_t *nq)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   WriteCPFile(ln, nq);
}

/* procedure 32 */
static ATL_cpnode_t *ReadCPFile(char *file)
/*
 * Reads in a standard ATLAS parsable CP index file, and returns a
 * list of all the kernels defined there.
 */
{
   ATL_cpnode_t *nq=NULL, *p;
   FILE *fpin;
   char *ln, *sp;
   int i, j, KeepOn, len;

   if (!file || !strcmp(file, "stdin"))
      fpin = stdin;
   else
      fpin = fopen(file, "r");
   if (!fpin)
      return(NULL);
   nq = p = GetCPNode();
   while (ln = GetJoinedLines(fpin))
   {
      if (ln[0] != '#' && ln[0] != '\0')
      {
         p->next = ParseCPLine(ln);
         p = p->next;
      }
   }
   fclose(fpin);
   return(KillCPNode(nq));
}

/* procedure 33 */
static ATL_cpnode_t *ReadCPFileWithPath
   (char pre, char *path, char *file)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   return(ReadCPFile(ln));
}


/* procedure 34 */
static ATL_cpnode_t *DelRepeatedCPKernels(ATL_cpnode_t *bp)
/*
 * Deletes any repeated IDs
 */
{
   ATL_cpnode_t *prev, *p, *np;
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
            prev->next = KillCPNode(np);
      }
      while (np);
   }
   return(bp);
}

/* procedure 35 */
static ATL_cpnode_t *DelBadArchCPKernels(ATL_cpnode_t *bp)
/*
 * Weeds out kernels that require SSE/assembly that we haven't got
 */
{
   int asmb=0, die;
   ATL_cpnode_t *p, *prev;
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
      if (die)
      {
         if (p == bp)
            bp = p = KillCPNode(p);
         else
            prev->next = p = KillCPNode(p);
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
static char CopyGetPre(int flag)
{
   char pre;
   if (flag & (1<<CPF_REAL))
      pre = (flag & (1<<CPF_SINGLE)) ? 's':'d';
   else
      pre = (flag & (1<<CPF_SINGLE)) ? 'c':'z';
   return(pre);
}

/* procedure 37 */
static char *CopyGetCompType(int flag)
{
   char *ctyp;
   if (flag & (1<<CPF_REAL))
      ctyp = (flag & (1<<CPF_SINGLE)) ? "SREAL":"DREAL";
   else
      ctyp = (flag & (1<<CPF_SINGLE)) ? "SCPLX":"DCPLX";
   return(ctyp);
}

/* procedure 38 */
static char CopyGetUpr(int flag)
{
   return((flag & (1<<CPF_SINGLE)) ? 's':'d');
}

/* procedure 39 */
static char CopyGetTrans(int flag)
{
   if (flag & (1<<CPF_TRANS))
      return((flag & (1<<CPF_CONJ))?'H':'T');
   return((flag & (1<<CPF_CONJ))?'C':'N');
}

/* procedure 39 */
static char *CopyGetDirect(int flag)
/*
 * In names, direction is in reference to col-major, not block!
 */
{
   return((flag&(1<<CPF_TOBLK)) ? "From" : "Into");
}

/* procedure 40 */
static int CopyEncodeScal
(
   unsigned int flag,  /* flag encoding at least matrix type */
   int ialp,           /* for mat='A', [1,-1,X], for mat='C': [0,1,-1,X] */
   int ibet            /* for mat='A', ignored, else [0,1,-1,X] */
)
{
   flag &= ~(CPF_ALLALP | CPF_ALLBET); /* clear alpha,beta bits */
   if (ialp == 1)
      flag |= (1L<<CPF_AL1);
   else if (ialp == -1)
      flag |= (1L<<CPF_ALN);
   else if (ialp != -2)
      flag |= (1L<<CPF_ALX);
   if (flag & (1L<<CPF_CBLK))
   {
      if (ibet == 1)
         flag |= (1L<<CPF_BE1);
      else if (ibet == 0)
         flag |= (1L<<CPF_BE0);
      else if (ibet == -1)
         flag |= (1L<<CPF_BEN);
      else if (ibet != -2)
         flag |= (1L<<CPF_BEX);
   }
   return(flag);
}

/* procedure 41 */
void CopyEncodeAllScal
(
   ATL_cpnode_t *cb,
   unsigned int flag,  /* flag encoding at least matrix type */
   int ialp,           /* for mat='A', [0,1,-1,2,-2] */
   int ibet            /* for mat='A', ignored, else [0,1,-1,2,-2] */
)
{
   ATL_cpnode_t *cp;

   for (cp=cb; cp; cp = cp->next)
      cp->flag = CopyEncodeScal(cp->flag, ialp, ibet);
}

/* procedure 42 */
static int CopyEncode
(              /* all args take 0, which means don't set */
   char pre,  /* [s,d,c,z] */
   char dir,  /* 'I': copy Into block, 'F': copy From block */
   char mat,  /* [C,A] */
   char TA    /* [N,T,C,H] */
)
/*
 * returns arginfo as ATL_cpnode_t flag
 */
{
   int flag;
   if (pre)
   {
      flag = (pre == 'd' || pre == 's') ? (1<<CPF_REAL) : 0;
      flag |= (pre == 's' || pre == 'c') ? (1<<CPF_SINGLE) : 0;
   }
   if (dir)
   {
      if ((dir == 'I' || dir == 'i'))
         flag |= (1<<CPF_TOBLK);
   }
   if (mat)
   {
      if (mat == 'C' || mat == 'c')
         flag |= (1<<CPF_CBLK);
      else /* A or B */
      {
         flag |= (mat == 'a' || mat == 'A') ? (1<<CPF_ABLK) : 0;
         if (TA)
         {
            if (TA == 't' || TA == 'T')
               flag |= (1<<CPF_TRANS);
            else if (TA == 'c' || TA == 'C')
               flag |= (1<<CPF_CONJ);
            else if (TA == 'H' || TA == 'h')
               flag |= (1<<CPF_TRANS) | (1<<CPF_CONJ);
         }
      }
   }
   return(flag);
}

/* procedure 43 */
static int CopyGetBetaI(int flag)
{
   int iret;
   if (flag&(1<<CPF_BE1))
      iret = 1;
   else if (flag&(1<<CPF_BEN))
      iret = -1;
   else if (flag&(1<<CPF_BE0))
      iret = 0;
   else
      iret = 2;
   return(iret);
}

/* procedure 44 */
static char CopyGetBetaC(int flag)
{
   char ch;
   if (flag&(1<<CPF_BE1))
      ch = '1';
   else if (flag&(1<<CPF_BEN))
      ch = 'N';
   else if (flag&(1<<CPF_BE0))
      ch = '0';
   else if (flag&(1<<CPF_BEX))
      ch = 'X';
   else
      assert(0);
   return(ch);
}
/* procedure 45 */
static int CopyGetAlphaI(int flag)
{
   int iret;
   if (flag&(1<<CPF_AL1))
      iret = 1;
   else if (flag&(1<<CPF_ALN))
      iret = -1;
   else
      iret = 2;
   return(iret);
}

/* procedure 46 */
static char CopyGetAlphaC(int flag)
{
   char ch;
   if (flag&(1<<CPF_AL1))
      ch = '1';
   else if (flag&(1<<CPF_ALN))
      ch = 'N';
   else if (flag&(1<<CPF_ALX))
      ch = 'X';
   else
      assert(0);
   return(ch);
}

/* procedure 47 */
char *GetCopyGenStr(ATL_cpnode_t *p)
{
   char *gs=NULL, *vecd = p->kvec ? "kdim" : "no";
   const int flag=p->flag;
   int sz, ialp, ibet=0;
   char pre;

   pre = CopyGetPre(flag);
   ialp = CopyGetAlphaI(flag);
   if (flag & (1<<CPF_CBLK))
   {
      char *gen;
      ibet = CopyGetBetaI(flag);
      if (p->vlen > 1)
      {
         assert((flag&(1<<CPF_SYRK)) == 0);
         gen = (flag&(1<<CPF_TOBLK)) ? "C2blk":"blk2C";
         if (flag&(1<<CPF_SINGLE))
            vecd = (p->vlen == 8) ? "avx" : "sse";
         else
            vecd = (p->vlen == 4) ? "avx" : "sse";
         sz = 45 + NumDecDigits(p->mu);
         sz += NumDecDigits(p->nu);
         sz += NumDecDigits(ialp);
         sz += NumDecDigits(ibet);
         sz += strlen(p->rout);
         gs = malloc(sz);
         assert(gs);
         sprintf(gs, "make %s_%s mu=%d nu=%d alpha=%d beta=%d rt=%s",
                 gen, vecd, p->mu, p->nu, ialp, ibet, p->rout);
         return(gs);
      }
      if (flag & (1<<CPF_TOBLK))
         gen = (flag&(1<<CPF_SYRK)) ? "syC2blk":"C2blk";
      else
         gen = (flag&(1<<CPF_SYRK)) ? "syblk2C":"blk2C";
      sz = 74 + NumDecDigits(p->mu);
      sz += NumDecDigits(p->nu);
      sz += NumDecDigits(p->kvec);
      sz += strlen(p->rout);
      sz += strlen(gen) + strlen(vecd);
      gs = malloc(sz);
      assert(gs);
      sprintf(gs,
"make gen_%s pre=%c vlen=%d mu=%d nu=%d cpvlen=1 alpha=%d beta=%d vec=%s rt=%s",
              gen, pre, p->kvec, p->mu, p->nu, ialp, ibet, vecd, p->rout);
      return(gs);
   }
   else
   {
      int sz, i;
      const char *frm;
      sz = 17 + 24;
      if (flag&(1<<CPF_TRMM)) /* make target gen_tAN2blk */
         sz++;
      sz += strlen(p->rout);
      sz += NumDecDigits(p->kvec);
      sz += NumDecDigits(p->nu);
      sz += NumDecDigits(p->mu);
      sz += NumDecDigits(ialp)+1;
      gs = malloc(sz*sizeof(char));
      assert(gs);
      if (flag&(1<<CPF_REAL))
      {
         if (flag&(1<<CPF_TRMM))
            frm = (flag&(1<<CPF_TOBLK)) ? "make gen_tA%c2blk":"make gen_blk2A%c";
         else
            frm = (flag&(1<<CPF_TOBLK)) ? "make gen_A%c2blk":"make gen_blk2A%c";
      }
      else /* complex */
      {
         if (flag&(1<<CPF_TRMM))
            frm = (flag&(1<<CPF_TOBLK)) ? "make gen_ctA%c2blk":"make gen_cblk2A%c";
         else
            frm = (flag&(1<<CPF_TOBLK)) ? "make gen_cA%c2blk":"make gen_cblk2A%c";
      }
      i = sprintf(gs, frm, (flag&(1<<CPF_TRANS))?'T':'N');
      i += sprintf(gs+i, " rt=%s kmaj=%u UR=%u ku=%u alpha=%d", p->rout,
                   p->kvec, p->nu, p->mu, ialp);
      assert(i < sz);
   }
   return(gs);
}

/* procedure 48 */
void SetAllUnsetCopyRout(ATL_cpnode_t *cb, char *rt)
{
   ATL_cpnode_t *cp;
   for (cp=cb; cp; cp = cp->next)
      if (!cp->rout)
         cp->rout = DupString(rt);
}

/* procedure 49 */
void SetAllCopyGenStr(ATL_cpnode_t *cb)
/*
 * For ID=0, ", and genstr
 */
{
   ATL_cpnode_t *cp;
   for (cp=cb; cp; cp = cp->next)
   {
      if (!cp->ID)
      {
         assert(cp->rout);
         if (cp->genstr)
            free(cp->genstr);
         cp->genstr = GetCopyGenStr(cp);
      }
   }
}

/* procedure 50 */
char *CopyFlag2Str(unsigned int flag)
/*
 * Encodes TOBLK,CBLK,TRANS,CONJ,KERN,alp,bet into standard string
 */
{
   static char nm[16];
   if (flag&(1<<CPF_TOBLK))
   {
      nm[0] = 'F'; nm[1] = 'r'; nm[2] = 'o'; nm[3] = 'm';
   }
   else
   {
      nm[0] = 'I'; nm[1] = 'n'; nm[2] = 't'; nm[3] = 'o';
   }
   nm[4] = (flag&(1<<CPF_CBLK)) ? 'C' : 'A';
   if (flag&(1<<CPF_TRANS))
      nm[5] = (flag&(1<<CPF_CONJ)) ? 'H' : 'T';
   else
      nm[5] = (flag&(1<<CPF_CONJ)) ? 'C' : 'N';
   if (!(flag&CPF_ALLKERN))
      nm[6] = 'g';
   else if (flag&(1<<CPF_SYRK))
      nm[6] = 'r';
   else if (flag&(1<<CPF_SYMM))
      nm[6] = 'y';
   else if (flag&(1<<CPF_TRMM))
      nm[6] = 't';
   else if (flag&(1<<CPF_TRSM))
      nm[6] = 's';
   else
      assert(0);

   nm[7] = '_';
   nm[8] = 'a';
   if (flag & (1<<CPF_AL1))
      nm[9] = '1';
   else
      nm[9] = (flag & (1<<CPF_ALN)) ? 'N' : 'X';
   if (flag&(1<<CPF_CBLK))
   {
      nm[10] = 'b';
      if (flag & (1<<CPF_BEN))
         nm[11] = 'N';
      else if (flag & (1<<CPF_BEX))
         nm[11] = 'X';
      else
         nm[11] = (flag & (1<<CPF_BE1)) ? '1' : '0';
      nm[12] = '\0';
   }
   else
      nm[10] = '\0';
   return(nm);
}

/* procedure 51 */
char *ConjCopyName(ATL_cpnode_t *cp)
{
   if (!(cp->flag&(1<<CPF_CBLK)))
   {
      char *sp = cp->rout + 12;
      char ch;
      assert(cp->rout);
      ch = *sp;
      assert(ch == 'T' || ch == 'N');
      *sp = (ch == 'T') ? 'H' : 'C';
   }
}

char *UnConjCopyName(ATL_cpnode_t *cp)
{
   if (!(cp->flag&(1<<CPF_CBLK)))
   {
      char *sp = cp->rout + 12;
      char ch;
      assert(cp->rout);
      ch = *sp;
      assert(ch == 'H' || ch == 'C');
      *sp = (ch == 'H') ? 'T' : 'N';
   }
}

/* procedure 52 */
char *GetCopyName(ATL_cpnode_t *p, int exlen)
{
   char *sp=NULL, *dir, *flgstr;
   const int ID=p->ID, flag=p->flag;
   int len, i=0;
   char pre, cal;

   pre = CopyGetPre(flag);
   dir = CopyGetDirect(flag);
   flgstr = CopyFlag2Str(flag);
   len = strlen(flgstr);
/*
 * For C format routine naming scheme is:
 *    ATL_<pre>cp[Into,From]C[N,T][g,k,y,r,s]_aXbX_<mu>x<nu>x<blksz>
 */
   cal = CopyGetAlphaC(flag);
   if (flag & (1<<CPF_CBLK))
   {
      char ck='g', cbe;
      unsigned int blksz = p->mu * p->nu;
      if (p->kvec > 1)
         blksz = ((blksz+p->kvec-1)/p->kvec)*p->kvec;
      len += 12+NumDecDigits(p->mu)+NumDecDigits(p->nu)+NumDecDigits(blksz);
      sp = malloc(len+exlen);
      assert(sp);
      i = sprintf(sp, "ATL_%ccp%s_%ux%ux%u", pre, flgstr, p->mu, p->nu, blksz);
   }
/*
 * For A format routine naming scheme is:
 *    ATL_<pre>cp[Into,From]A[N,T][g,k,y,r,s]_aX_<ku>x<nu>_<kvec>
 */
   else
   {
      len += 12 + NumDecDigits(p->nu) + NumDecDigits(p->mu) +
             NumDecDigits(p->kvec);
      sp = malloc(len+exlen);
      assert(sp);
      i = sprintf(sp, "ATL_%ccp%s_%ux%u_%u", pre, flgstr, p->mu, p->nu,p->kvec);
   }
   assert(i < len);
   return(sp);
}

/* procedure 53 */
static int CopyAreDiff(ATL_cpnode_t *c0, ATL_cpnode_t *c1)
/*
 * RETURNS: 1 if c0 & c1 have different functionality, 0 if equivalent
 */
{
   if (c0->STGID != c1->STGID)
      return(1);
   if (c0->STGID)  /* for now, each STGID matches with any other same ID */
      return(0);   /* this may change when we actually support */
   if (c0->nu != c1->nu)
      return(1);
   if (c0->mu != c1->mu)
      return(1);
   if (c0->flag & (1<<CPF_CBLK))
   {
      int msk = ~((1<<CPF_ASM)|(1<<CPF_TMP));
      if ((c0->flag&msk) != (c1->flag&msk))
         return(1);
      if (c0->kvec || c1->kvec)
      {
         int bsz0=c0->mu*c0->nu, bsz1=c1->mu*c1->nu;
         if (c0->kvec)
            bsz0 = ((bsz0+c0->kvec-1)/c0->kvec)*c0->kvec;
         if (c1->kvec)
            bsz1 = ((bsz1+c1->kvec-1)/c1->kvec)*c1->kvec;
         if (bsz0 != bsz1)
            return(1);
      }
   }
   else
   {
      int msk = ~((1<<CPF_ABLK)|(1<<CPF_ASM)|(1<<CPF_TMP)|CPF_ALLBET);
      if ((c0->flag&msk) != (c1->flag&msk))
         return(1);
      if (c0->kvec != c1->kvec)
         return(1);
   }
   return(0);
}

/* procedure 54 */
ATL_cpnode_t *FindEquivCopy(ATL_cpnode_t *cb, ATL_cpnode_t *dup)
/*
 * RETURNS: NULL if no copy rout equivalent to dup is in cb, else ptr to dup
 */
{
   ATL_cpnode_t *cp;
   if (dup)
   {
      for (cp=cb; cp; cp = cp->next)
         if (!CopyAreDiff(cp, dup))
            return(cp);
   }
   return(NULL);
}

/* procedure 55 */
ATL_cpnode_t *FindEquivUserCopy(ATL_cpnode_t *cb, ATL_cpnode_t *dup)
/*
 * Searches list of user-supplied kerns for one that matches dup.
 * dup is a fully-qualified case, where CBLK,ABLK,TRANS are set to the exact
 * case required, while cb is coming from a user cases index file.
 * Therefore cb TRANS setting in reference to B matrix. I.e., TRANS=1
 * means NoTrans A, Trans B, while TRANS=0 means Trans A, NoTrans B.
 * RETURNS: NULL if no copy rout equivalent to dup is in cb, else ptr to dup
 */
{
   ATL_cpnode_t *cp;
   if (dup)
   {
      const int dflag=dup->flag, dalp=dflag&CPF_ALLALP, dbet=dflag&CPF_ALLBET;
      const int CBLK=dup->flag&(1<<CPF_CBLK), dab=dflag&(CPF_ALLBET|CPF_ALLALP);
      const int TRANS=(dflag>>CPF_TRANS)&1, ABLK=(dflag>>CPF_ABLK)&1;
      for (cp=cb; cp; cp = cp->next)
      {
         int flg0=cp->flag, ab;
         if (CBLK != (flg0&(1<<CPF_CBLK)))
             continue;
         if (!(flg0 & dalp))
            continue;
         if (CBLK && !(flg0 & dbet))
            continue;
         if (!CBLK)
         {
            const int TA0 = (flg0>>CPF_TRANS)&1;
            if (ABLK)
            {
               if (TA0 == TRANS)
                  continue;
            }
            else if (TA0 != TRANS)
               continue;
         }
         cp->flag = (cp->flag & (~(CPF_ALLBET|CPF_ALLALP))) | dab;
         if (!CopyAreDiff(cp, dup))
         {
            cp->flag = flg0;
            return(cp);
         }
         cp->flag = flg0;
      }
   }
   return(NULL);
}

/* procedure 56 */
static char *GetCpySumNam(int flag, char ch)
{
   char pre, calp;
   static char fnam[16];

   pre = CopyGetPre(flag);
   calp = CopyGetAlphaC(flag);
   if (flag & (1<<CPF_CBLK))
   {
      char c0='b', c1='C', cbet;
      cbet = CopyGetBetaC(flag);
      if (flag & (1<<CPF_TOBLK))
      { c0 = 'C'; c1 = 'b'; }
      sprintf(fnam, "%c%c%c2%c_a%cb%c.CPS", pre, ch, c0, c1, calp, cbet);
   }
   else
   {
      char c0='b', c1='A';
      if (flag & (1<<CPF_TOBLK))
      { c0 = 'A'; c1 = 'b'; }
      sprintf(fnam, "%c%c%c2%c_a%c.CPS", pre, ch, c0, c1, calp);
   }
   return(fnam);
}

/* procedure 57 */
void CopyApplyBlasRules(ATL_cpnode_t *cb)
{
   ATL_cpnode_t *cp;
/*
 * SYRK A/B copies are same as GEMM A/B copies.
 * TRMM is treated totally like GEMM; don't need one of A or B for TRMM,
 * but may want to use same gemm to avoid extra copy.
 */
   for (cp=cb; cp; cp = cp->next)
   {
      if ((cp->flag & (1<<CPF_SYRK)) && !(cp->flag & (1<<CPF_CBLK)))
         cp->flag ^= 1<<CPF_SYRK;
      if (cp->flag & (1<<CPF_TRMM))
         cp->flag ^= 1<<CPF_TRMM;
   }
}

/* procedure 58 */
ATL_cpnode_t *GetCopyMatches(ATL_cpnode_t *cb, int flag)
/*
 * RETURNS: cloned queue of all nodes in cb that match flag on:
 *    ALLKERN,TOBLK,CBLK,TRANS,CONJ,ALx,BEx(only for C).
 * NOTE: cb is not changed.
 */
{
   ATL_cpnode_t *cp, *cn=NULL;
   const unsigned int exactmsk=CPF_ALLKERN|(1<<CPF_CBLK)|(1<<CPF_TOBLK)
      |(1<<CPF_TRANS)|(1<<CPF_CONJ);
   for (cp=cb; cp; cp = cp->next)
   {
      ATL_cpnode_t *tp;
      if ((flag&exactmsk) != (cp->flag&exactmsk))
         continue;
      if (!(flag & (cp->flag&CPF_ALLALP)))
         continue;
      if (flag&(1<<CPF_CBLK))
         if (!(flag & (cp->flag&CPF_ALLBET)))
            continue;
      tp = CloneCPNode(cp);
      tp->next = cn;
      cn = tp;
   }
   return(cn);
}

/* procedure 59 */
ATL_cpnode_t *FindCopy_cohere(ATL_cpnode_t *cb, ATL_cpnode_t *dup)
/*
 * Searches cb for kern that match dup, except for flag, where the only
 * bits that matter are: CBLK,MVEC,NVEC,KERN
 * RETURNS: NULL if no copy rout equivalent to dup is in cb, else ptr to dup
 */
{
   const unsigned int flgD=dup->flag;
   const unsigned msk=((1<<CPF_CBLK)|(1<<CPF_MVEC)|(1<<CPF_NVEC)|CPF_ALLKERN);
   ATL_cpnode_t *cp;
   dup->flag &= msk;
   for (cp=cb; cp; cp = cp->next)
   {
      const unsigned int flg=cp->flag;
      cp->flag &= msk;
      if (!CopyAreDiff(cp, dup))
      {
         dup->flag = flgD;
         cp->flag = flg;
         return(cp);
      }
      cp->flag = flg;
   }
   dup->flag = flgD;
   return(NULL);
}

/* procedure 60 */
void CopyMarkDup_cohere(ATL_cpnode_t *b0, ATL_cpnode_t *b1)
/*
 * Set CPF_TMP for all kerns of b0 & b1 that appear in both
 */
{
   ATL_cpnode_t *p0;
   if (!b0 || !b1)
      return;
   for (p0=b0; p0; p0 = p0->next)
   {
      ATL_cpnode_t *p1;
      p1 = FindCopy_cohere(b1, p0);
      if (p1)
      {
         p0->flag |= (1<<CPF_TMP);
         p1->flag |= (1<<CPF_TMP);
      }
   }
}
/* procedure 61 */
ATL_cpnode_t *CopySortMarkedFirst(ATL_cpnode_t *cb)
/*
 * Ensures nodes with CPF_TMP set appear first in cb.
 */
{
   ATL_cpnode_t *cp=cb;
   while (cp)
   {
      ATL_cpnode_t *cs;
      while (cp && (cp->flag & (1<<CPF_TMP)))  /* skip set nodes */
         cp = cp->next;
/*
 *    If we find unset node, find next set node, and swap them
 */
      if (cp)
      {
         for (cs=cp->next; cs && !(cp->flag & (1<<CPF_TMP)); cs = cs->next);
         if (cs)
         {
            cb = SwapCPNodeInQ(cb, cp, cs);
            cp = cs->next;
         }
         else
            break;
      }
   }
   return(cb);
}

/* procedure 62 */
ATL_cpnode_t *CPfindMaskedFlag(ATL_cpnode_t *cb, const uint msk, uint bv)
/*
 * RETURNS: first node in cb whose (msk&->flag) == (bv&msk)
 */
{
   if (cb)
   {
      ATL_cpnode_t *mp;
      bv &= msk;
      for (mp=cb; mp && (mp->flag&msk)!=bv; mp = mp->next);
      return(mp);
   }
   return(NULL);
}

/* procedure 63 */
ATL_cpnode_t *CPfindMaskedBV(ATL_cpnode_t *cb, const uint msk, uint bv)
/*
 * RETURNS: first node in cb whose (msk&->bv) == (bv&msk)
 */
{
   if (cb)
   {
      ATL_cpnode_t *mp;
      bv &= msk;
      for (mp=cb; mp && ((mp->bv&msk) != bv); mp = mp->next);
      return(mp);
   }
   return(NULL);
}

/* procedure 64 */
ATL_cpnode_t *CPmoveMaskBVNewQ(const uint msk, const uint bv, ATL_cpnode_t **CB)
/*
 * Assumes set bits in msk are a superset of mask bits set in bv.
 * Removes any entries in cb that match exactly on msk with bv and returns
 * a new queue with those peeled off nodes.
 */
{
   ATL_cpnode_t *cb = *CB, *nb=NULL;
   if (cb)
   {
      ATL_cpnode_t *cp=cb;
      while(cp)
      {
         ATL_cpnode_t *nxt=cp->next;
         if (((cp->bv)&msk) == bv)
         {
            cb = RemoveCPNodeFromQ(cb, cp);
            cp->next = nb;
            nb = cp;
         }
         cp = nxt;
      }
   }
   *CB = cb;
   return(nb);
}

/* procedure 65 */
ATL_cpnode_t *CopyDictateOrder(ATL_cpnode_t *b0, ATL_cpnode_t *b1)
/*
 * Queue b1 is reordered to match corresponding nodes in b0.  A node
 * corresponds as dictated by FindCopy_cohere.
 * ASSUMES: len(b0) <= len(b1); b0 & b1 are arrs of same type excpt for scalars.
 * RETURNS: possibly new base of b1
 */
{
   ATL_cpnode_t *p0;
   unsigned int cnt0=0;
   for (p0=b0; p0; p0 = p0->next, cnt0++)
   {
      ATL_cpnode_t *p1;
      p1 = FindCopy_cohere(b1, p0);
/*
 *    I must swap p1 with cnt0 in b1
 */
      if (p1)
      {
         ATL_cpnode_t *p;
         unsigned int cnt1=0, cnt2=0;
/*
 *       Find node presently at cnt0 in b1, and swap it with p1
 */
         for (p=b1; cnt1 != cnt0 && p; cnt1++, p = p->next);
         assert(p);
         if (p != p1)
            b1 = SwapCPNodeInQ(b1, p1, p);
      }
   }
   return(b1);
}

/* procedure 66 */
ATL_cpnode_t *CopyNoRep(ATL_cpnode_t *cb, int minSz)
/*
 * Eliminates cpnodes with duplicate copy routines. In throwing away duplicates,
 * retain routine with mb*nb > minSz, if both > minSz, throw away biggest,
 * if both < minSz throw away smallest.  Also eliminates NULL cpy kerns,
 * which are signaled by mu & nu = 0
 */
{
   ATL_cpnode_t *cp, *dup, *prev;
/*
 * Eliminate all NULL cpy kerns in first pass thru kerns
 */
   while (!(cb->mu | cb->nu))
      cb = KillCPNode(cb);
   if (!cb)
      return(NULL);
   cp = cb->next;
   prev = cb;
   while (cp)
   {
      ATL_cpnode_t *nxt=cp->next;
      if (!(cp->mu | cp->nu))
         prev->next = KillCPNode(cp);
      else
         prev = cp;
      cp = nxt;
   }
/*
 * Now, eliminate duplicates from queue containing no NULLs
 */
   for (cp=cb; cp && cp->next; cp = cp->next)
   {
      unsigned int sz0 = cp->mb * cp->nb;
      ATL_cpnode_t *dup;
      while ( (dup = FindEquivCopy(cp->next, cp)) )
      {
         unsigned int sz1 = dup->mb * dup->nb;;
         if (sz0 > minSz && sz1 > minSz) /* elim largest */
         {
            if (sz0 > sz1)
            {
               sz0 = sz1;
               cp->mb = dup->mb;
               cp->nb = dup->nb;
            }
         }
         else /* elim smallest */
         {
            if (sz0 < sz1)
            {
               sz0 = sz1;
               cp->mb = dup->mb;
               cp->nb = dup->nb;
            }
         }
         KillCPNodeFromQ(cb, dup);
      }
   }
   return(cb);
}

/* procedure 67 */
void CopyFixTransByMtx(char mtx, char ta, ATL_cpnode_t *cb)
/*
 * This function takes a list of copy functions and applies the Tranpose
 * rules necessary for A or B storage.  In NoTrans the K dim is contiguous,
 * while in Transpose K is strided.  Therefore, A (naturally MxK) must always
 * swap the transpose setting, while B (KxN) keeps the input setting.
 * More explicitly, the rules are:
 *    AN:cpT, AC:cpH, AT:cpN, AH:cpC;
 *    BN:cpN, BC:cpC, BT:cpyT, BH:cpyH
 */
{
   ATL_cpnode_t *cp;
   const int nmsk=~((1<<CPF_TRANS)|(1<<CPF_CONJ)|(1<<CPF_ABLK));
   int msk;
   if (mtx == 'A' || mtx == 'a')
   {
      msk = 1<<CPF_ABLK;
      if (ta == 'N')
         msk |= (1<<CPF_TRANS);
      else if (ta == 'C')
         msk |= (1<<CPF_TRANS) | (1<<CPF_CONJ);
      else if (ta == 'H')
         msk |= (1<<CPF_CONJ);
   }
   else
   {
      if (ta == 'T')
         msk = (1<<CPF_TRANS);
      else if (ta == 'C')
         msk = (1<<CPF_CONJ);
      else if (ta == 'H')
         msk = (1<<CPF_TRANS) | (1<<CPF_CONJ);
      else /* no-trans */
         msk = 0;
   }
   for (cp=cb; cp; cp = cp->next)
   {
      cp->flag = (cp->flag & nmsk) | msk;
      if (cp->rout)
         free(cp->rout);
      cp->rout = GetCopyName(cp, 0);
   }
}

#endif  /* end atlas_cpparse.h guard */
