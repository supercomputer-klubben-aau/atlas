#ifndef ATLAS_MMTESTTIME_H
   #define ATLAS_MMTESTTIME_H

#include "atlas_mmparse.h"
#include "atlas_gentesttime.h"

#define ATL_WANT_ILCM 1
#include "atlas_iopt.h"
#ifndef Mmax
   #define Mmax(x_, y_) ( ((x_) >= (y_)) ? (x_):(y_) )
#endif
/* procedure 1 */
double MMMflop2Time(char pre, int M, int N, int K, double mf)
/*
 * Translates mf back into time, assuming gemm timing of M,N,K
 */
{
   return((((2.0*M)*N)*K) / (1000000.0*mf));
}

/* procedure 2 */
double MMTime2Mflop(char pre, int M, int N, int K, double time)
/*
 * Translates mf back into time, assuming gemm timing of M,N,K
 */
{
   return((((2.0*M)*N)*K) / (1000000.0*time));
}

/* procedure 3 */
static double CopyMflop2Time(char pre, int M, int N, double mf)
/*
 * Translates mflop rate mf into time, assuming M*N flops
 */
{
   return(((2.0*M)*N) / (1e6*mf));
}

/* procedure 4 */
char *MMGetGenName(char pre, int kb, ATL_mmnode_t *mp)
{
#ifdef ATL_GENERATE  /* generating funcs get unique name */
   char *nm;
   switch(mp->blask)
   {
   case ATL_KSYRK:
      nm = "syrk";
      break;
   case ATL_KSYMM:
      nm = "symm";
      break;
   case ATL_KTRMM:
      nm = "trmm";
      break;
   case ATL_KTRSM:
      nm = "trsm";
      break;
   case ATL_KGECP2A:
      nm = "cp2A";
      break;
   case ATL_KGECPFA:
      nm = "cpFA";
      break;
   case ATL_KGECP2C:
      nm = "cp2C";
      break;
   case ATL_KGECPFC:
      nm = "cpFC";
      break;
   case ATL_KSKCP2C:
      nm = "cp2SK";
      break;
   case ATL_KSKCPFC:
      nm = "cpFSK";
      break;
   case ATL_KGEMM:
   default:
      nm = "amm";
   }
   return(GetMMFilename(pre, nm, mp));
#else                       /* tuning gets default name that is overwritten */
   return(DupString("ATL_tmp.c")); /* avoids cluttering dir wt unused files */
#endif
}

/* procedure 5 */
char *MMGetCpGenString
(
   char pre,                    /* precision/type prefix */
   ATL_mmnode_t *mp,            /* mmkern ptr */
   int ialp,                    /* 0,-1,1,2 */
   int ibet                     /* 0,-1,1,2 */
)
{
   char *sp, *gen, *vec;
   size_t sz;

   if (mp->vlen < 2)
      vec = "no";
   else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      vec = "kdim";
   else
      vec = "mdim";

   if (mp->blask == ATL_KGECPFC || mp->blask == ATL_KGECP2C ||
       mp->blask == ATL_KSKCPFC || mp->blask == ATL_KSKCP2C)
   {
      switch (mp->blask)
      {
      case ATL_KGECPFC:
         gen = "C2blk";
         break;
      case ATL_KGECP2C:
         gen = "blk2C";
         break;
      case ATL_KSKCPFC:
         gen = "syC2blk";
         break;
      case ATL_KSKCP2C:
         gen = "syblk2C";
         break;
      }
      sz = 72 + NumDecDigits(mp->mu);
      sz += NumDecDigits(mp->nu);
      sz += NumDecDigits(mp->vlen);
      sz += strlen(mp->rout);
      sp = malloc(sz);
      assert(sp);
      sprintf(sp,
"make gen_%s pre=%c vl=%d mu=%d nu=%d cpvl=1 alpha=%d beta=%d vec=%s rt=%s",
              gen, pre, mp->vlen, mp->mu, mp->nu, ialp, ibet, vec, mp->rout);
      return(sp);
   }
   sz = 53 + NumDecDigits(Mmax(mp->mu, mp->nu));
   sz += strlen(mp->rout);
   sp = malloc(sz);
   assert(sp);
   sprintf(sp, "make gen_A%c2blk UR=%d alpha=1 kmaj=%d vec=%s rt=%s",
           mp->TA==AtlasTrans?'T':'N', Mmax(mp->mu, mp->nu),
           FLAG_IS_SET(mp->flag,MMF_KVEC), vec, mp->rout);
   return(sp);
}

/* procedure 6 */
char *MMGetFkoGenString
(
   char pre,                    /* precision/type prefix */
   ATL_mmnode_t *mp             /* mmkern ptr */
)
{
   char *frm="make gen_%s_fko pre=%c rt=%s vec=%s vlen=%d mu=%d nu=%d ku=%d"
             " KB=%d bcast=%d pf=%d pfLS=%d";
   char *vec="mdim";
   char *gs, *sp;
   int len;
   int ku=mp->ku;
   int kb, bc;

   if (mp->blask >= ATL_KGECPFA && mp->blask <= ATL_KSKCP2C)
      return(MMGetCpGenString(pre, mp, 1, 1));
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         ku = ((mp->kbB+mp->vlen-1)/mp->vlen)*mp->vlen;
      else
         ku = mp->kbB;
   }
   len = strlen(frm) + strlen(mp->rout) + 4;
   gs = malloc(len);
   assert(gs);
   if (mp->vlen < 2)
      vec = "no";
   else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
   {
      vec = "kdim";
      if (!FLAG_IS_SET(mp->flag, MMF_KUISKB))
         assert(mp->ku % mp->vlen == 0);
   }
   else
      assert(mp->mu % mp->vlen == 0);
   switch(mp->blask)
   {/* 0:ammm, 1:syrk, 2:symm, 3:trmm */
   case ATL_KSYRK:
      sp = "amsyrk";
      break;
   case ATL_KSYMM:
      sp = "amsymm";
      break;
   case ATL_KTRMM:
      sp = "amtrmm";
      break;
   case ATL_KTRSM:
      sp = "amtrsm";
      break;
   case ATL_KGECPFA:
   case ATL_KGECP2A:
   case ATL_KGECPFC:
   case ATL_KGECP2C:
   case ATL_KSKCPFC:
   case ATL_KSKCP2C:
      fprintf(stderr, "Copy routines not supported by FKO gen, kernel type=%d\n",
              mp->blask);
      exit(1);
   default:
      sp = "amm";
   }
   if (FLAG_IS_SET(mp->flag, MMF_NOBCAST))
      bc = 0;
   else
      bc = FLAG_IS_SET(mp->flag, MMF_BREG1) ? 3 : 1;

   if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      kb = 0;
   else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      kb = ((mp->kbB+mp->vlen-1)/mp->vlen)*mp->vlen;
   else
      kb = mp->kbB;

   if (FLAG_IS_SET(mp->flag , MMF_FKO))
     SetExtBFkoRout(mp->rout);
   sprintf(gs, frm, sp, pre, mp->rout, vec, mp->vlen, mp->mu, mp->nu, ku, kb,
           bc, mp->pref, mp->pfLS);
   return(gs);
}

char *GetFKOComp(int i)
{
   return("../../../iFKO/fkoc");
}
char *GetFKOFlags(int i)
{
   return("-vec");
}

/* procedure 7 */
char *MMGetGenString
(
   char pre,                    /* precision/type prefix */
   ATL_mmnode_t *mp             /* mmkern ptr */
)
{
   char *frm="make gen_%s pre=%c rt=%s vec=%s vlen=%d mu=%d nu=%d ku=%d"
             " KB=%d bcast=%d pf=%d pfLS=%d";
   char *vec="mdim";
   char *gs, *sp;
   int len;
   int ku=mp->ku;
   int kb, bc, i;

   if (FLAG_IS_SET(mp->flag, MMF_FKO))
      return(MMGetFkoGenString(pre, mp));
   if (mp->blask >= ATL_KGECPFA && mp->blask <= ATL_KSKCP2C)
      return(MMGetCpGenString(pre, mp, 1, 1));
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
   {
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
         ku = ((mp->kbB+mp->vlen-1)/mp->vlen)*mp->vlen;
      else
         ku = mp->kbB;
   }
   len = strlen(frm) + strlen(mp->rout) + 4;
   if (mp->blask == ATL_KTRMM)
      len += 7;
   gs = malloc(len);
   assert(gs);
   if (mp->vlen < 2)
      vec = "no";
   else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
   {
      vec = "kdim";
      if (!FLAG_IS_SET(mp->flag, MMF_KUISKB))
      {
         if (mp->ku % mp->vlen)
         {
            fprintf(stderr, "WTF,KVEC: KU=%d, VLEN=%d\n", mp->ku, mp->vlen);
            assert(mp->ku % mp->vlen == 0);
         }
      }
   }
   else
      assert(mp->mu % mp->vlen == 0);
   switch(mp->blask)
   {/* 0:ammm, 1:syrk, 2:symm, 3:trmm */
   case ATL_KSYRK:
      sp = "amsyrk";
      break;
   case ATL_KSYMM:
      sp = "amsymm";
      break;
   case ATL_KTRMM:
      sp = "amtrmm";
      break;
   case ATL_KTRSM:
      sp = "amtrsm";
      break;
   case ATL_KGECPFA:
      if (mp->TA == AtlasTrans)
         sp = "AT2blk";
      else if (mp->TA == AtlasConjTrans)
         sp = "cAT2blk";
      else if (mp->TA == AtlasConj)
         sp = "cAN2blk";
      else
         sp = "AN2blk";
      break;
   case ATL_KGECP2A:
      if (mp->TA == AtlasTrans)
         sp = "blk2AT";
      else if (mp->TA == AtlasConjTrans)
         sp = "cblk2AT";
      else if (mp->TA == AtlasConj)
         sp = "cblk2AN";
      else
         sp = "blk2AN";
      break;
   case ATL_KGECPFC:
      sp = "blk2C";
      break;
   case ATL_KGECP2C:
      sp = "C2blk";
      break;
   case ATL_KSKCPFC:
      sp = "SYC2blk";
      break;
   case ATL_KSKCP2C:
      sp = "SYblk2C";
      break;
   default:
      sp = "amm";
   }
   if (FLAG_IS_SET(mp->flag, MMF_NOBCAST))
      bc = 0;
   else
      bc = FLAG_IS_SET(mp->flag, MMF_BREG1) ? 3 : 1;

   if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
      kb = 0;
   else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      kb = ((mp->kbB+mp->vlen-1)/mp->vlen)*mp->vlen;
   else
      kb = mp->kbB;
   i = sprintf(gs, frm, sp, pre, mp->rout, vec, mp->vlen, mp->mu, mp->nu, ku,
               kb, bc, mp->pref, mp->pfLS);
   if (mp->blask == ATL_KTRMM)
   {
      int itr;
      if (FLAG_IS_SET(mp->flag, MMF_RIGHT))
         itr = FLAG_IS_SET(mp->flag, MMF_UPPER) ? 4 : 3;
      else /* Left */
         itr = FLAG_IS_SET(mp->flag, MMF_UPPER) ? 2 : 1;
      sprintf(gs+i, " TRMM=%1.1u", itr);
   }
   return(gs);
}

/* procedure 8 */
void MMFillInGenStrings(char pre, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;

   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';

   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->ID == 0 && !mp->genstr)
         mp->genstr = MMGetGenString(pre, mp);
   }
}

/* procedure 9 : assumed kbB for routname */
void MMFillInGenNames(char pre, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';

   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->ID == 0 && !mp->rout)
         mp->rout = MMGetGenName(pre, mp->kbB, mp);
   }
}

/* procedure 10 : fix genstrs and routs that may have gone stale */
void MMRefreshGenNodes(char pre, ATL_mmnode_t *mmb)
{
   ATL_mmnode_t *mp;
   const char upr=(pre == 'c' || pre == 's') ? 's' : 'z';

   for (mp=mmb; mp; mp = mp->next)
   {
      if (!mp->ID)
      {
         if (mp->rout)
            free(mp->rout);
         if (mp->genstr)
            free(mp->genstr);
         if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
         {
            if (FLAG_IS_SET(mp->flag, MMF_KVEC))
            {
               mp->ku = mp->kbmax = ((mp->kbB+mp->vlen-1)/mp->vlen)*mp->vlen;
               mp->kbmin = mp->kbmax - mp->vlen + 1;
            }
            else
               mp->kbmin = mp->kbmax = mp->kbB;
         }
         else if (mp->ku > 1)
            mp->kbmin = mp->ku;
         mp->rout = MMGetGenName(pre, mp->kbB, mp);
         mp->genstr = MMGetGenString(pre, mp);
      }
   }
}
/* procedure 11 : fix genstrs and routs that may have gone stale */
void WriteRefreshedMMFileWithPath
   (char pre, char *path, char *file, ATL_mmnode_t *nq)
{
   MMRefreshGenNodes(pre, nq);
   WriteMMFileWithPath(pre, path, file, nq);
}

/* procedure 12 */
ATL_mmnode_t *MMGetNodeGEN(char pre, int bv, int kb, int mu, int nu, int ku,
                           int vlen, int KVEC, int pf, int pfLS, char *rt)
/*
 * bv: 0:don't use broadcast for B, 1: use only 1 reg for B load,
 *     2:compile with FKO
 */
{
   ATL_mmnode_t *mp;

   mp = GetMMNode();
   mp->mu = mu;
   mp->nu = nu;
   mp->ku = ku;
   mp->vlen = vlen;
   mp->pref = pf;
   mp->pfLS = pfLS;
   if (bv&1)
      mp->flag |= (1<<MMF_NOBCAST);
   else if (bv&2)
      mp->flag |= (1<<MMF_BREG1);
   if (bv&262144) /*must be same as FKO_FLAG in gmmsearch */
      mp->flag |= (1<<MMF_FKO);
   if (!kb)
      mp->flag |= (1<<MMF_KRUNTIME);
   else
      mp->kbB = kb;
   mp->vlen = vlen;
   if (KVEC && vlen > 1)
      mp->flag |= (1<<MMF_KVEC);
   else if (vlen < 2)
      KVEC=1;
   if (!rt)
      rt = MMGetGenName(pre, kb, mp);
   mp->rout = rt;
   mp->genstr = MMGetGenString(pre, mp);
   return(mp);
}

/* procedure 13 */
int MMDoGenString(int verb, char *genstr)
{
   int err=0;
   if (verb > 2)
      printf("genstr='%s'\n", genstr);
   if (!genstr)
      err=1;
   else
      err = system(genstr);
   if (err)
      fprintf(stderr, "UNABLE TO GENERATE WITH COMMAND: '%s'\n",
              genstr ? genstr : "NULL");
   return(err);
}

/* procedure 14 */
char *MMGetTestString
(
   char pre,                    /* precision/type prefix */
   int mb, int nb, int kb,      /* dimensions to test */
   int beta,                    /* beta case to test */
   ATL_mmnode_t *mp            /* mmkern ptr */
)
{
   int i, k, Dd, Ud, Vd;
   int len;
   char *ln, *nm;

   if (FLAG_IS_SET(mp->flag,MMF_FKO))
   {
      if (!mp->comp)
         mp->comp = DupString(GetFKOComp(1));
      if (!mp->cflags)
         mp->cflags = DupString(GetFKOFlags(1));
   }
/*
 * Horrific code so we don't have to rely on snprintf, which some ancient
 * compiler may lack.  The len adjustments should be >= the following
 * sprintf calls (in same order, so they can be adjusted together)
 */
   switch(mp->blask)
   {
   case ATL_KSYRK:
      nm = "syrk";
      break;
   case ATL_KSYMM:
      nm = "symm";
      break;
   case ATL_KTRMM:
      nm = "trmm";
      break;
   default:   /* need to add copy kerns later */
      nm = "ammm";
   }
   len = 22 + 21;                              /* begin & end */
   len += strlen(mp->rout);
   if (mp->ID > 0)
      len += 9;                                /* "AMMCASES/" */
   len += 38 + NumDecDigits(mb) + NumDecDigits(nb) + NumDecDigits(kb);
   len += 16 + NumDecDigits(mp->mu)+NumDecDigits(mp->nu)+NumDecDigits(mp->ku);
      len += 4 + 2 + NumDecDigits(mp->vlen);
      len += 7 + 2 + NumDecDigits(mp->szExtra);
      len += 3 + 2 + NumDecDigits(mp->szC);
      len += 3 + 2 + NumDecDigits(mp->szB);
      len += 3 + 2 + NumDecDigits(mp->szA);
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      len += 5 + NumDecDigits(mp->vlen);
   if (beta != 1)
      len += (beta == -1) ? 21 : 6;
   if (mp->comp)
      len += 28 + strlen(mp->comp) + strlen(mp->cflags);

   ln = malloc(len);
   assert(ln);

   i = sprintf(ln, "make %c%stst mmrout=", pre, nm);
   if (mp->ID > 0)
      i += sprintf(ln+i, "AMMCASES/%s ", mp->rout);
   else
      i += sprintf(ln+i, "%s ", mp->rout);
   i += sprintf(ln+i, "pre=%c M=%d N=%d K=%d mb=0 nb=0 kb=%d ",
                pre, mb, nb, kb, FLAG_IS_SET(mp->flag, MMF_KRUNTIME)?0:kb);
   i += sprintf(ln+i, "mu=%d nu=%d ku=%d ", mp->mu, mp->nu, mp->ku);
    if (mp->vlen)
       i += sprintf(ln+i, "vlen=%d ", mp->vlen);
    if (mp->szExtra)
       i += sprintf(ln+i, "szExtra=%d ", mp->szExtra);
    if (mp->szC)
       i += sprintf(ln+i, "szC=%d ", mp->szC);
    if (mp->szB)
       i += sprintf(ln+i, "szB=%d ", mp->szB);
    if (mp->szA)
       i += sprintf(ln+i, "szA=%d ", mp->szA);
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      i += sprintf(ln+i, "kmaj=%d ", mp->vlen);
   if (beta != 1)
   {
      if (beta == -1)
         i += sprintf(ln+i, "beta=-1 betan=\"N1\" ");
      else if (beta == 0)
         i += sprintf(ln+i, "beta=%d ", beta);
      else
         assert(0);
   }
   if (mp->comp)
   {
      char ch = (pre == 'c' || pre == 's') ? 'S' : 'D';
      i += sprintf(ln+i, "%cMC=\"%s\" %cMCFLAGS=\"%s\" ",
                   ch, mp->comp, ch, mp->cflags);
   }
   i += sprintf(ln+i, "\n");
   assert(i < len);
   return(ln);
}

/* procedure 15 */
int MMKernelFailsTest
(
   char pre,                    /* precision/type prefix */
   int mb, int nb, int kb,      /* dimensions to test */
   int beta,                    /* beta case to test */
   ATL_mmnode_t *umm            /* mmkern ptr */
)
/*
 * RETURNS: 0 on success, non-zero on failure
 */
{
   char *ln;
   int i;
   char ch;

/*
 * If the file is generated, call generator to create it
 */
   if (umm->ID == 0 && !umm->genstr)
      umm->genstr = MMGetGenString(pre, umm);
   if (umm->genstr)
      assert(!MMDoGenString(0, umm->genstr));
   ln = MMGetTestString(pre, mb, nb, kb, beta, umm);
   i = system(ln);
   if (i)
   {
      fprintf(stderr, "%d of %s: FAILED COMMAND : %s\n",__LINE__,__FILE__,ln);
      if (umm->genstr)
         fprintf(stderr, "   genstr was = '%s'\n", umm->genstr);
   }
   free(ln);
   return(i);
}

/* procedure 16 */
int MMKernelFailsAnyBeta
(
   char pre,                    /* precision/type prefix */
   int mb, int nb, int kb,      /* dimensions to test */
   ATL_mmnode_t *umm            /* mmkern ptr */
)
/*
 * RETURNS: 0 on success, non-zero on failure
 */
{
   int i;

   for (i=-1; i < 2; i++)
      if (MMKernelFailsTest(pre, mb, nb, kb, i, umm))
        return(1);
   return(0);
}

/* procedure 17 */
static ATL_mmnode_t *DelBadMMKernels(char pre, int verb, ATL_mmnode_t *bp)
/*
 * Deletes all kernels that can't pass basic usage test
 * RETURNS: modifed bp queue wt failing kernels removed
 */
{
   ATL_mmnode_t *p, *prev;
   int die;

   if (verb > 0)
       printf("\nBEGIN BASIC MATMUL KERNEL TESTS:\n");

   prev = p = bp;
   while (p)
   {
      if (MMKernelFailsTest(pre, p->mbB, p->nbB, p->kbB, 0, p) ||
          MMKernelFailsTest(pre, p->mbB, p->nbB, p->kbB, 1, p) ||
          MMKernelFailsTest(pre, p->mbB, p->nbB, p->kbB, -1, p))
      {
         if (verb > 0)
            printf("   NUKING bad kernel %s(%d)\n", p->rout, p->ID);
         if (p == bp)
            bp = p = KillMMNode(p);
         else
            prev->next = p = KillMMNode(p);
      }
      else
      {
         if (verb > 0)
            printf("   Kernel %s(%d) passes basic tests\n", p->rout, p->ID);
         prev = p;
         p = p->next;
      }
   }
   printf("DONE BASIC KERNEL TESTS.\n\n");
   return(bp);
}
ATL_mmnode_t *GetWorkingUserCases(int verb, char pre)
/*
 * Reads/writes <pre>WORKING.sum, list of kernels that work on this machine.
 * RETURNS: all kernels that pass sanity check on this architecture
 */
{
   ATL_mmnode_t *mmb, *mmp;
   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';
   mmb = ReadMMFileWithPath(pre, "res", "WORKING.sum");
   if (mmb)
      return(mmb);
   mmb = ReadMMFileWithPath(pre, "AMMCASES", "amcases.idx");
   if (!mmb)
      return(mmb);
/*
 * Eliminate those kernels that can't work for any block size
 */
   for (mmp=mmb; mmp; mmp = mmp->next)
   {
      if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
         mmp->kbmin = mmp->kbmax = mmp->mbB = mmp->nbB = mmp->kbB = mmp->ku;
      else
      {
         int m;
         m = ATL_iLCM(mmp->mu, mmp->nu);
         m = ((60+m-1)/m)*m;
         mmp->mbB = mmp->nbB = mmp->kbB = m;
         if (mmp->kbmin)
            mmp->kbB = Mmax(mmp->kbB, mmp->kbmin);
         if (mmp->kbmax)
            mmp->kbB = Mmin(mmp->kbB, mmp->kbmax);
      }
   }
   mmb = DelBadMMKernels(pre, verb, mmb);
   WriteMMFileWithPath(pre, "res", "WORKING.sum", mmb);
   return(mmb);
}

/* procedure 18 */
/* 1st 2 bits have verb, next 4 have mov */
#define MMR_NOKFULL 6 /* 1st two bits have verb */
#define MMR_NOKCOMP 7
#define MMR_NOKVEC  8
#define MMR_NOMVEC  9
#define MMR_MKCOMP 10
ATL_mmnode_t *MMApplyRules(ATL_mmnode_t *ob, int flag, int maxKU, int reqVL)
/*
 * Eliminates mmkerns that don't abide by rules given in flag bits >= 6,
 * and sets 4 MMF move bits using bits 2-5.
 * flag: 1st two bits is verb, next 4 is MVS, rest selection rules (MMR_).
 */
{
   ATL_mmnode_t *nb=NULL;
   const unsigned int MVS = (((flag>>2)<<MMF_MVA) & MMF_MVSET);
   if ((flag & (1<<MMR_MKCOMP)) && (flag & (1<<MMR_NOKCOMP)))
         flag ^= 1<<MMR_MKCOMP; /* no MKCOMP, if KCOMP disallowed! */
   while(ob)
   {
      int KILL=0, REGEN=0, KILLFIX=0;
      if (reqVL)
         KILL = (ob->vlen != reqVL);

      if (!KILL&&(flag&(1<<MMR_NOKFULL))&&FLAG_IS_SET(ob->flag, MMF_KUISKB))
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* genned codes poss changed to runtime */
         {
            KILL = 0;
            ob->kbmin = ob->ku = FLAG_IS_SET(ob->flag, MMF_KVEC) ? ob->vlen : 4;
            ob->kbmax = 0;
            ob->flag ^= 1<<MMF_KUISKB;
            ob->flag |= 1<<MMF_KRUNTIME;
         }
      }
      if (!KILL && maxKU && ob->ku > maxKU)
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* for genned codes, see if we can reset KU */
         {
            if (FLAG_IS_SET(ob->flag, MMF_KVEC))
            {
               if (maxKU >= ob->vlen)
               {
                  KILL = 0;
                  ob->kbmin = ob->ku = (maxKU/ob->vlen)*ob->vlen;
               }
            }
            else
            {
               KILL = 0;
               ob->kbmin = ob->ku = Mmin(maxKU, 4);
            }
         }
      }
      if (!KILL&&(flag&(1<<MMR_NOKCOMP))&&!FLAG_IS_SET(ob->flag, MMF_KRUNTIME))
      {
         REGEN = KILL = 1;
         if (!ob->ID) /* genned codes poss changed to runtime */
         {
            KILL = 0;
            if (FLAG_IS_SET(ob->flag, MMF_KVEC))
            {
               if (ob->ku < ob->vlen)
                  ob->ku = ob->vlen;
               else
                  while (ob->ku % ob->vlen)
                     ob->ku--;
            }
            ob->kbmax = 0;
            ob->kbmin = ob->ku;
            ob->flag |= 1<<MMF_KRUNTIME;
         }
      }
      if ((flag&(1<<MMR_NOKVEC)) && FLAG_IS_SET(ob->flag, MMF_KVEC))
         KILL = 1;
      else if ((flag&(1<<MMR_NOMVEC)) && !FLAG_IS_SET(ob->flag, MMF_KVEC))
         KILL = 1;
      if (KILL)
         ob = KillMMNode(ob);
      else /* add to queue, and set it to use selected search type */
      {
         ATL_mmnode_t *nxt=ob->next;
         if (!ob->ID)
         {
            if (ob->genstr)
            {
               free(ob->genstr);
               ob->genstr = NULL;
            }
            if (ob->rout)
               free(ob->rout);
            ob->rout = DupString("ATL_tmp.c");
            ob->flag &= ~MMF_MVSET;
            ob->flag |= MVS;
            if (flag&(1<<MMR_MKCOMP) && !FLAG_IS_SET(ob->flag, MMF_KUISKB)
                && !ob->ID)
            { /* want to try run- & compile-time K */
               ATL_mmnode_t *np=NULL;
               if (FLAG_IS_SET(ob->flag, MMF_KVEC))
               {
                  if (ob->ku == ob->vlen)
                     np = CloneMMNode(ob);
               }
               else if (ob->ku == 1)
                  np = CloneMMNode(ob);
               if (np)
               {
                  if (FLAG_IS_SET(ob->flag, MMF_KRUNTIME))
                     np->flag ^= 1<<MMF_KRUNTIME;
                  else
                  {
                     np->kbmax = 0;
                     np->kbmin = np->ku;
                     np->flag |= 1<<MMF_KRUNTIME;
                  }
                  np->next = nb;
                  nb = np;
               }
            }
         }     /* end if this is generated kern */
         ob->next = nb;
         nb = ob;
         ob = nxt;
      }     /* end else for keeper case */
   }        /* end while(ob); */
   return(nb);
}


/* procedure 19 */
void MMFixGenK(char pre, ATL_mmnode_t *mb, int K)
{
   ATL_mmnode_t *mp;

   for (mp=mb; mp; mp = mp->next)
   {
      if (!mp->ID && FLAG_IS_SET(mp->flag, MMF_KUISKB) && mp->kbB != K)
      {
         mp->kbmax = mp->kbmin = mp->ku = mp->kbB = K;
         if (!mp->rout)
            mp->rout = DupString("ATL_tmp.c");
         if (mp->genstr)
            free(mp->genstr);
         mp->genstr = MMGetGenString(pre, mp);
      }
   }
}
/* procedure 20 */
double TimeMMKernel3F
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flag,                   /* 1: ignore any prior output file */
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush,                  /* >=0: size of cache flush, else ignored */
   int aFL,                     /* >=0: size of cache flush, else ignored */
   int bFL,                     /* >=0: size of cache flush, else ignored */
   int cFL                      /* >=0: size of cache flush, else ignored */
)
/*
 * flag - take actions the following actions if bit location is set:
 *    0 : ignore any prior output file on output
 *    1 : time in serial rather than parallel
 * NOTE: this kernel not updated to time copy codes.  AFAIK, doesn't need to
 */
{
   char fnam[128], ln[4096];  /* security from 1991 */
   const char *LO = FLAG_IS_SET(mmp->flag, MMF_AOUTER) ? "IJK": "JIK";
   const int vl=mmp->vlen;
   const int ku=mmp->ku;
   const int KB = (!FLAG_IS_SET(mmp->flag,MMF_KVEC)) ? kb : ((kb+vl-1)/vl)*vl;
   const int FORCETIME = flag&1, SERIAL=flag&2;
   int KU=mmp->ku;
   int DOTIME=1;
   int MV=3;  /* bit pattern on move CBA (C=4, B=2, A=1) */
   char *be, *gs0=mmp->genstr;
   int i, j;
   char ch;
   double *dp;
   MV = ((mmp->flag) >> MMF_MVA)&7;
/*
 * If it's a emit_mm generated file with a missing or bad genstring, make it
 */

   if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
      KU = KB;
   else if (ku > KB)
      KU = KB;
   if (mmp->ID == 0 && (KB != mmp->kbB || KU != mmp->ku || !mmp->genstr))
   {
      int kb0=mmp->kbB, ku0=mmp->ku;
      int blask = mmp->blask;
      mmp->kbB = KB;
      mmp->ku = KU;
      if (blask == ATL_KTRSM)
         mmp->blask = ATL_KGEMM;
      mmp->genstr = MMGetGenString(pre, mmp);
      mmp->blask = blask;
      mmp->kbB = kb0;
      mmp->ku = ku0;
   }
/*
 * If the file is generated, call generator to create it
 */
   if (mmp->genstr)
      assert(!MMDoGenString(verb, mmp->genstr));
   else if (verb > 2)
      printf("NO genstr\n");
/*
 * Figure out the name of the output file
 */
   if (FORCETIME)
      strcpy(fnam, "res/tmpout.ktim");
/*
 * PREammID_MBxNBxKB_MUxNUxKU_FLAG_v[M,K]VLENbBETA_CFLUSH
 */
   else
   {
      if (mmp->vlen < 2)
         ch = 'S';
      else
         ch = (FLAG_IS_SET(mmp->flag,MMF_KVEC)) ? 'K' : 'M';

      sprintf(fnam,
      "res/%cammm%d_%dx%dx%d_%dx%dx%d_%x_v%c%db%d_pf%xx%d_%d.ktim",
              pre, mmp->ID, mb, nb, KB, mmp->mu, mmp->nu, KU,
              mmp->flag, ch, mmp->vlen, beta, mmp->pref, mmp->pfLS, cflush);
      if (mmp->blask == ATL_KSYRK)
      {
         fnam[5] = 's';
         fnam[6] = 'y';
         fnam[7] = 'r';
         fnam[8] = 'k';
      }
      else if (mmp->blask == ATL_KSYMM)
      {
         fnam[5] = 's';
         fnam[6] = 'y';
         fnam[7] = 'm';
         fnam[8] = 'm';
      }
      else if (mmp->blask == ATL_KSYRK)
      {
         fnam[5] = 's';
         fnam[6] = 'y';
         fnam[7] = 'r';
         fnam[8] = 'k';
      }
      else if (mmp->blask == ATL_KTRMM)
      {
         fnam[5] = 't';
         fnam[6] = 'r';
         fnam[7] = 'm';
         fnam[8] = 'm';
      }
   }
/*
 * If we actually need to do timing, must also construct timer call
 */
   if (!FORCETIME)
      DOTIME = !FileExists(fnam);
   #define RESACT 2   /* want average frm ReadResultsFile for parallel */
   if (DOTIME)
   {
      char *nm="amm";
      if (mmp->blask == ATL_KSYRK)
         nm = "syk";
      else if (mmp->blask == ATL_KSYMM)
         nm = "sym";
      else if (mmp->blask == ATL_KTRSM)
         nm = "trm";
      i = sprintf(ln, "make x%c%stime_pt3f mb=%d nb=%d kb=%d",
                  pre, nm, mb, nb, KB);
      if (mflop)
         i += sprintf(ln+i, " FMFS=\"-Rf %d\"", mflop);
      if (SERIAL)
         i += sprintf(ln+i, " NPROC=1");
      if (mmp->genstr)
         i += sprintf(ln+i, " mmrout=%s", mmp->rout);
      else
         i += sprintf(ln+i, " mmrout=AMMCASES/%s", mmp->rout);
      i += sprintf(ln+i, " mu=%d nu=%d ku=%d", mmp->mu, mmp->nu, KU);
    if (mmp->vlen)
       i += sprintf(ln+i, " vlen=%d", mmp->vlen);
    if (mmp->szExtra)
       i += sprintf(ln+i, " szExtra=%d", mmp->szExtra);
    if (mmp->szC)
       i += sprintf(ln+i, " szC=%d", mmp->szC);
    if (mmp->szB)
       i += sprintf(ln+i, " szB=%d", mmp->szB);
    if (mmp->szA)
       i += sprintf(ln+i, " szA=%d", mmp->szA);
      i += sprintf(ln+i, " mvA=%d mvB=%d mvC=%d", ((mmp->flag >> MMF_MVA)&1),
                   ((mmp->flag >> MMF_MVB)&1), ((mmp->flag >> MMF_MVC)&1));
      i += sprintf(ln+i, " kmoves=\"");
      if (FLAG_IS_SET(mmp->flag, MMF_MVA))
         i += sprintf(ln+i, " -DATL_MOVEA");
      if (FLAG_IS_SET(mmp->flag, MMF_MVB))
         i += sprintf(ln+i, " -DATL_MOVEB");
      if (FLAG_IS_SET(mmp->flag, MMF_MVC))
      i += sprintf(ln+i, " -DATL_MOVEC");
      i += sprintf(ln+i, "\"");
      /*if (cflush > 0 || aFL > 0 || bFL > 0 || cFL > 0) */
      {
         if (cflush > 0)
               i+= sprintf(ln+i, " CFLUSH=\"%d ", cflush);
            else
               i+= sprintf(ln+i, " CFLUSH=\"%d", -1);
            if (aFL > 0)
               i+= sprintf(ln+i, " -Fa %d ", aFL);
            if (bFL > 0)
               i+= sprintf(ln+i, " -Fb %d ", bFL);
            if (cFL > 0)
               i+= sprintf(ln+i, " -Fc %d ", cFL);
            i+= sprintf(ln+i, " \"");
      }
      if (mmp->pref)
      {
         if (mmp->pfLS)
            i += sprintf(ln+i, " pf=%d pfLS=%d", mmp->pref, mmp->pfLS);
         else
            i += sprintf(ln+i, " pf=%d", mmp->pref);
      }
      if (beta == 1 || beta == 0)
         i += sprintf(ln+i, " beta=%d", beta);
      else
         i += sprintf(ln+i, " beta=-1 betan=\"N1\"");
      ch = (pre == 'c' || pre == 's') ? 'S' : 'D';
      if (mmp->comp)
         i += sprintf(ln+i, " %cMC=\"%s\"", ch, mmp->comp);
      if (mmp->cflags)
         i += sprintf(ln+i, " %cMCFLAGS=\"%s\"", ch, mmp->cflags);
      i += sprintf(ln+i, " outF=\"-f %s\"", fnam);
   }
   if (FORCETIME || !FileExists(fnam))
   {
      i += sprintf(ln+i, "\n");
      if (verb > 1)
         fprintf(stdout, "SYSTEM: %s", ln);
      if (i=system(ln))
      {
         fprintf(stderr, "ERROR %d IN COMMAND: %s", i, ln);
         fprintf(stderr, "   PROPOSED FILENAME: %s\n", fnam);
         if (mmp->genstr)
            fprintf(stderr, "   GENSTR='%s'\n", mmp->genstr);
         sprintf(ln, "rm -f %s\n", fnam);
         assert(!system(ln));
         exit(-1);
      }
   }
   if (mmp->genstr != gs0)  /* put genstr back to original value */
   {
      free(mmp->genstr);
      mmp->genstr = gs0;
   }
   dp = ReadResultsFile(RESACT, 0, fnam);
   if (!dp)
   {
      fprintf(stderr, "\nEmpty file '%s'!\n", fnam);
      fprintf(stderr, "From command: '%s'\n", ln);
      fprintf(stderr, "DOTIME=%d, genstr='%s'\n", DOTIME,
              (mmp->genstr) ? mmp->genstr : "");
      exit(-1);
   }
   if (mmp->genstr && DOTIME)
   {
      sprintf(ln, "rm %s\n", mmp->rout);
      i = system(ln);  /* return value unused, just to shut gcc up */
   }
   if (kb != KB)
   {
      double mf;
      mf = *((double*)ReadResultsFile(RESACT, 0, fnam));
      mf = (mf / KB)*kb;
      return(mf);
   }
   return(*((double*)ReadResultsFile(RESACT, 0, fnam)));
   #undef RESACT
}  /* end TimeMMKernel3F */

/* procedure 21 */
double TimeMMKernel
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flag,                    /* 1: ignore any prior output file */
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * flag - take actions the following actions if specified bit position is set:
 *    0 : ignore any prior output file on output
 *    1 : time in serial rather than parallel
 *    2 : time Right case rather than Left
 *    3 : do C^T = B^T * A^T, rather than C=A*B
 */
{
   char fnam[128], ln[4096];  /* security from 1991 */
   const char *LO = FLAG_IS_SET(mmp->flag, MMF_AOUTER) ? "IJK": "JIK";
   const int vl=mmp->vlen;
   const int mu=mmp->mu, nu=mmp->nu, ku=mmp->ku;
   int MB, NB, KB = (!FLAG_IS_SET(mmp->flag,MMF_KVEC)) ? kb : ((kb+vl-1)/vl)*vl;
   const int FORCETIME = flag&1, SERIAL=flag&2;
   const int PUTC=FLAG_IS_SET(mmp->flag, MMF_CPC);
   const int RIGHT=(flag&4), TRANS=(flag&8), UPPER=(flag&16), TRANSA=(flag&32);
   const int NMFLAG = (mmp->blask == ATL_KTRSM) ? (flag&(~3)) : 0;
   int KU=mmp->ku;
   int DOTIME=1;
   char *be, *gs0=mmp->genstr;
   int i, j;
   char ch;
   double *dp;
   double mfmul;
/*
 * This routine can time any kernel legal from usage rules, so mb,nb,kb may
 * not be multiples of the unrolling factors, so choose nearest MB,NB,KB,
 * and then compute the percent of useless flops. mfmul can just be multiplied
 * by the return value to compute the useful flop rate.
 */
   if (mmp->blask == ATL_KTRSM || mmp->blask == ATL_KTRMM)
   {
      int uk, un;
      const int TR = (RIGHT) ? !TRANS : TRANS;
      if (!TR)
      {
         uk = mu;
         un = nu;
      }
      else
      {
         uk = nu;
         un = mu;
      }
      MB = 0;
      KB = ((kb+uk-1)/uk)*uk;
      NB = ((nb+un-1)/un)*un;
      mfmul = ((double)kb)*kb*nb;
      mfmul /= ((double)KB)*KB*NB;
   }
   else
   {
      MB = ((mb+mu-1)/mu)*mu;
      NB = ((nb+nu-1)/nu)*nu;
      KB = ((kb+ku-1)/ku)*ku;
      if (FLAG_IS_SET(mmp->flag,MMF_KVEC) && vl > ku)
         KB = ((kb+vl-1)/vl)*vl;
      mfmul = kb;
      mfmul = mfmul / KB;  /* ratio of useful to useless flops in K loop */
      mfmul *= ((double)mb)*nb / (((double)MB)*NB);
      if (mmp->blask >= ATL_KGECPFA && mmp->blask <= ATL_KSKCP2C)
         mfmul = 1.0;  /* don't use mfmul for copy routines */
   }
/*
 * If it's a emit_mm generated file with a missing or bad genstring, make it
 */

   if (FLAG_IS_SET(mmp->flag,MMF_FKO))
   {
      if (!mmp->comp)
         mmp->comp = DupString(GetFKOComp(1));
      if (!mmp->cflags)
         mmp->cflags = DupString(GetFKOFlags(1));
   }
   if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
      KU = KB;
   else if (ku > KB)
      KU = KB;
   if (mmp->blask < ATL_KGECPFA && mmp->ID == 0 &&
       (KB != mmp->kbB || KU != mmp->ku || !mmp->genstr))
   {
      int kb0=mmp->kbB, ku0=mmp->ku;
      int blask = mmp->blask;
      mmp->kbB = KB;
      mmp->ku = KU;
      if (blask == ATL_KTRSM)
         mmp->blask = ATL_KGEMM;
      mmp->genstr = MMGetGenString(pre, mmp);
      mmp->blask = blask;
      mmp->kbB = kb0;
      mmp->ku = ku0;
   }
/*
 * If the file is generated, call generator to create it
 */
   if (mmp->genstr)
   {
      if (MMDoGenString(verb, mmp->genstr))
      {
         fprintf(stderr, "\nUnable to generate with: '%s'!\n\n", mmp->genstr);
         assert(0);
      }
   }
   else if (verb > 2)
      printf("NO genstr\n");
/*
 * Figure out the name of the output file
 */
   if (FORCETIME)
      strcpy(fnam, "res/tmpout.ktim");
/*
 * PREammID_MBxNBxKB_MUxNUxKU_FLAG_v[M,K]VLENbBETA_CFLUSH
 */
   else
   {
      if (mmp->vlen < 2)
         ch = 'S';
      else
         ch = (FLAG_IS_SET(mmp->flag,MMF_KVEC)) ? 'K' : 'M';

      sprintf(fnam,
      "res/%cammm%d_%dx%dx%d_%dx%dx%d_%x_v%c%db%d_pf%xx%d_f%x_%d.ktim",
              pre, mmp->ID, MB, NB, KB, mmp->mu, mmp->nu, KU,
              mmp->flag, ch, mmp->vlen, beta, mmp->pref, mmp->pfLS,
              NMFLAG, cflush);
      switch(mmp->blask)
      {
      case ATL_KSYRK:
         fnam[5] = 's';
         fnam[6] = 'y';
         fnam[7] = 'r';
         fnam[8] = 'k';
         break;
      case ATL_KSYMM:
         fnam[5] = 's';
         fnam[6] = 'y';
         fnam[7] = 'm';
         fnam[8] = 'm';
         break;
      case ATL_KTRSM:
         fnam[5] = 't';
         fnam[6] = 's';
         fnam[7] = (RIGHT)?'R' : 'L';
         fnam[8] = (TRANS)?'T' : 'N';
         break;
      case ATL_KTRMM:
         fnam[5] = 't';
         fnam[6] = 'm';
         fnam[7] = 'm';
         fnam[8] = 'm';
         break;
      case ATL_KGECPFA:
         fnam[5] = 'c';
         fnam[6] = 'p';
         fnam[7] = 'F';
         fnam[8] = 'A';
         break;
      case ATL_KGECP2A:
         fnam[5] = 'c';
         fnam[6] = 'p';
         fnam[7] = '2';
         fnam[8] = 'A';
         break;
      case ATL_KGECPFC:
         fnam[5] = 'c';
         fnam[6] = 'p';
         fnam[7] = '2';
         fnam[8] = 'C';
         break;
      case ATL_KGECP2C:
         fnam[5] = 'c';
         fnam[6] = 'p';
         fnam[7] = 'F';
         fnam[8] = 'C';
         break;
      default:;
      }
      if (PUTC)  /* need a diff name if we incl cost of copy */
      {
         assert(mmp->blask != ATL_KTRSM);
         if (fnam[5] != 'c' || fnam[6] != 'p')
         {
            fnam[7]='c';
            fnam[8]='p';
         }
      }
   }
/*
 * If we actually need to do timing, must also construct timer call
 */
   if (!FORCETIME)
      DOTIME = !FileExists(fnam);
   #define RESACT 2   /* want average frm ReadResultsFile for parallel */
   if (DOTIME)
   {
      char *nm="amm";
      if (mmp->blask == ATL_KSYRK)
         nm = "syk";
      else if (mmp->blask == ATL_KSYMM)
         nm = "sym";
      else if (mmp->blask == ATL_KTRSM)
         nm = "trs";
      else if (mmp->blask == ATL_KTRMM)
         nm = "trm";
      ch = (pre == 'c' || pre == 's') ? 'S' : 'D';
      if (mmp->blask >= ATL_KGECPFA && mmp->blask <= ATL_KSKCP2C)
      {
         if (mmp->blask < ATL_KGECPFC)
            i = sprintf(ln, "make %ccpytime mb=%d nb=%d TA=%c", pre, mb, nb,
                        mmp->TA == AtlasTrans ? 'T' : 'N');
         else
            i = sprintf(ln, "make %ccpytimeC mb=%d nb=%d TA=N", pre, mb, nb);
         if (mmp->TB == AtlasTrans) /* time row-wise access */
            i += sprintf(ln+i, " MTDX='-D 8000 8000 8000 1'");
         i += sprintf(ln+i, " kfnam=%s", mmp->rout);
         if (mmp->comp)
            i += sprintf(ln+i, " %cKC=\"%s\"", ch, mmp->comp);
         if (mmp->cflags)
            i += sprintf(ln+i, " %cKCFLAGS=\"%s\"", ch, mmp->cflags);
         if (beta == 1 || beta == 0)
            i += sprintf(ln+i, " beta=%d", beta);
         else if (beta == -1)
            i += sprintf(ln+i, " beta=-1 betan=\"N1\"");
         else
            i += sprintf(ln+i, " beta=2 betan=X");
/*
 *       These routines take alpha from kb, if it is: 0,1,-1,2
 */
         if (kb == 0 || kb == 1)
            i += sprintf(ln+i, " alpha=%d", kb);
         else if (kb == -1)
            i += sprintf(ln+i, " alpha=-1 alphan=\"N1\"");
         else if (kb == 2)
            i += sprintf(ln+i, " alpha=2 alphan=X");
      }
      else
      {
         if (mmp->blask == ATL_KTRSM)
            i = sprintf(ln,
                        "make %ctrsmcase kb=%d rb=%d sd=%c up=%c ta=%c tALL=%c",
                        pre, KB, NB, RIGHT?'R':'L', UPPER?'U':'L',
                        TRANSA?'T':'N', TRANS?'T':'N');
         else if (mmp->blask == ATL_KTRMM)
            i = sprintf(ln,
                        "make %ctrmmcase kb=%d rb=%d sd=%c up=%c ta=%c tALL=%c",
                        pre, KB, NB, RIGHT?'R':'L', UPPER?'U':'L',
                        TRANSA?'T':'N', TRANS?'T':'N');
         else
            i = sprintf(ln, "make x%c%stime_pt mb=%d nb=%d kb=%d",
                        pre, nm, MB, NB, KB);
         if (mmp->genstr)
            i += sprintf(ln+i, " mmrout=%s", mmp->rout);
         else
            i += sprintf(ln+i, " mmrout=AMMCASES/%s", mmp->rout);
         if (mmp->flag & (1<<MMF_KVEC))
            i += sprintf(ln+i, " kmaj=%u", mmp->vlen);
         if (PUTC)
         {
            if (FLAG_IS_SET(mmp->flag,MMF_ALLTRANS))
               i += sprintf(ln+i, " extdefs=\"-DALLTRANS_=1 -DPUTC=1\"");
            else
               i += sprintf(ln+i, " extdefs=\"-DPUTC=1\"");
         }
         else if (FLAG_IS_SET(mmp->flag,MMF_ALLTRANS))
         {
            i += sprintf(ln+i, " extdefs=\"-DALLTRANS_=1\"");
         }
         if (mmp->comp)
            i += sprintf(ln+i, " %cMC=\"%s\"", ch, mmp->comp);
         if (mmp->cflags)
            i += sprintf(ln+i, " %cMCFLAGS=\"%s\"", ch, mmp->cflags);
         if (beta == 1 || beta == 0)
            i += sprintf(ln+i, " beta=%d", beta);
         else
            i += sprintf(ln+i, " beta=-1 betan=\"N1\"");
      }
      if (mflop)
         i += sprintf(ln+i, " FMFS=\"-Rf %d\"", mflop);
      if (SERIAL)
         i += sprintf(ln+i, " NPROC=1");
      i += sprintf(ln+i, " mu=%d nu=%d ku=%d", mmp->mu, mmp->nu, KU);
    if (mmp->vlen)
       i += sprintf(ln+i, " vlen=%d", mmp->vlen);
    if (mmp->szExtra)
       i += sprintf(ln+i, " szExtra=%d", mmp->szExtra);
    if (mmp->szC)
       i += sprintf(ln+i, " szC=%d", mmp->szC);
    if (mmp->szB)
       i += sprintf(ln+i, " szB=%d", mmp->szB);
    if (mmp->szA)
       i += sprintf(ln+i, " szA=%d", mmp->szA);
      i += sprintf(ln+i, " mvA=%d mvB=%d mvC=%d", ((mmp->flag >> MMF_MVA)&1),
                   ((mmp->flag >> MMF_MVB)&1), ((mmp->flag >> MMF_MVC)&1));
      i += sprintf(ln+i, " kmoves=\"");
      if (FLAG_IS_SET(mmp->flag, MMF_MVA))
         i += sprintf(ln+i, " -DATL_MOVEA");
      if (FLAG_IS_SET(mmp->flag, MMF_MVB))
         i += sprintf(ln+i, " -DATL_MOVEB");
      if (FLAG_IS_SET(mmp->flag, MMF_MVC))
      i += sprintf(ln+i, " -DATL_MOVEC");
      i += sprintf(ln+i, "\"");
      if (mmp->pref)
      {
         if (mmp->pfLS)
            i += sprintf(ln+i, " pf=%d pfLS=%d", mmp->pref, mmp->pfLS);
         else
            i += sprintf(ln+i, " pf=%d", mmp->pref);
      }
      i += sprintf(ln+i, " outF=\"-f %s\"", fnam);
   }
   if (FORCETIME || !FileExists(fnam))
   {
      i += sprintf(ln+i, "\n");
      if (verb > 1)
         fprintf(stdout, "SYSTEM: %s", ln);
      if (system(ln))
      {
         fprintf(stderr, "ERROR IN COMMAND: %s", ln);
         fprintf(stderr, "   PROPOSED FILENAME: %s\n", fnam);
         if (mmp->genstr)
            fprintf(stderr, "   GENSTR='%s'\n", mmp->genstr);
         sprintf(ln, "rm -f %s\n", fnam);
         assert(!system(ln));
         exit(-1);
      }
   }
   if (mmp->genstr != gs0)  /* put genstr back to original value */
   {
      free(mmp->genstr);
      mmp->genstr = gs0;
   }
   dp = ReadResultsFile(RESACT, 0, fnam);
   if (!dp)
   {
      fprintf(stderr, "\nEmpty file '%s'!\n", fnam);
      fprintf(stderr, "From command: '%s'\n", ln);
      fprintf(stderr, "DOTIME=%d, genstr='%s'\n", DOTIME,
              (mmp->genstr) ? mmp->genstr : "");
      exit(-1);
   }
   if (mmp->genstr && DOTIME)
   {
      sprintf(ln, "rm %s\n", mmp->rout);
      i = system(ln);  /* return value unused, just to shut gcc up */
   }
   return(mfmul*(*((double*)ReadResultsFile(RESACT, 0, fnam))));
   #undef RESACT
}  /* end TimeMMKernel */

/* procedure 22, times all 3 vectorization schemes */
double MMTimeRangeK
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flag,                    /* 1: ignore any prior output file */
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb,              /* dimensions to time */
   int K0, int KN, int incK,    /* inclusive range of K's to average */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
{
   double mfA=0.0;
   int k, n=0;
   assert(KN >= K0);
   for (k=K0; k <= KN; k += incK, n++)
   {
      mfA += TimeMMKernel(verb, flag, mmp, pre, mb, nb, k, beta, mflop, cflush);
   }
   return(mfA/n);
}

/* procedure 23, times all 3 vectorization schemes */
ATL_mmnode_t *MMTimeWithGenKB
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int imf,                     /* which mflop array to write to */
   int flag,
   ATL_mmnode_t *mmp,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * RETURNS: NULL if kernel cannot support this case, else newly allocated
 *          fully-qualified mmnode describing exact case
 */
{
   int kbOK;
   char *sp;
   ATL_mmnode_t *mp;
   if (!mmp->ID) /* genned: kvec OK wt any mul of vlen, mvec OK wt any K */
      kbOK = FLAG_IS_SET(mmp->flag, MMF_KVEC) ? (kb%mmp->vlen == 0):1;
   else
   {
      kbOK = (mmp->kbmin) ? (nb >= mmp->kbmin) : 1;
      if (kbOK && mmp->kbmax)
         kbOK = nb <= mmp->kbmax;
   }
   kbOK = kbOK && (kb%mmp->ku == 0);
   if (!kbOK || ((mb/mmp->mu)*mmp->mu != mb) || ((nb/mmp->nu)*mmp->nu != nb))
      return(NULL);
   sp = mmp->genstr;
   mmp->genstr=NULL;
   mp = CloneMMNode(mmp);  /* get local copy so we can change */
   mmp->genstr = sp;
   mmp = mp;
   mmp->mbB = mb; mmp->nbB = nb; mmp->kbB = kb;
/*
 * Generated files may need to get genstr and related info corrected
 */
   if (mmp->ID == 0)
   {
      if (FLAG_IS_SET(mmp->flag, MMF_KUISKB))
          mmp->kbmax = mmp->kbmin = mmp->ku = kb;
      if (!FLAG_IS_SET(mmp->flag, MMF_KRUNTIME))
          mmp->kbB = kb;
      mmp->genstr = MMGetGenString(pre, mmp);
   }
   mmp->mflop[imf] = TimeMMKernel(verb, flag, mmp, pre, mb, nb, kb,
                                  beta, mflop, cflush);
   return(mmp);
}

/* procedure 24, times all 3 vectorization schemes */
void MMPruneMflopTol
(
   ATL_mmnode_t *mmb,           /* ptr to mmkern struct queue */
   int imf,                     /* which mflop array to write to */
   float tol                    /* >1 allow slow kerns, < 1 penalize later */
)
/*
 * First entry always retained.  Later entries must be better than
 * prior best when scaled by tol.  So, tol > 1 will retain slower kernels,
 * while tol < 1 will penalize later kernels.
 */
{
   ATL_mmnode_t *mp;
   double mfB;

   if (tol <= 0.0 || !mmb)
      return;
   mfB = mmb->mflop[imf];

   mp = mmb->next;
   while (mp)
   {
      ATL_mmnode_t *next = mp->next;
      double mf = mp->mflop[imf];
      if (mf*tol < mfB)
      {
         mmb = RemoveMMNodeFromQ(mmb, mp);
         KillMMNode(mp);
      }
      else if (mf > mfB)
         mfB = mf;
      mp = next;
   }
}

/* procedure 25, times all 3 vectorization schemes */
ATL_mmnode_t *MMBestWithGenKB
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int imf,                     /* which mflop array to write to */
   int flag,
   ATL_mmnode_t *mmb,           /* ptr to mmkern struct queue */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * RETURNS: cloned mmnode of best-performing kern from mmB queue, or
 *          NULL if none works for this size
 */
{
   ATL_mmnode_t *mmB=NULL, *mmp;
/*
 * find first kernel that works, declare it best so far
 */
   for (mmp=mmb; mmp && !mmB; mmp = mmp->next)
      mmB = MMTimeWithGenKB(verb, imf, flag, mmp, pre, mb, nb, kb, beta,
                            mflop, cflush);
   if (verb && mmB)
      printf("      B=(%d,%d,%d): %d,%s, mf=%.2f\n", mb, nb, kb, mmB->ID,
             mmB->rout, mmB->mflop[imf]);
   if (!mmp)        /* if no other cases to consider */
      return(mmB);  /* return best found (may be NULL) */
/*
 * Now try remaining kerns, and always keep best-performing
 */
   for (mmp=mmp->next; mmp; mmp = mmp->next)
   {
      ATL_mmnode_t *mp;
      mp = MMTimeWithGenKB(verb, imf, flag, mmp, pre, mb, nb, kb, beta,
                           mflop, cflush);
      if (mp)
      {
         if (verb)
            printf("      B=(%d,%d,%d): %d,%s, mf=%.2f\n", mb, nb, kb, mmB->ID,
                   mmB->rout, mmB->mflop[imf]);
         if (mp->mflop[imf] > mmB->mflop[imf])
         {
            KillMMNode(mmB);
            mmB = mp;
         }
         else
            KillMMNode(mp);
      }
   }
   return(mmB);
}

/* procedure 26, times all 3 vectorization schemes */
int MMTimeAllVecGen
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flag,                    /* 1: ignore any prior output file */
   ATL_mmnode_t *mp,            /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * Times all 3 generator-provided vectorization schemes.
 * RETURNS: index of best-performing variant (0-2)
 *
 * OVERWRITES:
 *   ->mflop[0] : time for M-vec wt bcast
 *   ->mflop[1] : 0 if nu not a multiple of vlen, else time for M-vec wt splat
 *   ->mflop[2] : time for K-vec
 *
 * NOTE: mu & nu should be the number of vec regs to use (eg., for M-vec,
 *       actual M unrolling will be mp->mu * mp->vlen).
 */
{
   char *gs0=mp->genstr, *rt0=mp->rout;
   int flg0=mp->flag, mb0=mp->mbB, nb0=mp->nbB, mu0=mp->mu;
   int kb0=mp->kbB, ku0=mp->ku, kbmax0=mp->kbmax, kbmin0=mp->kbmin;
   int iB=0;
   double mfB;

   if (!rt0)
      mp->rout = DupString("ATL_tmp.c");
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB))
      mp->kbmax = mp->kbmin = mp->ku = mp->kbB = kb;
   mp->kbB = kb;  mp->nbB = nb;  mp->mbB = mb;  /* override blking factor */
/*
 * Try M-vectorized using broadcast
 */
   mp->flag &= ~(1<<MMF_KVEC);    /* ask for M-vectorized kernel */
   mp->flag &= ~(1<<MMF_NOBCAST); /* clear bcast bit */
   mp->mu = mu0 * mp->vlen;
   mp->genstr = MMGetGenString(pre, mp);
   mfB = TimeMMKernel(verb, flag, mp, pre, mb, nb, kb, beta, mflop, cflush);
   mp->mflop[0] = mfB;
   free(mp->genstr);
/*
 * If legal, try M-vec with splat
 */
   if ((mp->nu % mp->vlen == 0) && !FLAG_IS_SET(mp->flag, MMF_FKO))
   {
      mp->flag |= (1<<MMF_NOBCAST);
      mp->genstr = MMGetGenString(pre, mp);
      mp->mflop[1] = TimeMMKernel(verb, flag, mp, pre, mb, nb, kb, beta, mflop,
                                  cflush);
      free(mp->genstr);
      if (mp->mflop[1] > mfB)
      {
         iB = 1;
         mfB = mp->mflop[1];
      }
      mp->flag &= ~(1<<MMF_NOBCAST); /* clear bcast bit */
   }
   else
      mp->mflop[1] = 0;
/*
 * Try K-vectorized
 */
   mp->mu = mu0;
   mp->flag |= (1<<MMF_KVEC);    /* ask for K-vectorized kernel */
   mp->ku = ((mp->ku+mp->vlen-1)/mp->vlen)*mp->vlen;
   mp->genstr = MMGetGenString(pre, mp);
   mp->mflop[2] = TimeMMKernel(verb, flag, mp, pre, mb, nb, kb, beta, mflop,
                               cflush);
   free(mp->genstr);
   if (mp->mflop[2] > mfB)
      iB = 2;
/*
 * Restore changed mp-> values & return
 */
   if (!rt0)
   {
      free(mp->rout);
      mp->rout = NULL;
   }
   mp->genstr = gs0;
   mp->mbB = mb0; mp->nbB = nb0; mp->kbB = kb0;
   mp->kbmin = kbmin0; mp->kbmax = kbmax0; mp->ku = ku0;
   return(iB);
}

/* procedure 27 */
int kernWorksThisCaseK(ATL_mmnode_t *mp, int kb)
/*
 * RETURNS: KB if mp can be used according to standard rules, else 0.
 */
{
   const unsigned int mu=mp->mu, nu=mp->nu, ku=mp->ku, vlen=mp->vlen;
   unsigned int KB;

   if (mp->kbmin && kb < mp->kbmin)
      return(0);
   if (mp->kbmax && kb > mp->kbmax )
      return(0);
   KB = ((kb+ku-1)/ku)*ku;
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
   {
      if (KB-kb > mp->vlen)
         return(0);
   }
   else if (mp->ku > 4 && KB != kb)
      return(0);
   return(KB);
}

/* procedure 28 */
ATL_mmnode_t *bestNearMNK(char pre, int verb, ATL_mmnode_t *tb,
                          int mb, int nb, int kb, int tmflag, int beta,
                          double runbon)
/*
 * Tries all mm in tb (try base), at sizes near mb,nb,kb.
 * OVERWRITES: mbB,nbB,kbB,mflop[0] for all kerns!
 * RETURNS: clone of best kernel's MMNODE
 */
{
   ATL_mmnode_t *mp, *mpB=NULL;
   int KRUNB=0;
   const double cmppen=1.0 / runbon;
   double mfB=0.0;

   for (mp=tb; mp; mp = mp->next)
   {
      double mf, mfmul;
      const unsigned int mu=mp->mu, nu=mp->nu, ku=mp->ku;
      int MB, NB, KB;
      char *nm;

      if (!kernWorksThisCaseK(mp, kb))  /* don't try illegal kerns */
         continue;
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME)) /* get runtime bonus */
         mfmul = (KRUNB) ? 1.0 : runbon;
      else
         mfmul = (KRUNB) ? cmppen : 1.0;
      mp->mbB = MB = (mb > mu) ? (mb/mu)*mu : mu;
      mp->nbB = NB = (nb > nu) ? (nb/nu)*nu : nu;
      mp->kbB = KB = ((kb+ku-1)/ku)*ku;
      mf = TimeMMKernel(verb, tmflag, mp, pre, MB, NB, KB, beta, 0, -1);
      mp->mflop[0] = mf;
      printf("      ID=%d '%s': B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), MB, NB, KB, mf);
      if (mf*mfmul > mfB)
      {
         KRUNB = FLAG_IS_SET(mp->flag, MMF_KRUNTIME);
         mfB = mf;
         mpB = mp;
      }
   }
   if (mpB)
   {
      const unsigned int mu=mpB->mu, nu=mpB->nu, ku=mpB->ku;
      mpB = CloneMMNode(mpB);
      mpB->kbB = kb;
      mpB->mbB = (mb > mu) ? ((mb/mu)*mu) : mu;
      mpB->nbB = (nb > nu) ? ((nb/nu)*nu) : nu;
      mpB->mflop[0] = mfB;
   }
   GetMMLabelName(pre, NULL);
   return(mpB);
}

/* procedure 29 */
ATL_mmnode_t *ElimSlowKern(char pre, int verb, ATL_mmnode_t *tb,
                           int mb, int nb, int kb, int tmflag, int beta,
                           double tol)
/*
 * Times all kernels at size close to mb,nb,kb, the eliminates all kernels
 * not within tolerance of best
 */
{
   ATL_mmnode_t *mp;
   double mfB;

   mp = bestNearMNK(pre, verb, tb, mb, nb, kb, tmflag, beta, 1.0);
   mfB = mp->mflop[0];
   KillMMNode(mp);
   mp = tb;
   while (mp)
   {
      ATL_mmnode_t *nxt = mp->next;
      if (mp->mflop[0] < mfB*tol)
      {
         tb = RemoveMMNodeFromQ(tb, mp);
         KillMMNode(mp);
      }
      mp = nxt;
   }
   return(tb);
}

/* procedure 30 */
ATL_mmnode_t *bestNearSquare(char pre, int verb, ATL_mmnode_t *tb,
                             int kb, int tmflag, int beta, double runbon)
/*
 * Tries all mm in tb (try base), at near-square sizes around kb.
 * RETURNS: clone of best kernel's MMNODE
 */
{
   ATL_mmnode_t *mp, *mpB=NULL;
   const double cmppen=1.0 / runbon;
   double mfB=0.0;
   int MU=kb, NU=kb, KRUNB=0;
/*
 * Despite name, will use M/N B at least the size of largest unrolling,
 * otherwise large unrollings get huge MFLOP advantage in head-to-head
 * timings for small problems.
 */
   for (mp=tb; mp; mp = mp->next)
   {
      MU = Mmax(mp->mu, MU);
      NU = Mmax(mp->nu, NU);
   }
   for (mp=tb; mp; mp = mp->next)
   {
      double mf, mfmul;
      const unsigned int mu=mp->mu, nu=mp->nu, ku=mp->ku;
      int MB, NB, KB;
      char *nm;

      if (!kernWorksThisCaseK(mp, kb))  /* don't try illegal kerns */
         continue;
      if (FLAG_IS_SET(mp->flag, MMF_KRUNTIME)) /* get runtime bonus */
         mfmul = (KRUNB) ? 1.0 : runbon;
      else
         mfmul = (KRUNB) ? cmppen : 1.0;
      MB = ((MU+mu-1)/mu)*mu;
      NB = ((NU+nu-1)/nu)*nu;
      KB = ((kb+ku-1)/ku)*ku;
      mf = TimeMMKernel(verb, tmflag, mp, pre, MB, NB, kb, beta, 0, -1);
      printf("      ID=%d '%s': B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), MB, NB, KB, mf);
      if (mf*mfmul > mfB)
      {
         KRUNB = FLAG_IS_SET(mp->flag, MMF_KRUNTIME);
         mfB = mf;
         mpB = mp;
      }
   }
   if (mpB)
   {
      const unsigned int mu=mpB->mu, nu=mpB->nu, ku=mpB->ku;
      mpB = CloneMMNode(mpB);
      mpB->kbB = kb;
      mpB->mbB = ((MU+mu-1)/mu)*mu;
      mpB->nbB = ((NU+nu-1)/nu)*nu;
      mpB->mflop[0] = mfB;
   }
   GetMMLabelName(pre, NULL);
   return(mpB);
}

/* procedure 31 */
int MMVaryDim(int verb, char pre, int tflag, int MDIM, int MD, int beta,
                ATL_mmnode_t *mp)
/*
 * Tries improving performance by increasing a M or N dimension up to MD.
 * UPDATES: mflop[0] with time ([0] not overwritten in case needed).
 * RETURNS: best blocking for dimension.
 */
{
   int mb = mp->mbB, nb=mp->nbB, kb=mp->kbB;
   int U=(MDIM)?mp->mu:mp->nu, dB, d;
   double mfB=0.0;
   if (U == 3)
      U = 6;
   else if (U < 4)
      U = 4;
   dB = U;
   for (d=U; d <= MD; d += U)
   {
      double mf;
      if (MDIM)
         mb = d;
      else
         nb = d;
      mf = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kb, beta, 0, -1);
      printf("      ID=%d '%s': B=(%d,%d,%d) mf=%.2f\n", mp->ID,
             GetMMLabelName(pre, mp), mb, nb, kb, mf);
      if (mf > mfB)
      {
         mfB = mf;
         dB = d;
      }
   }
   mp->mflop[1] = mfB;
   return(dB);
}

/* procedure 32 */
#define NFUT  24   /* max # of futile expansion tries */
#define BLK1  60
#define B1PEN 1.0  /* no penalty for D <= BLK1 */
#define BLK2  120
#define B2PEN 0.998 /* 0.2% penalty for BLK1 < D <= BLK2 */
#define BLK3  240
#define B3PEN 0.995 /* 0.5% penalty for BLK2 < D <= BLK3 */
#define BLK4  480
#define B4PEN 0.990 /* 1.0% penalty for BLK3 < D <= BLK4 */
#define B5PEN 0.970 /* 3.0% penalty for D > BLK4 */

double MMExpandMN_K
   (int verb, char pre, int tflag, int beta, unsigned long sz,
    int (*wrkSetOK)(unsigned long, ATL_mmnode_t*,int,int,int),
    const int maxB, ATL_mmnode_t *mp)
/*
 * Assume U=LCM(MU,NU), this search startes at MB=NB=U,KB=KU and then tries
 * increasing M,N,K and picks whichever provides the greatest speedup.
 * Continues expanding MB,NB,KB in this way until perf drops by around 3%.
 * maxB must be > 0, and it sets inclusive limit on growth of MB,NB.
 * Once we have stopped expanding each by U, consider all MU,NU,KU < U as well.
 * RETURNS: mflop for best MB,NB,KB which is given in mp->(mbB,nbB,kbB).
 */
{
   const unsigned int mu=mp->mu, nu=mp->nu;
   unsigned int MU, NU, KU, ku;
   unsigned int maxKB = Mmax(240,maxB+maxB);
   int mbB, nbB, kbB, mb, nb, kb, maxDB, maxD;
   int FASTER=0, nfut;
   double mfB, penB=B1PEN;

   if (mp->ID && mp->kbmax && mp->kbmax < maxKB)
      maxKB = mp->kbmax;
   NU = ATL_iLCM(mu, nu);
   ku = mp->ku;
   if (!FLAG_IS_SET(mp->flag, MMF_KUISKB) || !mp->ID) /* gen can adapt KB */
   {  /* genned KUISKB or non-KUISKB can always adapt KU as needed */
      KU = ATL_iLCM(NU, ku);
      if (KU > 16 || KU > maxKB) /* too big, make smaller if possible */
      {
         if (ku <= 8)
         {
            for (KU=ku; KU < 8; KU += ku);
         }
         else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
            KU = mp->vlen;
         else
            KU = NU;
         kbB = kb = KU;
      }
      else
         NU = kbB = kb = KU;
   }
   else if (mp->kbmin > 1)
   {                               /* reach here --> user KUISKB kern */
      if (mp->kbmax == mp->kbmin)
         kbB = kb = KU = mp->kbmax;
      else
      {
         unsigned const int inc=NU;
         KU = ATL_iLCM(NU, ku);
         if (KU > 16 || KU > maxKB) /* too big, make smaller if possible */
         {
            for (KU=ku; KU < Mmax(8, mp->kbmin); KU += ku);
            kbB = kb = KU; /* NOT IN LOOP! */
         }
         else
            NU = kbB = kb = KU;
      }
   }
   else
   {
      ku = 1;
      kbB = kb = KU = NU;
   }
   maxDB = mb = nb = mbB = nbB = MU = NU;

   mfB = TimeMMKernel(verb, tflag, mp, pre, mbB, nbB, kb, beta, 0, -1);
   printf("      ID=%d : B=(%u,%u,%u) maxB=%d mf=%.2f\n", mp->ID,
          mbB, nbB, kbB, maxB, mfB);

   AGAIN:
      nfut = 0;
      do   /* loop to check that best of M/N expansion improved performance */
      {
         double mfM, mfN, mfK, mf;
         int mbN=mb+MU, nbN=nb+NU, kbN=kb+KU;
         int MOK, NOK, KOK;
         int ku0=mp->ku, kbmin0=mp->kbmin, kbmax0=mp->kbmax;
         double penM, penN, penK;

         maxD = Mmax(nbB, kbB);
         maxD = Mmax(maxD, mbN);
         if (maxD <= BLK1)
            penM = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penM = B2PEN;
            else if (maxD <= BLK3)
               penM = B3PEN;
            else if (maxD <= BLK4)
               penM = B4PEN;
            else
               penM = B5PEN;
            if (penM == penB)
               penM = 0.998; /* 0.2% penalty for any expansion */
            else
               penM /= penB;
         }
         maxD = Mmax(mbB, kbB);
         maxD = Mmax(maxD, nbN);
         if (maxD <= BLK1)
            penN = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penN = B2PEN;
            else if (maxD <= BLK3)
               penN = B3PEN;
            else if (maxD <= BLK4)
               penN = B4PEN;
            else
               penN = B5PEN;
            if (penN == penB)
               penN = 0.998; /* 0.2% penalty for any expansion */
            else
               penN /= penB;
         }
         maxD = Mmax(mbB, nbB);
         maxD = Mmax(maxD, kbN);
         if (maxD <= BLK1)
            penK = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penK = B2PEN;
            else if (maxD <= BLK3)
               penK = B3PEN;
            else if (maxD <= BLK4)
               penK = B4PEN;
            else
               penK = B5PEN;
            if (penK == penB)
               penK = 0.998; /* 0.2% penalty for any expansion */
            else
               penK /= penB;
         }
         mfM = mfN = mfK = mf = 0.0;
         if (sz && wrkSetOK)
         {
            MOK = wrkSetOK(sz, mp, mbN, nb, kb);
            NOK = wrkSetOK(sz, mp, mb, nbN, kb);
            KOK = wrkSetOK(sz, mp, mb, nb, kbN);
            if (!(MOK|NOK|KOK))
               break;
         }
         else
            MOK = NOK = KOK = 1;
         if (maxKB && kbN > maxKB)
         {
            if (mbN > maxB && nbN > maxB)
               break;
            mfK = 0.0;
         }
         else if (!NOK)
            mfK = 0.0;
         else
         {
            mfK = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kbN, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",mp->ID,mb,nb,kbN,mfK);
            if (!maxB || mbN <= maxB || nbN <= maxB)
               printf("\n");
         }
         if (maxB && nbN > maxB)
            mfN = 0.0;
         else
         {
            mfN = TimeMMKernel(verb, tflag, mp, pre, mb, nbN, kb, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",mp->ID,mb,nbN,kb,mfN);
            if (!maxB || mbN <= maxB)
               printf("\n");
         }
         if (maxB && mbN > maxB)
            mfM = 0.0;
         else
         {
            mfM = TimeMMKernel(verb, tflag, mp, pre, mbN, nb, kb, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",  mp->ID,mbN,nb,kb,mfM);
         }

         FASTER = 0;
         if (mfK*penK-mfB >= mfM*penM-mfB &&
             mfK*penK-mfB >= mfN*penN-mfB)
         {
            mf = mfK;
            kb = kbN;
            if (penK*mfK >= mfB)
            {
               printf(" --> INC KB TO %u!\n", kbN);
               kbB = kbN;
               mfB = mfK;
               penB = penK;
               maxDB = Mmax(maxDB, kbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else if (mfN*penN-mfB >= mfM*penM-mfB)
         {
            mf = mfN;
            nb = nbN;
            if (penN*mfN >= mfB)
            {
               printf(" --> INC NB TO %u!\n", nbN);
               nbB = nbN;
               mfB = mfN;
               penB = penN;
               maxDB = Mmax(maxDB, nbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else
         {
            mf = mfM;
            mb = mbN;
            if (penM*mfM >= mfB)
            {
               printf(" --> INC MB TO %u!\n", mbN);
               mbB = mbN;
               mfB = mfM;
               penB = penM;
               maxDB = Mmax(maxDB, mbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         if (FLAG_IS_SET(mp->flag, MMF_KUISKB) && mp->kbmin != kbB)
         {
            if (!mp->ID && mp->ku != 1)
               mp->ku = kbB;
            mp->kbmin = mp->kbmax = kbB;
            if (mp->genstr)
            {
               free(mp->genstr);
               mp->genstr = MMGetGenString(pre, mp);
            }
         }
         if (!FASTER)
         {
            mf = Mmax(mfM, mfN);
            mf = Mmax(mf, mfK);
            FASTER = (nfut < NFUT);
            printf(" --> SLOWDOWN of %.4f\n", mf/mfB);
            nfut++;
         }
      }
      while (FASTER);
   if (MU != mu || NU != nu || KU != ku)
   {
      printf("\n   REFINEMENT SEARCH B=(%d,%d,%d), mf=%.2f\n", mbB, nbB, kbB,
             mfB);
      mb = mbB;
      nb = nbB;
      kb = kbB;
      MU = mu;
      NU = nu;
      KU = ku;
      goto AGAIN;
   }

   mp->mbB = mbB;
   mp->nbB = nbB;
   mp->kbB = kbB;
   mp->mflop[0] = mfB;
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB) && !mp->ID)
   {
      mp->ku = kbB;
      if (mp->genstr)
      {
         free(mp->genstr);
         mp->genstr = NULL;
      }
   }
   assert(!MMKernelFailsTest(pre, mbB, nbB, kbB, beta, mp));
   return(mfB);
}
double MMExpandMNK
   (int verb, char pre, int tflag, int beta, unsigned long sz,
    int (*wrkSetOK)(unsigned long, ATL_mmnode_t*,int,int,int),
    const int maxB, ATL_mmnode_t *mp)
/*
 * Assume U=LCM(MU,NU), this search startes at MB=NB=U,KB=KU and then tries
 * increasing M,N,K and picks whichever provides the greatest speedup.
 * Continues expanding MB,NB,KB in this way until perf drops by around 3%.
 * If maxB is non-zero, then it sets inclusive limit on growth of MB,NB,KB.
 * Once we have stopped expanding each by U, consider all MU,NU,KU < U as well.
 * RETURNS: mflop for best MB,NB,KB which is given in mp->(mbB,nbB,kbB).
 */
{
   const unsigned int mu=mp->mu, nu=mp->nu;
   unsigned int MU, NU, KU, ku;
   unsigned int maxKB = maxB;
   int mbB, nbB, kbB, mb, nb, kb, maxDB, maxD;
   int FASTER=0, nfut;
   double mfB, penB=B1PEN;

   if (mp->ID && mp->kbmax && mp->kbmax < maxKB)
      maxKB = mp->kbmax;
   NU = ATL_iLCM(mu, nu);
   ku = mp->ku;
   if (!FLAG_IS_SET(mp->flag, MMF_KUISKB) || !mp->ID) /* gen can adapt KB */
   {  /* genned KUISKB or non-KUISKB can always adapt KU as needed */
      KU = ATL_iLCM(NU, ku);
      if (KU > 16 || KU > maxKB) /* too big, make smaller if possible */
      {
         if (ku <= 8)
         {
            for (KU=ku; KU < 8; KU += ku);
         }
         else if (FLAG_IS_SET(mp->flag, MMF_KVEC))
            KU = mp->vlen;
         else
            KU = NU;
         kbB = kb = KU;
      }
      else
         NU = kbB = kb = KU;
   }
   else if (mp->kbmin > 1)
   {                               /* reach here --> user KUISKB kern */
      if (mp->kbmax == mp->kbmin)
         kbB = kb = KU = mp->kbmax;
      else
      {
         unsigned const int inc=NU;
         KU = ATL_iLCM(NU, ku);
         if (KU > 16 || KU > maxKB) /* too big, make smaller if possible */
         {
            for (KU=ku; KU < Mmax(8, mp->kbmin); KU += ku);
            kbB = kb = KU; /* NOT IN LOOP! */
         }
         else
            NU = kbB = kb = KU;
      }
   }
   else
   {
      ku = 1;
      kbB = kb = KU = NU;
   }
   maxDB = mb = nb = mbB = nbB = MU = NU;

   mfB = TimeMMKernel(verb, tflag, mp, pre, mbB, nbB, kb, beta, 0, -1);
   printf("      ID=%d : B=(%u,%u,%u) maxB=%d mf=%.2f\n", mp->ID,
          mbB, nbB, kbB, maxB, mfB);

   AGAIN:
      nfut = 0;
      do   /* loop to check that best of M/N expansion improved performance */
      {
         double mfM, mfN, mfK, mf;
         int mbN=mb+MU, nbN=nb+NU, kbN=kb+KU;
         int MOK, NOK, KOK;
         int ku0=mp->ku, kbmin0=mp->kbmin, kbmax0=mp->kbmax;
         double penM, penN, penK;

         maxD = Mmax(nbB, kbB);
         maxD = Mmax(maxD, mbN);
         if (maxD <= BLK1)
            penM = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penM = B2PEN;
            else if (maxD <= BLK3)
               penM = B3PEN;
            else if (maxD <= BLK4)
               penM = B4PEN;
            else
               penM = B5PEN;
            if (penM == penB)
               penM = 0.998; /* 0.2% penalty for any expansion */
            else
               penM /= penB;
         }
         maxD = Mmax(mbB, kbB);
         maxD = Mmax(maxD, nbN);
         if (maxD <= BLK1)
            penN = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penN = B2PEN;
            else if (maxD <= BLK3)
               penN = B3PEN;
            else if (maxD <= BLK4)
               penN = B4PEN;
            else
               penN = B5PEN;
            if (penN == penB)
               penN = 0.998; /* 0.2% penalty for any expansion */
            else
               penN /= penB;
         }
         maxD = Mmax(mbB, nbB);
         maxD = Mmax(maxD, kbN);
         if (maxD <= BLK1)
            penK = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penK = B2PEN;
            else if (maxD <= BLK3)
               penK = B3PEN;
            else if (maxD <= BLK4)
               penK = B4PEN;
            else
               penK = B5PEN;
            if (penK == penB)
               penK = 0.998; /* 0.2% penalty for any expansion */
            else
               penK /= penB;
         }
         mfM = mfN = mfK = mf = 0.0;
         if (sz && wrkSetOK)
         {
            MOK = wrkSetOK(sz, mp, mbN, nb, kb);
            NOK = wrkSetOK(sz, mp, mb, nbN, kb);
            KOK = wrkSetOK(sz, mp, mb, nb, kbN);
            if (!(MOK|NOK|KOK))
               break;
         }
         else
            MOK = NOK = KOK = 1;
         if (maxKB && kbN > maxKB)
         {
            if (mbN > maxB && nbN > maxB)
               break;
            mfK = 0.0;
         }
         else if (!NOK)
            mfK = 0.0;
         else
         {
            mfK = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kbN, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",mp->ID,mb,nb,kbN,mfK);
            if (!maxB || mbN <= maxB || nbN <= maxB)
               printf("\n");
         }
         if (maxB && nbN > maxB)
            mfN = 0.0;
         else
         {
            mfN = TimeMMKernel(verb, tflag, mp, pre, mb, nbN, kb, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",mp->ID,mb,nbN,kb,mfN);
            if (!maxB || mbN <= maxB)
               printf("\n");
         }
         if (maxB && mbN > maxB)
            mfM = 0.0;
         else
         {
            mfM = TimeMMKernel(verb, tflag, mp, pre, mbN, nb, kb, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f",  mp->ID,mbN,nb,kb,mfM);
         }

         FASTER = 0;
         if (mfK*penK-mfB >= mfM*penM-mfB &&
             mfK*penK-mfB >= mfN*penN-mfB)
         {
            mf = mfK;
            kb = kbN;
            if (penK*mfK >= mfB)
            {
               printf(" --> INC KB TO %u!\n", kbN);
               kbB = kbN;
               mfB = mfK;
               penB = penK;
               maxDB = Mmax(maxDB, kbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else if (mfN*penN-mfB >= mfM*penM-mfB)
         {
            mf = mfN;
            nb = nbN;
            if (penN*mfN >= mfB)
            {
               printf(" --> INC NB TO %u!\n", nbN);
               nbB = nbN;
               mfB = mfN;
               penB = penN;
               maxDB = Mmax(maxDB, nbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else
         {
            mf = mfM;
            mb = mbN;
            if (penM*mfM >= mfB)
            {
               printf(" --> INC MB TO %u!\n", mbN);
               mbB = mbN;
               mfB = mfM;
               penB = penM;
               maxDB = Mmax(maxDB, mbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         if (FLAG_IS_SET(mp->flag, MMF_KUISKB) && mp->kbmin != kbB)
         {
            if (!mp->ID && mp->ku != 1)
               mp->ku = kbB;
            mp->kbmin = mp->kbmax = kbB;
            if (mp->genstr)
            {
               free(mp->genstr);
               mp->genstr = MMGetGenString(pre, mp);
            }
         }
         if (!FASTER)
         {
            mf = Mmax(mfM, mfN);
            mf = Mmax(mf, mfK);
            FASTER = (nfut < NFUT);
            printf(" --> SLOWDOWN of %.4f\n", mf/mfB);
            nfut++;
         }
      }
      while (FASTER);
   if (MU != mu || NU != nu || KU != ku)
   {
      printf("\n   REFINEMENT SEARCH B=(%d,%d,%d), mf=%.2f\n", mbB, nbB, kbB,
             mfB);
      mb = mbB;
      nb = nbB;
      kb = kbB;
      MU = mu;
      NU = nu;
      KU = ku;
      goto AGAIN;
   }

   mp->mbB = mbB;
   mp->nbB = nbB;
   mp->kbB = kbB;
   mp->mflop[0] = mfB;
   if (FLAG_IS_SET(mp->flag, MMF_KUISKB) && !mp->ID)
   {
      mp->ku = kbB;
      if (mp->genstr)
      {
         free(mp->genstr);
         mp->genstr = NULL;
      }
   }
   assert(!MMKernelFailsTest(pre, mbB, nbB, kbB, beta, mp));
   return(mfB);
}


/* procedure 33 */
void MMExpandMorN(int verb, char pre, int tflag, int beta, ATL_mmnode_t *mp,
                  int EXN, unsigned int maxD)
{
   const unsigned int mu=EXN?0:mp->mu, nu=EXN?mp->nu:0, kb=mp->kbB;
   unsigned int mb=mp->mbB, nb=mp->nbB, nbB, U, i;
   double mfB;

   U = (EXN) ? mp->nu : mp->mu;
   U = (U==1 || U == 2) ? 4 : U;
   mb = (EXN) ? mb : U;
   nb = (EXN) ? U : nb;
   nbB = (EXN) ? nb : mb;
   mfB = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kb, beta, 0, -1);
   printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID,mb, nb, kb, mfB);
   for (i=U+U; i <= maxD; i += U)
   {
      double mf;
      if (EXN)
         nb = i;
      else
         mb = i;
      mf = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kb, beta, 0, -1);
      printf("      ID=%d : B=(%d,%d,%d) mf=%.2f", mp->ID,mb, nb, kb, mf);
      if (mf > mfB)
      {
         printf(" --> %.5f faster\n", mf/mfB);
         mfB = mf;
         if (EXN)
	    nbB = nb;
         else
	    nbB = mb;
      }
      else
         printf(" --> %.5f slower\n", mf/mfB);
   }
   mp->mflop[0] = mfB;
   if (EXN)
      mp->nbB = nbB;
   else
      mp->mbB = nbB;
}
/*
/* procedure 34 */
int MMExpandK(int verb, char pre, int tflag, int beta, ATL_mmnode_t *mp,
              unsigned int maxK)
/*
 * Times expanding K until perf drops off by around 2%.  If maxK != 0, also
 * stops expanding once K is reached.
 */
{
   double mfB, mf;
   const unsigned int ku = mp->ku, mb=mp->mbB, nb=mp->nbB;
   unsigned int kbB, k, U;
   if (!maxK)
      maxK = 1024;
   if (mp->kbmax)
      maxK = Mmin(maxK, mp->kbmax);
/*
 * In initial search, recursive doubling
 */
   printf("   %d %s: FINDING KB FOR B=(%d,%d):\n", mp->ID,
          GetMMLabelName(pre, mp), mb, nb);
   mfB = TimeMMKernel(verb, tflag, mp, pre, mb, nb, ku, beta, 0, -1);
   printf("      KB=%u, mf=%.2f\n", ku, mfB);
   k = kbB = ku;
   do
   {
      k += k;
      mf = TimeMMKernel(verb, tflag, mp, pre, mb, nb, k, beta, 0, -1);
      printf("      KB=%u, mf=%.2f\n", k, mf);
      if (mf > mfB)
      {
         mfB = mf;
         kbB = k;
      }
   }
   while (k+k <= maxK && mf*1.02 > mfB);
   if (ku < 4)
   {
      U = ATL_iLCM(mp->nu, ku);
      if (U > 8)
         U = ATL_iLCM(mp->mu, ku);
      if (U > 8)
         U = ku+ku;
   }
   else
      U = ku;
   maxK++;
   k = kbB+kbB;
   maxK = Mmin(maxK, k);
   printf("\n   INITIAL KB=%u, refine [%u,%u) steps of %u:\n", kbB, kbB>>1,
          maxK, U);
   for (k=(kbB>>1)+U; k < maxK; k += U)
   {
      mf = TimeMMKernel(verb, tflag, mp, pre, mb, nb, k, beta, 0, -1);
      printf("      KB=%u, mf=%.2f\n", k, mf);
      if (mf > mfB)
      {
         mfB = mf;
         kbB = k;
      }
   }
   mp->mflop[0] = mfB;
   mp->kbB = kbB;
   printf("   %d %s, B=(%d,%d,%d), mf=%.2f\n", mp->ID, GetMMLabelName(pre, mp),
          mb, nb, kbB, mfB);
   GetMMLabelName(pre, NULL);
   return(kbB);
}

/* procedure 35 */
double MMExpandMN(int verb, char pre, int tflag, int beta, int maxB,
                  ATL_mmnode_t *mp)
/*
 * Assume U=LCM(MU,NU), this search startes at MB=NB=U, and then tries
 * increasing M and N, and picks whichever provides the greatest speedup.
 * Continues expanding MB,NB in this way until performance drops by around 3%.
 * Once we have stopped expanding each by U, consider all MU,NU < U as well.
 * KB is not tuned: all timings use KB=mp->kbB.
 * RETURNS: mflop for best MB,NB, which is given in mp->mbB, mp->nbB.
 */
{
   const unsigned int mu=mp->mu, nu=mp->nu, MB0=mp->mbB, NB0=mp->nbB;
   unsigned int MU, NU;
   int kb=mp->kbB, mbB, nbB, mb, nb, maxD, maxDB;
   int FASTER=0, nfut;
   double mfB, penB=B1PEN;

   maxDB = mb = nb = mbB = nbB = MU = NU = ATL_iLCM(mu, nu);
   mfB = TimeMMKernel(verb, tflag, mp, pre, mbB, nbB, kb, beta, 0, -1);
   printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID, mbB, nbB, kb, mfB);

   AGAIN:
      nfut = 0;
      do   /* loop to check that best of M/N expansion improved performance */
      {
         double mfM, mfN, mf=0.0, penM, penN;
         int mbN=mb+MU, nbN=nb+NU;

         maxD = Mmax(nbB, kb);
         maxD = Mmax(maxD, mbN);
         if (maxD <= BLK1)
            penM = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penM = B2PEN;
            else if (maxD <= BLK3)
               penM = B3PEN;
            else if (maxD <= BLK4)
               penM = B4PEN;
            else
               penM = B5PEN;
            if (penM == penB)
               penM = 0.998; /* 0.2% penalty for any expansion */
            else
               penM /= penB;
         }
         maxD = Mmax(mbB, kb);
         maxD = Mmax(maxD, nbN);
         if (maxD <= BLK1)
            penN = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penN = B2PEN;
            else if (maxD <= BLK3)
               penN = B3PEN;
            else if (maxD <= BLK4)
               penN = B4PEN;
            else
               penN = B5PEN;
            if (penN == penB)
               penN = 0.998; /* 0.2% penalty for any expansion */
            else
               penN /= penB;
         }

         if (maxB && nbN > maxB && mbN > maxB)
            break;
         if (maxB && nbN > maxB)
            mfN = 0.0;
         else
         {
            mfN = TimeMMKernel(verb, tflag, mp, pre, mb, nbN, kb, beta, 0, -1);
            printf("      ID=%d : B=(%d,%d,%d) mf=%.2f", mp->ID,mb, nbN,
                   kb, mfN);
            if (!maxB || mbN <= maxB)
               printf("\n");
         }
         if (maxB && mbN > maxB)
            mfM = 0.0;
         else
         {
            mfM = TimeMMKernel(verb, tflag, mp, pre, mbN, nb, kb, beta, 0, -1);
            printf("      ID=%d ': B=(%d,%d,%d) mf=%.2f",  mp->ID,mbN,nb,
                   kb, mfM);
         }

         FASTER = 0;
         if (mfN*penN-mfB >= mfM*penN-mfB)
         {
            mf = mfN;
            nb = nbN;
            if (mfN*penN > mfB)
            {
               printf(" --> INC NB TO %u!\n", nbN);
               nbB = nbN;
               mfB = mfN;
               penB = penN;
               maxDB = Mmax(maxDB, nbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else
         {
            mf = mfM;
            mb = mbN;
            if (mfM*penM >= mfB)
            {
               printf(" --> INC MB TO %u!\n", mbN);
               mbB = mbN;
               mfB = mfM;
               penB = penM;
               maxDB = Mmax(maxDB, mbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         if (!FASTER)
         {
            nfut++;
            FASTER = (nfut < 5);
            printf(" --> SLOWDOWN of %.4f\n", mf/mfB);
         }
      }
      while (FASTER);
   if (MU != mu || NU != nu)
   {
      printf("\n");
      mb = mbB;
      nb = nbB;
      MU = mu;
      NU = nu;
      goto AGAIN;
   }
/*
 * Now retime original MB/NB if they make sense
 */
   if ((MB0/mu)*mu == MB0 && (NB0/nu)*nu == NB0)
   {
      double mf;
      mf = TimeMMKernel(verb, tflag|1, mp, pre, MB0, NB0, kb, beta, 0, -1);
      printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID,MB0,NB0,kb,mf);
      if (mf > mfB)
      {
         mfB = mf;
         mbB = MB0;
         nbB = NB0;
      }
   }
   mp->mbB = mbB;
   mp->nbB = nbB;
   mp->mflop[0] = mfB;
   return(mfB);
}

/* procedure 36 */
double MMExpandNK(int verb, char pre, int tflag, int beta, ATL_mmnode_t *mp)
/*
 * Expand either N=K or M=K, with other dims unrestrained until max perf
 * is reached.  If most sig bit in tflag is set, expand N=K, else M=K.
 * RETURNS: mflop for best MB,NB, which is given in mp->mbB, mp->nbB.
 */
{
   const unsigned int mu=mp->mu, nu=mp->nu;
   unsigned int MB0=mp->mbB, NB0=mp->nbB, KB0=mp->kbB, ku=mp->ku;
   unsigned int KBMAX = mp->kbmax ? mp->kbmax : 512;
   unsigned int KBMIN = mp->kbmin ? mp->kbmin : 1;
   unsigned int ADJKUKB=0;
   unsigned int MU, NU;
   const int NK=(tflag>>31)&1;
   int kb=mp->kbB, mbB, nbB, kbB, mb, nb, maxD, maxDB;
   int FASTER=0, nfut=0;
   double mfB, penB=B1PEN;
   if (mp->flag & (1<<MMF_KUISKB))
   {
      if (!mp->ID)
         ADJKUKB= ku = 1;
      else if (!mp->kbmax || mp->kbmax > mp->kbmin)
         ADJKUKB = 1;
   }
   if (NK)
   {

      if (ADJKUKB)
      {
         if (!mp->ID)
         {
            kb = nu;
            mp->kbB = mp->ku = nu;
         }
         else /* user contributed kernel */
         {
            if (mp->kbmin)
            {
               ku = ((mp->kbmin+nu-1)/nu)*nu;
               kb = ku;
            }
            else
            {
               ku = 1;
               kb = nu;
            }
            if (mp->kbmax)
               assert(ku < mp->kbmax);
         }
      }
      else
         kb = ATL_iLCM(nu, ku);
      NU = nb = kb;
      MU = mb = mu;
      if (NB0 != KB0)
         MB0 = NB0 = KB0 = 0;
   }
   else
   {
      if (ADJKUKB)
      {
         if (!mp->ID)
         {
            kb = nu;
            mp->kbB = mp->ku = nu;
         }
         else /* user contributed kernel */
         {
            if (mp->kbmin)
            {
               ku = ((mp->kbmin+nu-1)/nu)*nu;
               kb = ku;
            }
            else
            {
               ku = 1;
               kb = nu;
            }
            if (mp->kbmax)
               assert(ku < mp->kbmax);
         }
      }
      else
         kb = ATL_iLCM(mu, ku);
      MU = mb = kb;
      NU = nb = nu;
      if (MB0 != KB0)
         MB0 = NB0 = KB0 = 0;
   }
   mbB = mb;
   nbB = nb;
   kbB = kb;
   maxDB = Mmax(mbB, nbB);
   maxDB = Mmax(maxD, kbB);

   mfB = TimeMMKernel(verb, tflag, mp, pre, mb, nb, kb, beta, 0, -1);
   printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID, mbB, nbB, kb, mfB);

      do   /* loop to check that best of M/N expansion improved performance */
      {
         double mfM, mfN, mf=0.0, penM, penN;
         int mbN=mb+MU, nbN=nb+NU;
         const unsigned int kbm = (NK) ? nb : mbN, kbn = (NK) ? nbN : mb;

         if (kbn > KBMAX && kbm > KBMAX)
             break;
         maxD = Mmax(nbB, kbm);
         maxD = Mmax(maxD, mbN);
         if (maxD <= BLK1)
            penM = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penM = B2PEN;
            else if (maxD <= BLK3)
               penM = B3PEN;
            else if (maxD <= BLK4)
               penM = B4PEN;
            else
               penM = B5PEN;
            if (penM == penB)
               penM = 0.998; /* 0.2% penalty for any expansion */
            else
               penM /= penB;
         }
         maxD = Mmax(mbB, kbn);
         maxD = Mmax(maxD, nbN);
         if (maxD <= BLK1)
            penN = B1PEN;
         else
         {
            if (maxD <= BLK2)
               penN = B2PEN;
            else if (maxD <= BLK3)
               penN = B3PEN;
            else if (maxD <= BLK4)
               penN = B4PEN;
            else
               penN = B5PEN;
            if (penN == penB)
               penN = 0.998; /* 0.2% penalty for any expansion */
            else
               penN /= penB;
         }
         if (kbn <= KBMAX)
         {
            if (ADJKUKB)
               mp->kbB = mp->ku = kbn;
            mfN = TimeMMKernel(verb, tflag, mp, pre, mb, nbN, kbn, beta, 0, -1);
         }
         else
            mfN = 0;
         printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID,mb, nbN,kbn,mfN);
         if (kbm <= KBMAX)
         {
            if (ADJKUKB)
               mp->kbB = mp->ku = kbn;
            mfM = TimeMMKernel(verb, tflag, mp, pre, mbN, nb, kbm, beta, 0, -1);
         }
         else
            mfM = 0;
         printf("      ID=%d ': B=(%d,%d,%d) mf=%.2f",  mp->ID,mbN,nb, kbm,mfM);

         FASTER = 0;
         if (mfN*penM-mfB >= mfM*penM-mfB)
         {
            mf = mfN;
            nb = nbN;
            if (mfN*penN >= mfB)
            {
               printf(" --> INC NB TO %u!\n", nbN);
               mbB = mb;
               nbB = nbN;
               mfB = mfN;
               kbB = kbn;
               penB = penN;
               maxDB = Mmax(maxDB, nbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         else
         {
            mf = mfM;
            mb = mbN;
            if (mfM*penM >= mfB)
            {
               printf(" --> INC MB TO %u!\n", mbN);
               mbB = mbN;
               nbB = nb;
               mfB = mfM;
               kbB = kbm;
               penB = penM;
               maxDB = Mmax(maxDB, mbB);
               FASTER = 1;
               nfut = 0;
            }
         }
         if (!FASTER)
         {
            nfut++;
            FASTER = (mf*1.03 > mfB) && (nfut < NFUT);
            printf(" --> SLOWDOWN of %.4f\n", mf/mfB);
         }
      }
      while (FASTER);
/*
 * Now retime original MB/NB/KB if they make sense
 */
   if (MB0 && NB0 && KB0 && (MB0/mu)*mu == MB0 && (NB0/nu)*nu == NB0 &&
       (KB0/ku)*ku == KB0)
   {
      double mf;
      if (ADJKUKB)
         mp->kbB = mp->ku = KB0;
      mf = TimeMMKernel(verb, tflag|1, mp, pre, MB0, NB0, KB0, beta, 0, -1);
      printf("      ID=%d : B=(%d,%d,%d) mf=%.2f\n", mp->ID,MB0,NB0,KB0,mf);
      if (mf > mfB)
      {
         mfB = mf;
         mbB = MB0;
         nbB = NB0;
         kbB = KB0;
      }
   }
   if (mbB && nbB && kbB)
   {
      mp->mbB = mbB;
      mp->nbB = nbB;
      mp->kbB = kbB;
      if (ADJKUKB)
         mp->ku = kbB;
      mp->mflop[0] = mfB;
   }
   return(mfB);
}

#undef NFUT
#undef BLK1
#undef BLK2
#undef BLK3
#undef FSTR0
#undef FSTR1
#undef FSTR2
#undef FSTR3
#undef FSTR4
/* procedure 37 */
int MMKernCanHandleCase(ATL_mmnode_t *mp, int mb, int nb, int kb)
/*
 * RETURNS: 0 if mp cannot handle GEMM of size mbxnbxkb
 *          1 if it can handle w/o extra computation
 *          2 if handling it requires extra comptutation (KVEC only: on kb)
 */
{
   if (mb%mp->mu == 0 && nb%mp->nu == 0)
   {
      const int KRUN=FLAG_IS_SET(mp->flag, MMF_KRUNTIME), ku=mp->ku;
      if (mp->kbB == kb)
         return(1);
      if (KRUN && kb%ku == 0)
         return(1);
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      {
         const int vl = mp->vlen, kbB=mp->kbB;
         if (KRUN && kb%ku < vl)
            return(2);
         if (kbB > kb && (kbB-kb < vl))
            return(2);
         if (kbB < kb && (kb-kbB < vl))
            return(2);
      }
   }
   return(0);
}

/* procedure 38 */
ATL_mmnode_t *MMTimeKernWithAlt
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int imf,                     /* which mflop entry to store res in */
   int flag,                    /* 1: ignore any prior output file */
   ATL_mmnode_t *mp0,           /* ptr to mmkern struct */
   ATL_mmnode_t *mp1,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
{
   int i;
   double mf0=0.0, mf1=0.0;
   i = MMKernCanHandleCase(mp0, mb, nb, kb);
   if (i)
   {
      mf0 = TimeMMKernel(verb, flag, mp0, pre, mb, nb, kb, beta, mflop, cflush);
      if (i == 2)
      {
         int k = mp0->vlen;
         k = ((kb+k-1)/k)*k;
         mf0 *= kb;
         mf0 /= k;
      }
      mp0->mflop[imf] = mf0;
   }
/*
 * If we couldn't use first kernel, or if first kernel needed extra flops
 */
   if (!i || i == 2)
   {
      int j;
      j = MMKernCanHandleCase(mp1, mb, nb, kb);
      if (j)
      {
         mf1 = TimeMMKernel(verb, flag, mp1, pre, mb, nb, kb,beta,mflop,cflush);
         if (j == 2)
         {
            int k = mp1->vlen;
            k = ((kb+k-1)/k)*k;
            mf1 *= kb;
            mf1 /= k;
         }
         mp1->mflop[imf] = mf1;
      }
      else if (!i)
         return(NULL);
   }
   if (mf1 > mf0)
      return(mp1);
   return(mp0);
}

/* procedure 39 */
ATL_mmnode_t *MMTimeAllKernsWithAlt
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int imf,
   int flag,                   /* 1: ignore any prior output file */
   ATL_mmnode_t *mmb,           /* ptr to mmkern struct */
   ATL_mmnode_t *mmA,           /* ptr to mmkern struct */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb, int nb, int kb,      /* dimensions to time */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
/*
 * For given problem size, times all kernels in mmb.  If any entry in mmb
 * cannot handle this exact problem (eg., unrolling mismatch), or requires
 * extra flops to handle the problem, the corresponding
 * entry in the alternate list mmA will be tried.
 * RETURNS: ptr to fastest timed kernel, or NULL if none worked.
 * NOTE: neither list is changed by this function (mflop is overwritten!),
 *       but ptr into original listis returned.
 */
{
   ATL_mmnode_t *mpB=NULL, *mp0, *mp1;
   double mfB=0.0;
   for (mp0=mmb, mp1=mmA; mp0; mp0 = mp0->next, mp1 = mp1->next)
   {
      ATL_mmnode_t *mp;
      mp = MMTimeKernWithAlt(verb, imf, flag, mp0, mp1, pre, mb, nb, kb,
                             beta, mflop, cflush);
      if (mp)
      {
         double mf;
         mf = mp->mflop[imf];
         if (mf > mfB)
         {
            mfB = mf;
            mpB = mp;
         }
      }
   }
   return(mpB);
}

/* procedure 40 */
int TimeNegMMKernels            /* RET: 0 if no retiming required */
(
   int imf,                     /* index of mflop array to check/set */
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flag,                    /* 1: ignore any prior output file */
   ATL_mmnode_t *mmb,           /* ptr to mmkern struct queue */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >0: size of cache flush, else ignored */
)
{
   ATL_mmnode_t *mp;
   int RETIME=0;
   for (mp=mmb; mp; mp = mp->next)
   {
      if (mp->mflop[imf] <= 0.0)
      {
         RETIME++;
         mp->mflop[imf] = TimeMMKernel(verb, flag, mp, pre,
                                       mp->mbB, mp->nbB, mp->kbB,
                                       beta, mflop, cflush);
         if (verb)
            printf("  %cID=%u, B(%u,%u,%u), mf=%.2f\n", pre, mp->ID, mp->mbB,
                   mp->nbB, mp->kbB, mp->mflop[imf]);
      }
   }
   return(RETIME);
}


/* procedure 41 */
static ATL_mmnode_t *TimeMMFile
(
   char pre,
   char *file,
   int imf,                    /* index of mflop array to check/set */
   int verb,                   /* 0: no output, 1 min output, 2: full output */
   int flag,                   /* flag for TimeMMKernel */
   int beta,                   /* beta to time */
   int mflop,                  /* >0: force mflop MFLOPs in each time interv */
   int cflush                  /* >0: size of cache flush, else ignored */
)
{
   ATL_mmnode_t *mmb;
   mmb = ReadMMFile(file);
   if (!mmb)
      return(NULL);
   MMFillInGenStrings(pre, mmb);
   printf("\nRETIMING %s:\n", file);
   if (TimeNegMMKernels(imf, verb, flag, mmb, pre, beta, mflop, cflush))
      WriteMMFile(file, mmb);
   printf("DONE.\n");
   return(mmb);
}

/* procedure 42 */
static ATL_mmnode_t *TimeMMFileWithPath
(
   char pre,
   char *path,
   char *file,
   int imf,                    /* index of mflop array to check/set */
   int verb,                   /* 0: no output, 1 min output, 2: full output */
   int flag,
   int beta,                   /* beta to time */
   int mflop,                  /* >0: force mflop MFLOPs in each time interv */
   int cflush                  /* >0: size of cache flush, else ignored */
)
{
   ATL_mmnode_t *mb;
   int CHNG=0;

   mb = ReadMMFileWithPath(pre, path, file);
   if (!mb)
      return(NULL);
   printf("\nRETIMING %s:\n", file);
   MMFillInGenStrings(pre, mb);
   if (TimeNegMMKernels(imf, verb, flag, mb, pre, beta, mflop, cflush))
      WriteMMFileWithPath(pre, path, file, mb);
   printf("DONE.\n");
   return(mb);
}

/* procedure 43 */
void TimeAllMMKernels
(
   int itime,                   /* index of mflop array to set */
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int FORCETIME,               /* 1: ignore any prior output file */
   ATL_mmnode_t *mmb,           /* ptr to mmkern struct queue */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int beta,                    /* beta to time */
   int mflop,                   /* >0: force mflop MFLOPs in each time interv */
   int cflush                   /* >0: size of cache flush, else ignored */
)
{
   ATL_mmnode_t *mmp;
   for (mmp=mmb; mmp; mmp = mmp->next)
      mmp->mflop[itime] = TimeMMKernel(verb, FORCETIME?1:0, mmp, pre,
                                       mmp->mbB, mmp->nbB, mmp->kbB,
                                       beta, mflop, cflush);
}

/* procedure 44 */
double TimeTSKernel
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int FORCETIME,               /* 1: ignore any prior output file */
   char pre,                    /* type/prec prefix: z,c,d,s */
   int mb,                      /* triangle is mbxmb (kb=mb) */
   int nb,                      /* NRHS */
   int mflop                    /* >0: force mflop MFLOPs in each time interv */
)
{
   int i, DOTIME=1;
   char ln[2048], resf[256];

   if (FORCETIME)
      strcpy(resf, "res/tmpout.ktim");
   else
   {
      sprintf(resf, "res/%ctrsm%dx%d_F%d.ktim", pre, mb, nb, mflop);
      DOTIME = !FileExists(resf);
   }

   if (DOTIME)
   {
      i = sprintf(ln, "make %ctrsmKtime mb=%d nb=%d outF=\" -f %s \"",
                  pre, mb, nb, resf);
      if (mflop)
         i += sprintf(ln+i, " FMF=%d", mflop);

      sprintf(ln+i, "\n");
      if (verb > 1)
         printf("SYSTEM: %s", ln);
      if (system(ln))
      {
         fprintf(stderr, "ERROR IN COMMAND: %s", ln);
         fprintf(stderr, "   PROPOSED FILENAME: %s\n", resf);
         sprintf(ln, "rm -f %s\n", resf);
         assert(!system(ln));
         exit(-1);
      }
   }
   return(*((double*)ReadResultsFile(0, 0, resf)));
}

#endif  /* end guard around atlas_mmtesttime.h */
