#ifndef ATLAS_CPTESTTIME_H
   #define ATLAS_CPTESTTIME_H

#include "atlas_cpparse.h"
#include "atlas_gentesttime.h"

/* procedure 1 */
ATL_cpnode_t *CopyApplyTypeScal2Hand(ATL_cpnode_t *cb, unsigned int flag)
/*
 * Get rid of cb not matching flag.  For survivors, change their optional
 * alpha/beta cases to only match the one in flag
 */
{
   ATL_cpnode_t *cp=cb;
   const unsigned int typ = flag & ((1<<CPF_SINGLE) | (1<<CPF_REAL));
   const unsigned int msk = (flag&(1<<CPF_CBLK)) ?
                            (CPF_ALLBET|CPF_ALLALP) : CPF_ALLALP, nmsk=~msk;
   while (cp)
   {
      ATL_cpnode_t *nxt = cp->next;
      cp->flag |= typ;
      if (flag&msk&cp->flag)
         cp->flag = (cp->flag&(nmsk))|(flag&msk);
      else
         cb = KillCPNodeFromQ(cb, cp);
      cp = nxt;
   }
   return(cb);
}

/* procedure 2 */
static int DoGenString(int verb, char *genstr)
{
   int err=0;
   if (verb > 2)
      printf("genstr='%s'\n", genstr);
   if (!genstr)
      err=1;
   else if (verb < 3) /* want to redirect output */
   {
      char *ln;
      ln = NewMergedString(genstr, " > /dev/null 2>&1");
      err = system(ln);
      free(ln);
   }
   else
      err = system(genstr);
   if (err)
      fprintf(stderr, "UNABLE TO GENERATE WITH COMMAND: '%s'\n",
              genstr ? genstr : "NULL");
   return(err);
}

/* procedure 3 */
double TimeCPKernel
(
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int flg,                     /* 1: ignore any prior output file */
   ATL_cpnode_t *cp,            /* ptr to cpkern struct */
   int mb, int nb,              /* rows & cols to time */
   int alpha,                   /* alpha to time [1,-1,X] */
   int beta,                    /* beta to time [0,1,-1,X] */
   unsigned long nrep,          /* >0: # of repetitions to force */
   int cflush                   /* >=0: size of cache flush, else ignored */
)
{
   char *fnam, ln[4096];  /* security from 1991 */
   const unsigned int flag=cp->flag;
   unsigned int i;
   const int FORCETIME = flg&1;
   double *dp;
/*
 * Figure out the name of the output file
 */
   if (FORCETIME)
      fnam = DupString("res/tmpout.ktim");
   else
   {
      int k;
      k = NumDecDigits(mb) + NumDecDigits(nb) + NumDecDigits(cp->vlen) + 3 + 5;
      fnam = GetCopyName(cp, k);
      fnam[0] = 'r';
      fnam[1] = 'e';
      fnam[2] = 's';
      fnam[3] = '/';
      for (k=4; fnam[k]; k++);
      if (cp->ID)
         sprintf(fnam+k, "_%ux%u_%u.ktim", mb, nb, cp->ID);
      else
         sprintf(fnam+k, "_%ux%ux%u.ktim", mb, nb, cp->vlen);
      dp = ReadResultsFile(2, 0, fnam);
      if (dp)
      {
         free(fnam);
         return(*dp);
      }
   }
   if (cp->genstr)
   {
      sprintf(ln, "%s > /dev/null 2>&1", cp->genstr);
      if (system(ln))
      {
         fprintf(stderr, "UNABLE TO GENERATE WITH='%s'\n", ln);
         assert(0);
      }
   }
   else
      assert(cp->ID);
   if (flag & (1<<CPF_CBLK))
      i = sprintf(ln, "make %ccpytimeC vlen=%d mu=%d nu=%d",
                  CopyGetPre(flag), cp->kvec?cp->kvec:1, cp->mu, cp->nu);
   else
      i = sprintf(ln, "make %ccpytime TA=%c MTDX='-D 8000 8000 8000 1'",
                  CopyGetPre(flag), CopyGetTrans(flag));
   i += sprintf(ln+i, " kfnam=%s outF=\"-f %s\"", cp->rout, fnam);
   i += sprintf(ln+i, " mb=%d nb=%d TOBLK=%u", mb, nb, (flag>>CPF_TOBLK)&1);
   if (beta == 0 || beta == 1)
      i += sprintf(ln+i, " beta=%d", beta);
   else if (beta == -1)
      i += sprintf(ln+i, " beta=-1 betan=\"N1\"");
   else
      i += sprintf(ln+i, " beta=2 betan=X");
   if (alpha == 0 || alpha == 1)
      i += sprintf(ln+i, " alpha=%d", alpha);
   else if (alpha == -1)
      i += sprintf(ln+i, " alpha=-1 alphan=\"N1\"");
   else
      i += sprintf(ln+i, " alpha=2 alphan=X");
   if (nrep)
      i += sprintf(ln+i, " FMFS=\"-r %lu\"", nrep);
   if (verb < 3)
      i += sprintf(ln+i, " > /dev/null 2>&1");
   if (system(ln))
   {
      fprintf(stderr, "ERROR IN COMMAND: '%s'\n", ln);
      if (cp->genstr)
         fprintf(stderr, "   GENSTR='%s'\n", cp->genstr);
      sprintf(ln, "rm -f %s\n", fnam);
      assert(!system(ln));
      exit(-1);
   }
   dp = ReadResultsFile(2, 0, fnam);
   free(fnam);
   return(*dp);
}

/* procedure 4 */
static int TimeNegCPKernels     /* RET: 0 if no retiming required */
(
   int imf,                     /* index of mflop array to check/set */
   int verb,                    /* 0: no output, 1 min output, 2: full output */
   int FRCTM,                   /* 1: ignore any prior output file */
   ATL_cpnode_t *cb             /* queue to time */
)
{
   ATL_cpnode_t *cp;
   int ntim=0;
   for (cp=cb; cp; cp = cp->next)
   {
      if (cp->mflop[imf] <= 0.0)
      {
         int ialp, ibet=0;
         ntim++;
         ialp = CopyGetAlphaI(cp->flag);
         if (cp->flag & (1<<CPF_CBLK))
            ibet = CopyGetBetaI(cp->flag);
         cp->mflop[imf] = TimeCPKernel(verb, FRCTM?1:0, cp, cp->mb, cp->nb,
                                       ialp, ibet, 0, -1);
      }
   }
   return(ntim);
}

/* procedure 5 */
static ATL_cpnode_t *TimeCPFile
(
   char pre,
   char *file,
   int imf,
   int verb,
   int FRCTIM   /* 1: ignore prior output file */
)
{
   ATL_cpnode_t *cb;
   cb = ReadCPFile(file);
   if (!cb)
      return(cb);
   if (TimeNegCPKernels(imf, verb, FRCTIM, cb))
      WriteCPFile(file, cb);
   return(cb);
}

/* procedure 5 */
static ATL_cpnode_t *TimeCPFileWithPath
(
   char pre,
   char *path,
   char *file,
   int imf,
   int verb,
   int FRCTIM   /* 1: ignore prior output file */
)
{
   ATL_cpnode_t *cb;
   cb = ReadCPFileWithPath(pre, path, file);
   if (!cb)
      return(cb);
   if (TimeNegCPKernels(imf, verb, FRCTIM, cb))
      WriteCPFileWithPath(pre, path, file, cb);
   return(cb);
}

/* procedure 6 */
int CPKernelFailsTestMM(char pre, int mb, int nb, int ialp, int ibet,
                        ATL_cpnode_t *cp)
/*
 * This tester only works for some alpha/beta combos.
 * RETURNS: 0 on success, non-zero on failure
 */
{
   int i, sz, mu=cp->mu, nu=cp->nu, ku=(cp->kvec)?cp->kvec:1, kb;
   char *ln, *tst="ammmcpytst";

   kb = (mb >= nb) ? mb : nb;
   kb = ((kb+ku-1)/ku)*ku;
   if (!(cp->flag & (1<<CPF_AL1)))  /* this tester can only test alpha=1 */
      return(0);                    /* so say anything else passes for now */
   ialp = 1;
   if (cp->flag & (1<<CPF_SYRK))
   {
      tst = "syrkcpytst";
      if (cp->flag & (1<<CPF_CBLK))
         assert(cp->mu == cp->nu);
   }
   if (cp->flag & (1<<CPF_CBLK))
   {  /* can only test alpha=1, beta=0 for these guys */
      if (!(cp->flag & (1<<CPF_BE0)))  /* this tester can only test beta=0 */
         return(0);                    /* so say anything else passes for now */
      ibet=0;
   }
   else /* can only test alpha=1 & TOBLK for A/B */
   {
      ibet = 0;  /* not used for this test */
      mu = nu;   /* so that A & B are same */
      mb = ((mb+mu-1) / mu)*mu;
      if (!(cp->flag & (1<<CPF_TOBLK)))
         return(0);
   }
   if (cp->genstr)
      assert(!DoGenString(0, cp->genstr));
   sz = 10 + 62 + NumDecDigits(nb) + NumDecDigits(mb) + NumDecDigits(kb);
   sz += NumDecDigits(mu) + NumDecDigits(nu) + NumDecDigits(ku);
   if (cp->kvec > 1)
      sz += 21 + 2*NumDecDigits(cp->kvec);
   else
      sz += 14;
   sz += 7 + strlen(cp->rout);
   sz += 17;
   ln = malloc(sz);
   assert(ln);
   i = sprintf(ln,
               "make %c%s alpha=%d beta=%d mu=%u nu=%u ku=%u mb=%u nb=%u kb=%u",
               pre, tst, ialp, ibet, mu, nu, ku, mb, nb, kb);
   if (cp->kvec > 1)
      i += sprintf(ln+i, " vec=kdim vlen=%u kmaj=%u", cp->kvec, cp->kvec);
   else
      i += sprintf(ln+i, " vec=no vlen=1");
   if (cp->flag & (1<<CPF_CBLK))
   {
      if (cp->flag & (1<<CPF_TOBLK))
         i += sprintf(ln+i, " C2blk=%s", cp->rout);
      else
         i += sprintf(ln+i, " blk2C=%s", cp->rout);
   }
   else
   {
      if (cp->flag & (1<<CPF_TRANS))
         i += sprintf(ln+i, " RM2blk=%s", cp->rout);
      else
         i += sprintf(ln+i, " CM2blk=%s", cp->rout);
   }
   i += sprintf(ln+i, " > /dev/null 2>&1");
   assert(i < sz);
   i = system(ln);
   if (i)
   {
      fprintf(stderr, "FAILED COMMAND: '%s'\n", ln);
      if (cp->genstr)
         fprintf(stderr, "   genstr= '%s'\n", cp->genstr);
   }
   free(ln);
   return(i);
}

/* procedure 7 */
static int CPKernelFailsTest
(
   char pre,    /* s,d,c,z */
   int TA,      /* 0: copy to/from col-major storage, else row-major */
   int mb,      /* if TA, # of cols, else # of rows */
   int nb,      /* if TA, # of rows, else # of cols */
   int ialp,    /* 0,1,-1,2 */
   int ibet,    /* 0,1,-1,2, unused if not C copy */
   ATL_cpnode_t *cp
)
/*
 * RETURNS: 0 on success, non-zero on failure
 */
{
   const unsigned int flag=cp->flag;
   int i, ierr, L, mu=cp->mu, nu=cp->nu, ku=(cp->kvec)?cp->kvec:1, kb, CPLX=0;
   unsigned int sz;
   char *ln;

   TA = (TA) ? 1 : 0;  /* make sure 0 or 1 only */
   L = NumDecDigits(mb)+NumDecDigits(nb)+NumDecDigits(cp->mu)+NumDecDigits(nu)
       + NumDecDigits(mb*nb*2)+NumDecDigits(cp->kvec)+NumDecDigits(cp->vlen);
   if (pre == 'z' || pre == 'c')
   {
      CPLX = 0;
      L += 14; /* " CNJ='-DConj_'" */
   }
   L += 52 + 54 + 10 + 1 + strlen(cp->rout);
   ln = malloc(L);
   assert(ln);
   if (cp->genstr)
      assert(!DoGenString(0, cp->genstr));
   if (flag & (1<<CPF_CBLK))
   {
      assert(TA == 0);
      if (cp->kvec > 1)
      {
         int vl = cp->kvec;
         sz = ((mu*nu+vl-1)/vl)*vl;
      }
      else
         sz = mu*nu;
      sz *= ((mb+mu-1)/mu)*((nb+nu-1)/nu)*sz;
      i = sprintf(ln, "make %ccopytest CPM=C beta=%d alpha=%d UR=%u ku=%u "
                  "vlen=%u mb=%u nb=%u szb=%u TO_BLK=%u TRANS=0 kvec=%u "
                  "kfnam='%s'",
                  pre, ibet, ialp, mu, nu, cp->vlen, mb, nb, sz,
                  ((flag>>CPF_TOBLK)&1), (cp->kvec > 1)?cp->kvec:0, cp->rout);
   }
   else /* A or B matrix */
   {
      sz = ((mb+nu-1)/nu)*((nb+mu-1)/mu)*mu*nu;
      i = sprintf(ln, "make %ccopytest CPM=A beta=%d alpha=%d UR=%u ku=%u "
                  " vlen=%u mb=%u nb=%u szb=%u TO_BLK=%u TRANS=%u kvec=%u "
                  "kfnam='%s'",
                  pre, ibet, ialp, nu, mu, cp->vlen, mb, nb, sz,
                  ((flag>>CPF_TOBLK)&1), TA, (cp->kvec > 1)?cp->kvec:0,
                  cp->rout);
   }
   assert(i < L);
   ierr = Sys2File(ln, "ATL_tmp.out");
   if (ierr)
   {
      #if 0
         fprintf(stderr, "\nOUTPUT OF FAILED COMMAND:\n"
         system("cat ATL_tmp.out");
      #endif
      fprintf(stderr, "\nFAILED (%d) COMMAND: '%s'\n", i, ln);
      if (cp->genstr)
         fprintf(stderr, "   genstr= '%s'\n", cp->genstr);
   }
   remove("ATL_tmp.out");
   if (!ierr && CPLX)
   {
      assert(i+14 < L);
      strcpy(ln+i, " CNJ='-DConj_'");
      ierr = Sys2File(ln, "ATL_tmp.out");
      if (ierr)
      {
         #if 0
            fprintf(stderr, "\nOUTPUT OF FAILED COMMAND:\n"
            system("cat ATL_tmp.out");
         #endif
         fprintf(stderr, "\nFAILED CNJ (%d) COMMAND: '%s'\n", i, ln);
         if (cp->genstr)
         fprintf(stderr, "   genstr= '%s'\n", cp->genstr);
      }
      remove("ATL_tmp.out");
   }
   free(ln);
   return(ierr);
}

/* procedure 8 */
ATL_cpnode_t *KillFailingCPNodes(char pre, ATL_cpnode_t *cb)
{
   ATL_cpnode_t *cp;
   int i;
   int nerr=0;
   for (i=0,cp=cb; cp; i++)
   {
      ATL_cpnode_t *next = cp->next;
      int err, ibet=0, ialp, mb=cp->mb, nb=cp->nb;

      ialp = CopyGetAlphaI(cp->flag);
      if (cp->flag & (1<<CPF_CBLK))
         ibet = CopyGetBetaI(cp->flag);
      if (!mb)
         mb = 5*cp->mu;
      if (!nb)
         nb = 7*cp->nu;
      if (cp->flag&(1<<CPF_SYRK))
         nb = mb;
      err = CPKernelFailsTest(pre, ((cp->flag>>CPF_TRANS)&1), mb, nb,
                              ialp, ibet, cp);
      if (cp->ID)
         printf("%3d. ID=%d (%s) : %s\n", i, cp->ID, cp->rout,
                (err) ? "FAIL!":"PASS.");
      else
         printf("%3d. B=(%u,%u), U=(%u,%u), VL=(%d,%d) : %s\n", i,
                cp->mb, cp->nb, cp->mu, cp->nu, cp->vlen, cp->kvec,
                (err) ? "FAIL!":"PASS.");
      if (err)
      {
         cb = KillCPNodeFromQ(cb, cp);
         nerr++;
      }
      cp = next;
   }
   if (nerr)
      printf("FAILED %u of %u tests!\n\n", nerr, i);
   else
      printf("PASSED ALL %u tests.\n\n", i);
   return(cb);
}

#endif
