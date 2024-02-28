#include "atlas_cpparse.h"
#include "atlas_cptesttime.h"
#include "atlas_mmparse.h"
#include "atlas_mmtesttime.h"
#include "atlas_mmgen.h"

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -u <ucase.idx>: (none) user cases to try\n");
   fprintf(stderr, "   -i <cpyin.CPS>: (stdin) list of copy to search\n");
   fprintf(stderr, "   -o <cpyout.CPS>: (stdout) fastest copy kernels\n");
   exit(ierr ? ierr : -1);
}

ATL_cpnode_t *GetFlags(int nargs, char **args, int *VERB, char **FOUT,
                       ATL_cpnode_t **UB)
/*
 * RETURNS: list of all kerns/blk factors to tune copies for
 */
{
   ATL_cpnode_t *cb=NULL, *cp, *ub=NULL;
   int i;

   *VERB = 0;
   *FOUT = NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'u':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        ub = ReadCPFile(args[i]);
        break;
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        cb = ReadCPFile(args[i]);
        assert(cb);
        break;
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *FOUT = args[i];
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!cb)
      cb = ReadCPFile(NULL);
   *UB = ub;
   return(cb);
}

unsigned long GetRep(unsigned long reps, unsigned int mb, unsigned int nb)
{
   double ratio;
   unsigned long nrep;
   ratio = 16.0 / (((double)mb)*nb);
   nrep = reps*ratio;
   nrep = (nrep > 0) ? nrep : 1;
   return(nrep);
}

ATL_cpnode_t *GetWorkingKerns(char pre, int verb, ATL_cpnode_t *cb)
/*
 * Right now, just reads in all kerns.  Later, should eliminate kerns that
 * won't work with any cb.
 */
{
   ATL_cpnode_t *hb=NULL;
   int DOTEST=0;
   if (!hb)
   {
      hb = ReadCPFileWithPath(pre, "CPYCASES", "copy.idx");
      #if 0
      if (flag&(1<<CPF_CBLK))
         DOTEST = (CopyGetAlphaI(flag) == 1 && CopyGetBetaI(flag) == 0);
      else
         DOTEST = CopyGetAlphaI(flag) == 1;
      #endif
   }

   PrefixStrAllNodes(hb, GetOffset(&hb->next, hb), GetOffset(&hb->rout, hb),
                     "CPYCASES/");
   if (DOTEST)
   {
      hb = KillFailingCPNodes(pre, hb);
      WriteCPFileWithPath(pre, "res", "cpywrk_a1b0.sum", hb);
   }
   return(hb);
}

ATL_cpnode_t *DoHandTimings
   (char pre, int verb, unsigned long rep4, ATL_cpnode_t *cb, ATL_cpnode_t *hb)
{
   ATL_cpnode_t *cp, *hp, *HB=NULL;
   unsigned long nrep;


   if (!hb)
   {
      printf("NO USER-SUPPLIED COPY ROUTINES!\n");
      return(NULL);
   }
   printf("FINDING BEST HAND-TUNED COPY:\n");
   for (cp=cb; cp; cp=cp->next)
   {
      double mf;
      ATL_cpnode_t *np=NULL;
      const unsigned int mb=cp->mb, nb=cp->nb, flag=cp->flag;
      int ialp, ibet=0;

      nrep = GetRep(rep4, mb, nb);
      ialp = CopyGetAlphaI(flag);
      if (flag & (1<<CPF_CBLK))
         ibet = CopyGetBetaI(flag);
/*
 *    Try any hand-tuned cases that work
 */
      hp = FindEquivUserCopy(hb, cp);
      if (hp)
      {
         printf("   TIMINGS B=(%d,%d), U=(%d,%d) TA=%c, reps=%lu:\n",
                cp->mb, cp->nb, cp->mu, cp->nu,
                (cp->flag&(1<<CPF_TRANS))?'T':'N', nrep);
         printf("      ID=%6d VL=%2d: mf=%.0f\n", 0, cp->vlen, cp->mflop[0]);
         do
         {
            double mf;

            mf = TimeCPKernel(verb, 0, hp, mb, nb, ialp, ibet, nrep, -1);
            printf("      ID=%6d VL=%2d: mf=%.0f\n", hp->ID, hp->vlen, mf);
            if (mf > cp->mflop[0])
            {
               if (np)
                  KillCPNode(np);
               np = CloneCPNode(hp);
               np->mflop[0] = mf;
               np->mb = mb;
               np->nb = nb;
               np->flag = CopyEncodeScal(flag, ialp, ibet);
            }
            hp = hp->next;
         }
         while( (hp = FindEquivUserCopy(hp, cp)) );
         printf("   DONE\n");
      }
      if (!np)
         np = CloneCPNode(cp);
      np->next = HB;
      HB = np;
   }
   if (hb)
      KillAllCPNodes(hb);
   return(ReverseCPQ(HB));
}

void DoGenTimings(char pre, int verb, unsigned long rep4, ATL_cpnode_t *cb)
{
   ATL_cpnode_t *cp;
   char *fnam;

   printf("FINDING BEST GENERATED COPY:\n");
   for (cp=cb; cp; cp = cp->next)
   {
      double mf, mfB;
      unsigned long nrep;
      int ialp, ibet=0, flag=cp->flag;
      ialp = CopyGetAlphaI(flag);
      if (flag & (1<<CPF_CBLK))
         ibet = CopyGetBetaI(flag);

      nrep = GetRep(rep4, cp->mb, cp->nb);
      printf("   TIMINGS B=(%d,%d), U=(%d,%d), TA=%c, reps=%lu:\n",
             cp->mb, cp->nb, cp->mu, cp->nu, (cp->flag&(1<<CPF_TRANS))?'T':'N',
             nrep);
      cp->vlen = 1;
      if (!cp->rout)
         cp->rout = DupString("ATL_tmp.c");
      if (cp->genstr)
         free(cp->genstr);
      cp->genstr = GetCopyGenStr(cp);
      mfB = TimeCPKernel(verb, 0, cp, cp->mb, cp->nb, ialp, ibet, nrep, -1);
      printf("      ID=     0 VL= 1: mf=%.0f\n", mfB);
      cp->mflop[0] = mfB;
/*
 *    See if we can use vectorized generator.
 *    Present vector generator only supports SSE/AVX for real GEMM C copy
 */
      #if defined(ATL_AVX) || defined(ATL_SSE2) || defined(ATL_SSE)
         if (!(cp->flag & CPF_ALLKERN) && (pre == 's' || pre == 'd') &&
             (cp->flag&(1<<CPF_CBLK)))
         {
            char *gens;
            int VL0=cp->vlen, VL = (pre == 's') ? 4 : 2;

         #if defined(ATL_SSE2) || defined(ATL_SSE1)
            #if defined(ATL_SSE2)
            if (cp->mu % VL == 0 && (pre == 'd' || pre == 's') &&
                !(cp->flag&CPF_ALLKERN))
            #else
            if (cp->mu % VL == 0 && pre == 's' && !(cp->flag&CPF_ALLKERN))
            #endif
            {
               gens = cp->genstr;
               cp->vlen = VL;
               cp->genstr = GetCopyGenStr(cp);
               mf = TimeCPKernel(verb, 0, cp, cp->mb, cp->nb, ialp, ibet,
                                 nrep, -1);
               printf("      ID=     0 VL=%2d: mf=%.0f\n", VL, mf);
               if (mf > mfB)
               {
                  free(gens);
                  cp->mflop[0] = mfB = mf;
               }
               else
               {
                  free(cp->genstr);
                  cp->genstr = gens;
                  cp->vlen = VL0;
               }
            } /* end SSE if */
         #endif
         #ifdef ATL_AVX
            VL0 = cp->vlen;
            VL <<= 1;
            if (cp->mu % VL == 0)
            {
               gens = cp->genstr;
               cp->vlen = VL;
               cp->genstr = GetCopyGenStr(cp);
               mf = TimeCPKernel(verb, 0, cp, cp->mb, cp->nb, ialp, ibet,
                                 nrep, -1);
               printf("      ID=     0 VL=%2d: mf=%.0f\n", VL, mf);
               if (mf > mfB)
               {
                  free(gens);
                  cp->mflop[0] = mfB = mf;
               }
               else
               {
                  free(cp->genstr);
                  cp->genstr = gens;
                  cp->vlen = VL0;
               }
            } /* end AVX if */
         #endif
         }    /* end if over trying old vec-basefile */
      #endif  /* end #if on Intel SSE/AVX */
   }
   printf("DONE GENERATED COPY SEARCH\n");
}

void DoTimings
   (char *fnam, char pre, int verb, unsigned long rep4, ATL_cpnode_t *cb,
    ATL_cpnode_t *ub)
{
   ATL_cpnode_t *cp, *fb;

   cp = TimeCPFile(pre, fnam, 0, verb, 0);
   if (cp)
   {
      KillAllCPNodes(cp);
      return;
   }
   DoGenTimings(pre, verb, rep4, cb);
   fb = DoHandTimings(pre, verb, rep4, cb, ub);
   if (!fb)
      fb = cb;
   for (cp=fb; cp; cp=cp->next)
   {
      if (cp->flag & (1L<<CPF_SYRK))
         cp->mflop[0] = 2e-6 / cp->mflop[0]; /* time per elt */
      else
         cp->mflop[0] = 1e-6 / cp->mflop[0]; /* time per elt */
   }
   WriteCPFile(fnam, fb);
   if (fb != cb)
      KillAllCPNodes(fb);
   printf("DONE, OUTPUT IN: %s.\n", fnam);
}

unsigned long GetReps4x4(int verb, double res)
/*
 * Find how many repititions to do for 4x4 copy to get timing interval to
 * around res seconds
 */
{
   ATL_cpnode_t *cp;
   double tim, mf;
   unsigned long nrep = 32;
   FILE *fp;

   fp = fopen("res/ncprep4x4.txt", "r");
   if (fp)
   {
      assert(fscanf(fp, " %lu", &nrep) == 1);
      fclose(fp);
      return(nrep);
   }
   cp = GetCPNode();
   cp->mu = 4;
   cp->nu = 1;
   cp->flag = (1<<CPF_BE0)|(1<<CPF_AL1)|(1<<CPF_CBLK)|
              (1<<CPF_SINGLE)|(1<<CPF_REAL);
   cp->rout = DupString("ATL_tmp.c");
   cp->genstr = GetCopyGenStr(cp);

   printf("FINDING NREP:\n");
   do
   {
      nrep += nrep;
      mf = TimeCPKernel(verb, 1, cp, 4, 4, 1, 0, nrep, -1);
      tim = (mf / 16e6) * nrep;
      printf("   NREP=%16lu, tim=%e\n", nrep, tim);
   }
   while (tim < res && nrep < 100000);
   printf("DONE.\n");
   KillCPNode(cp);
   fp = fopen("res/ncprep4x4.txt", "w");
   fprintf(fp, "%lu", nrep);
   fclose(fp);
   return(nrep);
}

int main(int nargs, char **args)
/*
 * For now read geAMMRES.sum, later add rkAMMRES.sum & SYRK.
 * Will use this info to find unique lists of AB & C copiers.
 * -> Later must extend for A & B using possibly different copiers.
 * It then finds the fastest working copy routine for each case,
 * trying each of following:
 * (1) Codes generated from atlas-mmkg.base (all standard formats)
 * (2) atlas-mmg.base (vector ops, but only mu%vlen Cblk formats)
 * (3) User-contributed cases
 */
{
   ATL_cpnode_t *cb, *cp, *cb0, *ub;
   char *fout;
   unsigned long nrep;
   char pre;
   int verb, minSz, ialp, ibet=0;
   unsigned int flag;

   cb = GetFlags(nargs, args, &verb, &fout, &ub);
   pre = CopyGetPre(cb->flag);
   nrep = GetReps4x4(verb, 0.10);
   DoTimings(fout, pre, verb, nrep, cb, ub);

   KillAllCPNodes(cb);
   return(0);
}
