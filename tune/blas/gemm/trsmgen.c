/*
 * Actually generates TRSM micro & nano kernels and their Makefile, and
 * moves this to BLDdir/src/blas/ukernel/<pre>UTRSM
 * Actual performance data is handled by views.
 */
#include "atlas_mmgen.h"
#include "atlas_type.h"
ATL_mmnode_t *allT2bv(ATL_mmnode_t *mp, unsigned long *lp)
{
   unsigned long bv=0;
   unsigned int n;
   assert(sizeof(long) >= ATL_PSIZE);
   for (n=0; mp && n < (ATL_PSIZE<<3); n++, mp=mp->next)
      if (mp->TA == AtlasTrans && mp->TB == AtlasTrans)
         bv |= (1L<<n);

   *lp = bv;
   return(mp);
}
void GenTrsmALLT
   (char pre, char sd, char up, char ta, FILE *fp, ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;
   unsigned long flg=0;
   int k, L, n, ALL1=1, ALL0=1;

   for (mp=mb, n=0; mp; mp = mp->next, n++)
   {
      if (mp->TA == AtlasTrans && mp->TB == AtlasTrans)
         ALL0 = 0;
      else
         ALL1 = 0;
   }
   if (ALL0)
      fprintf(fp, "#define ATL_trsm%c_%c%c_allT(i_) 0\n", sd, up, ta);
   else if (ALL1)
      fprintf(fp, "#define ATL_trsm%c_%c%c_allT(i_) 1\n", sd, up, ta);
   else
   {
      mp = allT2bv(mb, &flg);
      if (!mp) /* only 1 constant needed */
         fprintf(fp, "#define ATL_trsm%c_%c%c_allT(i_) ((0x%lxL>>(i_))&1)\n",
                 sd, up, ta, flg);
      else
      {
         int ne = (ATL_PSIZE<<3);
         ATL_mmnode_t *cp;
         for (cp=mb, n=0; cp; cp = cp->next, n++);
         ne = (n+ne-1)/ne;
         fprintf(fp, "#ifdef ATL_DECL_\n   const ATL_iptr_t ATL_UTRSM%c_%c%c"
                 "_ALLT[%u]={0x%lxL", sd, up, ta, ne, flg);
         do
         {
            mp = allT2bv(mp, &flg);
            fprintf(fp, ",0x%lxL", flg);
         }
         while(mp);
         fprintf(fp, "};\n#else\n   extern const ATL_iptr_t ATL_UTRSM%c_%c%c"
                 "_ALLT[%u];\n#endif\n", sd, up, ta, ne);
         for (ne=(ATL_PSIZE<<3), n=0; (ne^(1L<<n)); n++);
         n++;
         fprintf(fp, "#define ATL_trsm%c_%c%c_allT(i_) \\\n", sd, up, ta);
         fprintf(fp, "   ((ATL_UTRSM%c_%c%c_ALLT[((i_)>>%u)]>>((i_)-"
                 "(((i_)>>%u)<<%u)))&1)\n", sd, up, ta, n, n, n);
      }
   }
   fprintf(fp, "#ifndef ATL_trsm_allT\n");
   fprintf(fp, "   #define ATL_trsm_allT ATL_trsm%c_%c%c_allT\n",
           sd, up, ta);
   fprintf(fp, "#endif\n");
}

char *GetTrsmDecor(char pre, char sd, char up, char ta, ATL_mmnode_t *mp)
{
   static char fn[64];
   char ALLT = (mp->TA == AtlasTrans && mp->TB == AtlasTrans) ? 'T':'N';
   if (ALLT == 'T')
      sprintf(fn, "%c_%c%c%ux%u%c", sd, up, ta, mp->nu, mp->mu, ALLT);
   else
      sprintf(fn, "%c_%c%c%ux%u%c", sd, up, ta, mp->mu, mp->nu, ALLT);
   return(fn);
}

void genTrsmHead(char pre, char *outd, char sd, char up, char ta,
                 ATL_mmnode_t *mb, ATL_mmnode_t *ub)
{
   char *of;
   FILE *fp;
   ATL_mmnode_t *mp;
   int k, L, n, muoff, nuoff, nxtoff, MULMU;
   L = strlen(outd) + 24;
   of = malloc(L);
   assert(of);

   k = sprintf(of, "%s/atlas_%cutrsm%c_%c%c.h", outd, pre, sd, up, ta);
   assert(k<L);
   fp = fopen(of, "w");
   assert(fp);
   fprintf(fp, "#ifndef ATLAS_%cTRSM%c_%c%c_H\n",
           toupper(pre), toupper(sd), toupper(up), toupper(ta));
   fprintf(fp, "   #define ATLAS_%cTRSM%c_%c%c_H 1\n\n",
           toupper(pre), toupper(sd), toupper(up), toupper(ta));

   fprintf(fp, "#include \"atlas_amm.h\"\n");
   GenTrsmALLT(pre, sd, up, ta, fp, mb);
   fprintf(fp, "/*\n * TRSM microkernel prototypes\n */\n");
   for (mp=ub; mp; mp = mp->next)
   {
      char *nm;
      nm = GetTrsmDecor(pre, sd, up, ta, mp);
      fprintf(fp, "void ATL_%cutrsm%s\n", pre, nm);
      fprintf(fp, "   (sminfo_t *ip, const enum ATLAS_DIAG Diag,\n"
              "    ATL_CINT N, ATL_CINT R, const SCALAR alpha, const TYPE *A,"
              " ATL_CINT lda,\n    TYPE *X, ATL_CINT ldx, TYPE *diag, "
              "TYPE *L, TYPE *RW, TYPE *w);\n");
   }
   fprintf(fp, "\n#ifndef ATL_NOLOOKUP\n");
   fprintf(fp, "#ifndef ATL_findutrsm\n");
   fprintf(fp, "   #define ATL_findutrsm ATL_findutrsm%c_%c%c\n", sd, up, ta);
   fprintf(fp, "#endif\n");
   fprintf(fp, "static INLINE void *ATL_findutrsm%c_%c%c(ATL_UINT ALLT, "
           "ATL_UINT mu, ATL_UINT nu)\n", sd, up, ta);
   fprintf(fp, "{\n");
   fprintf(fp, "   void *vp=NULL;\n");
   fprintf(fp, "\n");
   ub = CloneMMQueue(ub);
   muoff = GetOffset(&ub->mu, ub);
   nuoff = GetOffset(&ub->nu, ub);
   nxtoff = GetOffset(&mb->next, mb);
   ub = SortListByIval_G2L(ub, nxtoff, muoff);
   k = ub->mu;
   for (mp=ub->next; mp && mp->mu == k; mp = mp->next);
   MULMU = (mp != NULL);
   if (MULMU)
      fprintf(fp, "   switch(mu)\n   {\n");
   while(ub)
   {
      char *nm;
      int mu=ub->mu, nu=ub->nu;
      int n;
      ATL_mmnode_t *mub, *last;

      if (MULMU)
         fprintf(fp, "   case %u:\n", mu);
/*
 *    Remove all identical mu kernels from queue
 */
      last = mub = ub;
      for (n=1,mp=ub->next; mp && mp->mu == mu; mp = mp->next, n++)
         last = mp;
      ub = mp;
      last->next = NULL;
/*
 *    Find if there's another utrsm with same mu & nu, then one must be
 *    using transpose case while other does not
 */
      mp = FindNodeWithIval(mub->next, nxtoff, nuoff, mub->nu);
      if (mp)
      {
         mub = RemoveNodeFromList(mub, mp, nxtoff);
         n--;
      }
      nm = GetTrsmDecor(pre, sd, up, ta, mub);
      if (n == 1)  /* is there only 1 kernel with this mu? */
      {
         if (mp)
         {
            int ALLT = (mub->TA == AtlasTrans && mub->TB == AtlasTrans);
            fprintf(fp, "         vp = (%s) ? ATL_%cutrsm%s :",
                    ALLT ? "ALLT":"!ALLT", pre, nm);
            nm = GetTrsmDecor(pre, sd, up, ta, mp);
            fprintf(fp, " ATL_%cutrsm%s;\n", pre, nm);
            KillMMNode(mp);
         }
         else
            fprintf(fp, "         vp = ATL_%cutrsm%s;\n", pre, nm);
         mub = KillMMNode(mub);
      }
      else
      {
         n = CountListEntries(mub, nxtoff);
         if (n > 1)
            fprintf(fp, "      switch(nu)\n      {\n");
         while(mub)
         {
            if (n > 1)
               fprintf(fp, "      case %u:\n", mub->nu);
            if (mp)
            {
               int ALLT = (mub->TA == AtlasTrans && mub->TB == AtlasTrans);
               fprintf(fp, "         vp = (%s) ? ATL_%cutrsm%s :",
                       ALLT ? "ALLT":"!ALLT", pre, nm);
               nm = GetTrsmDecor(pre, sd, up, ta, mp);
               fprintf(fp, " ATL_%cutrsm%s;\n", pre, nm);
               KillMMNode(mp);
            }
            else
               fprintf(fp, "         vp = ATL_%cutrsm%s;\n", pre, nm);
            if (n > 1)
               fprintf(fp,    "         break;\n");
            mub = KillMMNode(mub);
            if (mub)
            {
               nm = GetTrsmDecor(pre, sd, up, ta, mub);
               mp = FindNodeWithIval(mub->next, nxtoff, nuoff, mub->nu);
               if (mp)
                  mub = RemoveNodeFromList(mub, mp, nxtoff);
            }
         }
         if (n > 1)
            fprintf(fp, "      } /* end switch on nu */\n");
      }
      if (MULMU)
         fprintf(fp, "      break;\n");
   }
   if (MULMU)
      fprintf(fp, "   } /* end switch on mu */\n");
   fprintf(fp, "   return(vp);\n");
   fprintf(fp, "} /* end findUtrsm */\n");
   fprintf(fp, "#endif /* end no-lookup guard */\n");

   fprintf(fp, "\n#endif /* end multiple inclusion guard */\n");
   fclose(fp);
   free(of);
}

ATL_mmnode_t *getUniqueTrsmCases(char pre, ATL_mmnode_t *ub, ATL_mmnode_t *mb)
/*
 * Adds all cases of unique cases in mb to ub; mb nodes are cloned, so
 * mb is unchanged.  ub is already assumed unique on entry, so it is only
 * added to.
 * RETURNS: possibly changed ub
 */
{
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
   {
      ATL_mmnode_t *p;
      int mu=mp->mu, nu=mp->nu;
      char ALLT=(mp->TA == AtlasTrans && mp->TB == AtlasTrans) ? 'T' : 'N';
      assert(mp->TA == mp->TB);
      for (p=ub; p; p = p->next)
      {
         int mu0=p->mu, nu0=p->nu;
         char ALLT0=(p->TA == AtlasTrans && p->TB == AtlasTrans) ? 'T' : 'N';
         if (mu0 == mu && nu0 == nu && ALLT0 == ALLT) /* already seen this */
            break;                                    /* case so stop */
      }
      if (p == NULL)  /* this case not already in unique queueu */
      {
         p = CloneMMNode(mp);
         p->next = ub;
         ub = p;
      }
   }
   return(ub);
}

ATL_mmnode_t *GenMicroTrsm(char pre, char *outd, char sd, char up, char ta,
                           ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp, *ub;
   int *mus, *nus;
   char *fn;
   int i, L;
   int nmu, nnu;
   char *ln;

   L = strlen(outd) + 256;
   ln = malloc(L);
   assert(ln);
   ub = getUniqueTrsmCases(pre, NULL, mb);
   for (mp=ub; mp; mp = mp->next)
   {
      int mu=mp->mu, nu=mp->nu;
      char *nm;
      char ALLT=(mp->TA == AtlasTrans && mp->TB == AtlasTrans) ? 'T' : 'N';

      assert(mu > 0 && nu > 0 && mu < 1000 && nu < 1000);
      nm = GetTrsmDecor(pre, sd, up, ta, mp);
      i = sprintf(ln, "make gen_utrsm tALL=%c mu=%u nu=%u sd=%c up=%c ta=%c"
                  " cx=\"%s\" cnj=N rt=\"%s/ATL_utrsm%s.c\"",
                  ALLT, mu, nu, sd, up, ta, (pre == 'c' || pre == 'z')?"c":"",
                  outd, nm);
      assert(i < L);
      assert(!Sys2File(ln, NULL));
   }
   free(ln);
   return(ub);
}

void GenMakeTrsm(char pre, char *outd, char up, char *sds, char *tas,
                 ATL_mmnode_t **MBs)
{
   int k, i;
   FILE *fp;
   char *fn, *typ, *comp;

   k = strlen(outd) + 10;
   fn = malloc(k);
   assert(fn);
   i = sprintf(fn, "%s/Makefile", outd);
   fp = fopen(fn, "w");
   assert(fp);
   free(fn);
   fprintf(fp, "include ../Make.inc\n\n");

   fprintf(fp, "lib : %clib.grd\n", pre);
   fprintf(fp, "all : %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n\n", pre);

   fprintf(fp, "obj = ");
   if (pre == 'c' || pre == 's')
      comp = "$(SKC)";
   else
      comp = "$(DKC)";
   if (pre == 'z')
      typ = "DCPLX";
   else if (pre == 'c')
      typ = "SCPLX";
   else if (pre == 's')
      typ = "SREAL";
   else
      typ = "DREAL";
   for (k=0; k < 2; k++)
   {
      int s;
      for (i=s=0; s < 2; s++)
      {
         char sd=sds[s];
         int it;

         for (it=0; it < 2; it++)
         {
            ATL_mmnode_t *mp;
            char ta = tas[it];
            const int NCNJ = (ta == 'T' && (pre == 'c' || pre == 'z')) ? 2:1;

            for (mp=MBs[i++]; mp; mp=mp->next)
            {
               int c;
               for (c=0; c < NCNJ; c++)
               {
                  char *nm;
                  nm = GetTrsmDecor(pre, sd, up, ta, mp);
                  if (c)
                     nm[3] = 'H';
                  if (!k)
                     fprintf(fp, "\\\n      ATL_%cutrsm%s.o ", pre, nm);
                  else
                  {
                     fprintf(fp, "ATL_%cutrsm%s.o : ", pre, nm);
                     if (c)
                        nm[3] = 'T';
                     fprintf(fp, "ATL_utrsm%s.c $(deps)\n", nm);
                     if (c)
                        fprintf(fp,
                "\t%s $(%cKCFLAGS) -c -o $@ -D%s=1 -DConj_=1 ATL_utrsm%s.c\n",
                           comp, pre, typ, nm);
                     else
                        fprintf(fp,
                           "\t%s $(%cKCFLAGS) -c -o $@ -D%s=1 ATL_utrsm%s.c\n",
                                comp, pre, typ, nm);
                  }
               }
            }
         }
      }
      if (!k)
      {
         fprintf(fp, "\n\n%clib : %clib.grd\n", pre, pre);
         fprintf(fp, "%clib.grd : $(obj)\n", pre);
         fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(obj)\n");
         fprintf(fp, "\t$(RANLIB) $(ATLASlib)\n\ttouch %clib.grd\n\n", pre);
      }
   }
   fclose(fp);
}

void GenAllTrsm(char pre, char *outd)
{
   int L, i, s;
   char sds[2] = {'L', 'R'};
   char tas[2] = {'N', 'T'};
   char up='L';
   char *od;
   char fn[16];
   ATL_mmnode_t *MBs[4];
/*
 * Delete and the re-create <pre>UTRSM subdir
 */
   L = strlen(outd) + 8 + 7;
   od = malloc(L);
   assert(od);
   i = sprintf(od, "rm -rf %s/%cUTRSM", outd, pre);
   assert(i < L);
   Sys2File(od, NULL);
   od[0]='m'; od[1]='k'; od[2]='d'; od[3]='i'; od[4]='r'; od[5]=' ';
   assert(!Sys2File(od, NULL));
   strcpy(fn, "trsmS_LN.sum");

   for (i=s=0; s < 2; s++)
   {
      char sd=sds[s];
      int it;

      for (it=0; it < 2; it++)
      {
         char ta = tas[it];
         ATL_mmnode_t *mb;
         fn[4] = sd;
         fn[7] = ta;
         mb = ReadMMFileWithPath(pre, "res", fn);
         MBs[i] = GenMicroTrsm(pre, od+7, sd, up, ta, mb);
         genTrsmHead(pre, od+7, sd, 'L', ta, mb, MBs[i]);
         KillAllMMNodes(mb);
         i++;
      }
   }
   GenMakeTrsm(pre, od+7, up, sds, tas, MBs);
   for (i=0; i < 4; i++)
      KillAllMMNodes(MBs[i]);
   free(od);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   fprintf(stderr,
"This will create the <pre>UTRSM subdir in <outdir>, and populate it with\n"
"all utrsm kernels and their Makefile.\n");
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr,"USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -o <path>: what directory to output all files to?\n");
   fprintf(stderr, "   -p [s,d]: set type/precision prefix (d) \n"
           "      s/d will generate for complex (c/z) as well\n");

   exit(ierr ? ierr : -1);
}

char *GetFlags(int nargs, char **args, char *PRE)
{
   FILE *fpin=stdin;
   ATL_mmnode_t *mb=NULL, *bkb=NULL, *mp;
   char *outd=NULL;
   int i;
   char pre='d';

   for (i=1; i < nargs; i++)
   {
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
         outd = DupString(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!outd)
      outd = DupString("tmp");

   *PRE = pre;
   return(outd);
}

int main(int nargs, char **args)
{
   char *outd;
   char pre;
   outd = GetFlags(nargs, args, &pre);
   GenAllTrsm(pre, outd);
   free(outd);
   return(0);
}
