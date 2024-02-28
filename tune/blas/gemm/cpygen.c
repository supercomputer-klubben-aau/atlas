#include "atlas_cpparse.h"
#include "atlas_mmgen.h"
#include "atlas_iopt.h"
void GenMake(ATL_cpnode_t *cb, char *path)
{
   FILE *fp;
   ATL_cpnode_t *cp;
   int len, CPLX;
   char *fnam;
   char pre;
   char UPR;
   char dcomp[7] = {'$','(', 'D', 'K', 'C',')', '\0'};
   char dflag[12] = {'$','(','D','K','C', 'F','L','A','G','S', ')', '\0'};

   if (!cb)
      return;
   pre = CopyGetPre(cb->flag);
   CPLX = pre == 'c' || pre == 'z';
   len = strlen(path) + 10;  /* /Makefile */
   fnam = malloc(len);
   assert(fnam);
   sprintf(fnam, "%s/Makefile", path);
   assert(!Sys2File("fgrep 'BLDdir =' Make.inc", fnam));
   fp = fopen(fnam, "a");
   assert(fp);
   free(fnam);
   fprintf(fp, "include $(BLDdir)/Make.inc\n");

   UPR = (pre == 'c' || pre == 's') ? 'S' : 'D';
   dflag[2] = dcomp[2] = UPR;
   fprintf(fp, "\n   objs =");
   for (cp=cb; cp; cp = cp->next)
   {
      fprintf(fp, " \\\n          %s.o", cp->rout);
      if (CPLX && !(cp->flag&(1<<CPF_CBLK)))  /* Need conjugate target? */
      {
         ConjCopyName(cp);
         fprintf(fp, " \\\n          %s.o", cp->rout);
         UnConjCopyName(cp);
      }
   }
   PrintMakeTargs(fp, pre);
/*
 * Print compilation rules
 */
   for (cp=cb; cp; cp = cp->next)
   {
      char ch = (cp->flag&(1<<CPF_ASM)) ? 'S' : 'c';
      if (cp->comp || cp->cflags)
      {
         fprintf(fp, "%s.o : %s.%c\n", cp->rout, cp->rout, ch);
         fprintf(fp, "\t%s $(CDEFS) %s -c %s.%c", cp->comp ? cp->comp:dcomp,
                 cp->cflags ? cp->cflags : dflag, cp->rout, ch);
      }
      if (CPLX && !(cp->flag&(1<<CPF_CBLK)))  /* Need conjugate target? */
      {
         ConjCopyName(cp);
         fprintf(fp, "%s.o :", cp->rout);
         UnConjCopyName(cp);
         fprintf(fp, " %s.%c\n", cp->rout, ch);
         fprintf(fp, "\t%s $(CDEFS) %s -DConj_=1 -c %s.%c\\\n",
                 cp->comp ? cp->comp:dcomp, cp->cflags ? cp->cflags:dflag,
                 cp->rout, ch);
         ConjCopyName(cp);
         fprintf(fp, "           -o %s.o\n", cp->rout);
         UnConjCopyName(cp);
      }
   }
/*
 * Print default rules
 */
   fprintf(fp, ".S.o:\n");
   fprintf(fp, "\t%s $(CDEFS) %s -c $*.S\n", dcomp, dflag);
   fprintf(fp, ".s.o:\n");
   fprintf(fp, "\t%s $(CDEFS) %s -c $*.s\n", dcomp, dflag);
   fprintf(fp, ".c.o:\n");
   fprintf(fp, "\t%s $(CDEFS) %s -c $*.c\n", dcomp, dflag);

   fclose(fp);
}

void GenCopy(ATL_cpnode_t *cb, char *path)
{
   ATL_cpnode_t *cp;
   char *sp, *typs, *untyp;
   int len, plen;
   char pre;

   if (!cb)
      return;

   plen = strlen(path);
   pre = CopyGetPre(cb->flag);
   untyp = "#ifdef SCALAR\n   #undef SCALAR\n#endif\n"
           "#ifdef TYPE\n   #undef TYPE\n#endif\n"
           "#ifdef SREAL\n   #undef SREAL\n#endif\n"
           "#ifdef DREAL\n   #undef DREAL\n#endif\n"
           "#ifdef SCPLX\n   #undef SCPLX\n#endif\n"
           "#ifdef DCPLX\n   #undef DCPLX\n#endif\n"
           "#ifdef ALPHA1\n  #undef ALPHA1\n#endif\n"
           "#ifdef ALPHAN\n  #undef ALPHAN\n#endif\n"
           "#ifdef ALPHAN1\n  #undef ALPHAN1\n#endif\n"
           "#ifdef ALPHAX\n  #undef ALPHAX\n#endif\n"
           "#ifdef BETA1\n  #undef BETA1\n#endif\n"
           "#ifdef BETAN\n  #undef BETAN\n#endif\n"
           "#ifdef BETAN1\n  #undef BETAN1\n#endif\n"
           "#ifdef BETAX\n  #undef BETAX\n#endif\n";
   if (pre == 'z')
      typs = "#define DCPLX 1\n#define TYPE double\n#define SCALAR TYPE*\n";
   else if (pre == 'c')
      typs = "#define SCPLX 1\n#define TYPE float\n#define SCALAR TYPE*\n";
   else if (pre == 'd')
      typs = "#define DREAL 1\n#define TYPE double\n#define SCALAR TYPE\n";
   else if (pre == 's')
      typs = "#define SREAL 1\n#define TYPE float\n#define SCALAR TYPE\n";

   for (cp=cb; cp; cp = cp->next)
   {
      char *nam, *fnam=NULL, *rt=cp->rout;
      char calp, cbet=0, fe = 'c';
      FILE *fp;

      calp = CopyGetAlphaC(cp->flag);
      if (cp->flag&(1<<CPF_CBLK))
         cbet = CopyGetBetaC(cp->flag);
      nam = GetCopyName(cp, 0);
      if (cp->ID)
      {
         char *onam;
         assert(rt);
         for (onam=rt; *onam; onam++);
         fe = onam[-1];
         if (fe != 'c' && fe != 'S' && fe != 's')
            fe = 'c';

         len = 10 + strlen(rt);
         cp->rout = malloc(len);
         assert(cp->rout);
         if (strncmp(rt, "CPYCASES/", 9))
            sprintf(cp->rout, "CPYCASES/%s", rt);
         else
            strcpy(cp->rout, rt);
      }
      else
      {
         char *gens=cp->genstr;
         cp->rout = "ATL_tmp.c";
         cp->genstr = GetCopyGenStr(cp);
         MMDoGenString(0, cp->genstr);
         free(cp->genstr);
         cp->genstr = gens;
      }
      if (rt)
         free(rt);
      if (fe == 's' || fe == 'S')
         cp->flag |= 1L<<CPF_ASM;
      len = strlen(path) + 1 + strlen(nam) + 2 + 1;
      fnam = malloc(len);
      sprintf(fnam, "%s/%s.%c", path, nam, fe);
      fp = fopen(fnam, "w");
      assert(fp);
      fprintf(fp, untyp);
      fprintf(fp, typs);
      if (cp->flag&(1<<CPF_CBLK))
      {
         fprintf(fp, "#define BETA%c\n", cbet);
         if (cbet == 'N')
            fprintf(fp, "#define BETAN1 1\n");
         fprintf(fp, "#define ATL_MU %u\n#define ATL_NU %u\n", cp->mu, cp->nu);
      }
      else
         fprintf(fp, "#define ATL_NU %u\n", cp->nu);
      fprintf(fp, "#define ALPHA%c 1\n", calp);
      if (calp == 'N')
         fprintf(fp, "#define ALPHAN1 1\n");
      fprintf(fp, "#define ATL_USERCPMM %s\n", nam);
      if (!(cp->flag&(1<<CPF_REAL)) && !(cp->flag&(1<<CPF_CBLK)))
      {  /* if complex A/B copy, need conj_! */
         char *rt = cp->rout;
         cp->rout = nam;
         fprintf(fp, "#ifdef Conj_\n   #undef ATL_USERCPMM\n");
         ConjCopyName(cp);
         fprintf(fp, "   #define ATL_USERCPMM %s\n#endif\n", nam);
         UnConjCopyName(cp);
         cp->rout = rt;
      }
      fprintf(fp, "\n");
      fclose(fp);
      len = 10 + strlen(cp->rout) + strlen(fnam);
      sp = malloc(len);
      assert(sp);
      sprintf(sp, "cat %s >> %s", cp->rout, fnam);
      free(fnam);
      if (system(sp))
      {
         fprintf(stderr, "ERROR IN CMD: '%s'\n", sp);
         assert(0);
      }
      free(sp);
      if (cp->ID)
         free(cp->rout);
      cp->rout = nam;
   }
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -i cpyout.sum: repeat for multiple\n");
   fprintf(stderr, "   -o <output path>\n");
   exit(ierr ? ierr : -1);
}

ATL_cpnode_t *GetFlags(int nargs, char **args, char **OUTD)
{
   ATL_cpnode_t *cb=NULL, *cp;
   char *outd=NULL;
   int i;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'i':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        if (!cb)
           cb = ReadCPFile(args[i]);
        else
        {
           cp = ATL_LastCPNode(cb);
           cp->next = ReadCPFile(args[i]);
        }
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
   if (!cb)
      cb = ReadCPFile(NULL);
   if (!outd)
      outd = DupString("tmp/");
   *OUTD = outd;
   return(cb);
}

int main(int nargs, char **args)
{
   ATL_cpnode_t *cb;
   char *outd;

   cb = GetFlags(nargs, args, &outd);
   assert(cb);
   GenCopy(cb, outd);   /* sets info for Protos&Make! */
   GenMake(cb, outd);
   free(outd);
   KillAllCPNodes(cb);
   return(0);
}
