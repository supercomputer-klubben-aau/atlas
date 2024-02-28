#include "atlas_mmgen.h"
#define uint unsigned int
void GetAllCopy(char pre, ATL_view_t *vb, int minSz)
{
   int d;
   ATL_cpnode_t *cb=NULL, *cbC=NULL, *cp;
   ATL_mmnode_t *ML, *mp;
   char upr;

   if (pre == 'c')
      upr = 's';
   else if (pre == 'z')
      upr = 'd';
   else
      upr = pre;
   ML = ReadMMFileWithPath(upr, "res", "FNLK1.sum");
   assert(ML);
   for (d=0; d < 2; d++)  /* 0: M2BLK, 1: BLK2M */
   {
      int a;
      for (a=0; a < 3; a++) /* 0:alp=1, 1:alp=-1, 2:alp=X */
      {
         ATL_view_t *vp;
         for (vp=vb; vp; vp = vp->next) /* scope all input views */
         {
            const int vflg=vp->flag;
            int flg;
            int b;
            ATL_mmnode_t *MLV = (vp->flag & (1<<CPV_SELFK)) ? NULL : ML;

            flg = (d) ? 0 : (1<<CPF_TOBLK);
            flg |= 1<<(CPF_AL1+a);
            if (pre == 'd')
               flg |= (1<<CPF_REAL);
            else if (pre == 's')
               flg |= (1<<CPF_REAL) | (1<<CPF_SINGLE);
            else if (pre == 'c')
               flg |= (1<<CPF_SINGLE);
/*
 *          Add copies A needs from all views (TRANS added later)
 */
            if ((vflg&(1<<(CPV_A2BLK+d))) && (vflg&(1<<(CPV_AL1A+a))))
            {
               ATL_mmnode_t *mb;
               ATL_cpnode_t *bb;
               mb = ReadMMFile(vp->fnam);
               if (!mb)
               {
                  fprintf(stderr, "NO KERNS FROM VIEW FILE '%s'!\n", vp->fnam);
                  assert(mb);
               }
               bb = GetAllCopyNodesFromMM(flg, mb, MLV);
               if (bb)
               {
                  cp = ATL_LastCPNode(bb);
                  cp->next = cb;
                  cb = bb;
                  bb = GetAllCopyNodesFromMM(flg|(1<<CPF_ABLK), mb, MLV);
                  cp = ATL_LastCPNode(bb);
                  cp->next = cb;
                  cb = bb;
               }
               KillAllMMNodes(mb);
            }

            if (!(vflg&(1<<(CPV_AL1C+a))))
               continue;
            if (!(vflg&(1<<(CPV_C2BLK+d))))
               continue;
            flg |= 1<<CPF_CBLK;
            for (b=0; b < 4; b++) /* 0:beta=1, 1:beta=N, 2:beta=X, 3:beta=0 */
            {
               flg &= ~CPF_ALLBET;
               flg |= 1<<(CPF_BE1+b);
               if (vflg & (1<<(CPV_BE1C+b)))
               {
                  ATL_mmnode_t *mb;
                  ATL_cpnode_t *bb;
                  mb = ReadMMFile(vp->fnam);
                  assert(mb);
                  bb = GetAllCopyNodesFromMM(flg, mb, MLV);
                  if (bb)
                  {
                     cp = ATL_LastCPNode(bb);
                     cp->next = cbC;
                     cbC = bb;
                  }
                  KillAllMMNodes(mb);
               }
            }
         }  /* end search thru views */
      }
   }
/*
 * Having added all copies from views that matched, get rid of dups
 */
   if (cb)
   {
      ATL_cpnode_t *b, *last;
      cb = CopyNoRep(cb, minSz);
      last = b = CloneCPQueue(cb);
      for (cp=b; cp; cp = cp->next)
      {
         cp->flag |= 1<<CPF_TRANS;
         last = cp;
      }
      last->next = cb;
      cb = b;
   }
   cbC = CopyNoRep(cbC, minSz);
   if (cb)
   {
      cp = ATL_LastCPNode(cb);
      cp->next = cbC;
   }
   else
      cb = cbC;
   for (cp=cb; cp; cp = cp->next)
   {
      if (!cp->ID && cp->rout)  /* generated filenames should be */
      {                         /* regenerated & set to ouput name */
         free(cp->rout);
         cp->rout = NULL;
      }
   }
/*
 * Check if current views same as last time
 */
   cbC = ReadCPFileWithPath(pre, "res", "cpylst.CPS");
   if (cbC)
   {
      int SAME;
      SAME = ATL_CountNumberOfCPNodes(cb) == ATL_CountNumberOfCPNodes(cbC);
      if (SAME)  /* if same # of entries, same list if */
      {          /* every kernel in cbC matches a kernel in cb */
         for (cp=cb; cp; cp = cp->next)
         {
            ATL_cpnode_t *km;
            km = FindEquivCopy(cbC, cp);
            if (km)
            {
               cbC = RemoveCPNodeFromQ(cbC, km);
               KillCPNode(km);
            }
            else
            {
               SAME = 0;
               break;
            }
         }
         if (cbC)
         {
            KillAllCPNodes(cbC);
            SAME = 0;
         }
      }
      if (!SAME)
      {
         char ln[18];
         WriteCPFileWithPath(pre, "res", "cpylst.CPS", cb);
         printf("UPDATED: res/%ccpylst.CPS!\n\n", pre);
         strcpy(ln, "res/XcpyPERF.CPS");
         ln[4] = pre;
         remove(ln);
      }
      else
         printf("UNCHANGED: res/%ccpylst.CPS!\n\n", pre);
      KillAllCPNodes(cbC);
   }
   else
   {
      WriteCPFileWithPath(pre, "res", "cpylst.CPS", cb);
      printf("CREATED NEW res/%ccpylst.CPS!\n\n", pre);
   }
   KillAllCPNodes(cb);
   KillAllMMNodes(ML);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr,
"This routine takes a list of views, and finds the minimal list of copy\n"
"routines needed to provide copy for all kernels, output in Xcpylist.CPS\n"
"which can then be read by xcpysearch for tuning\n\n");

   fprintf(stderr,"USAGE: %s [flags:\n", name);
   fprintf(stderr,"   -p [d,s,c,z]\n");
   fprintf(stderr,"   -z <minsz>: (24*24) demand mb*nb>minsz\n");
   fprintf(stderr,
   "   -V [C,A]a=1,N,X Cb=0,1,N,X [A,C]d=I,F S=A,C,M,U <name> mmview.sum:\n");
   fprintf(stderr,
"      If args to left of = are omitted, these vals don't appear in the view.\n"
"      Args to right of = may be omitted, which means all values are set.\n"
"      <name> and <mmview> must appear in that order.\n"
"      A/Cd=[I,F]: copy Into or From col-major (I,F means both directions).\n"
"      A/Ca/b give the list of needed alpha/beta for that matrix copy.\n"
"      S=A,C,M,U: don't store [A,C] copyI, matmulI, or unrollings, resp.\n"
"      Suppressed copies will still be added to master list, just not in view\n"
"      <name> is a unique string that will be used in all header files.\n"
"      mmview.sum: list of amm kerns demanding the copies.\n");
   exit(ierr ? ierr : -1);
}

ATL_view_t *GetFlags(int nargs, char **args, char *PRE, int *MINSZ)
{
   ATL_view_t *vb=NULL, *vp;
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
      case 'z':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         minSz = atoi(args[i]);
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
   *PRE = pre;
   *MINSZ = minSz;
   return(vb);
}

int main(int nargs, char **args)
{
   ATL_view_t *vb;
   int minSz;
   char pre;

   vb = GetFlags(nargs, args, &pre, &minSz);
   GetAllCopy(pre, vb, minSz);
   KillAllViews(vb);
   return(0);
}
