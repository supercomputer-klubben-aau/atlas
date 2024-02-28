#include "atlas_mmgen.h"
#include "atlas_type.h"
void PrintProtos(FILE *fp, ATL_cpnode_t *cb)
{
   ATL_cpnode_t *cp;
   int i;
   char *args;

   if (!cb)
      return;
   cb = CloneCPQueue(cb);
   cb = CopyNoRep(cb, 24*24);
   if (cb->flag & (1<<CPF_CBLK))
   {
      if (cb->flag & (1<<CPF_TOBLK)) /* copy from C into blk */
      {
         if (cb->flag & (1<<CPF_REAL))
            args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,"
                   "const SCALAR, TYPE *";
         else
            args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,"
                   "const SCALAR, TYPE*,TYPE*";
      }
      else /* copy from block into C */
      {
         if (cb->flag & (1<<CPF_REAL))
            args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,"
                   "const SCALAR,TYPE*,ATL_CSZT";
         else
            args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,"
                   "const SCALAR,TYPE*,ATL_CSZT";
      }
   }
   else if (cb->flag & (1<<CPF_TOBLK)) /* A2blk */
   {
      if (cb->flag & (1<<CPF_REAL))
         args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*";
      else
         args = "ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,"
                "TYPE*,TYPE*";
   }
   else /* blk2A */
   {
      if (cb->flag & (1<<CPF_REAL))
         args = "ATL_CSZT,ATL_CSZT,const SCALAR, const TYPE*,ATL_CSZT,TYPE*";
      else
         args = "ATL_CSZT,ATL_CSZT,const SCALAR, const TYPE*,ATL_CSZT,"
                "TYPE*,TYPE*";
   }
   for (cp=cb; cp; cp = cp->next)
      fprintf(fp, "void %s\n   (%s);\n", cp->rout, args);
   KillAllCPNodes(cb);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   fprintf(stderr,
      "This program generates copy header files for lists of kernels.\n"
      "It will generate seperate header file for each alpha/beta combo.\n");
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -p [d,s,c,z]\n");
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

ATL_view_t *GetFlags(int nargs, char **args, char *PRE, char **PTH)
{
   ATL_view_t *vb=NULL;
   char *pth=NULL, *nm=NULL;
   int i;
   char pre='d', dir='I';

   for (i=1; i < nargs; i++)
   {
      ATL_view_t *vp;
      char *nam;
      int k, flag=0;

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
   *PRE = pre;
   *PTH = pth;
   return(vb);
}

void GenHead
(
   char pre,        /* s,d,c,z */
   char mtx,        /* C,A,B */
   char dir,        /* 'F': FromC, i.e. C2blk, 'I': blk2C */
   char calp,       /* 1,N,X */
   char cbet,       /* C:0,1,N,X or A:T,N */
   char *nm,        /* view name */
   char *pnam,      /* name with path */
   char *nam,       /* name without path */
   ATL_cpnode_t *cb /* queue of copies to produce header for */
)
{
   ATL_cpnode_t *cp;
   FILE *fpp, *fpd;
   char *sp, *arrnm;
   int i, flag, N;

   assert(cb);
   sp = strstr(nam, "ZZZZ");
   assert(sp);
   sp[0] = 'p';
   sp[1] = 'e';
   sp[2] = 'r';
   sp[3] = 'f';
   fpp = fopen(pnam, "w");
   assert(fpp);
   PrintBegIfdef(fpp, nam);
   sp[0] = 'd';
   sp[1] = 'e';
   sp[2] = 'c';
   sp[3] = 'l';
   fpd = fopen(pnam, "w");
   assert(fpd);
   PrintBegIfdef(fpd, nam);
   sp[0] = 'Z';
   sp[1] = 'Z';
   sp[2] = 'Z';
   sp[3] = 'Z';

   fprintf(fpd, "\n#include \"atlas_amm.h\"\n");
   fprintf(fpd, "/*\n * Prototypes\n */\n");
   PrintProtos(fpd, cb);
   fprintf(fpd, "/*\n * Index array\n */\n");
   N = ATL_CountNumberOfCPNodes(cb);
   if (mtx == 'C')
   {
      fprintf(fpd, "static const %s_t ATL_AMM_%s_a%c_b%c[%d] =\n{\n",
              (dir=='F')?"cmat2ablk":"ablk2cmat", (dir=='F')?"BLK2C":"C2BLK",
              calp, cbet, N);
      fprintf(fpp, "static const float ATL_AMM_%sTIME_a%cb%c[%d] = \n{\n",
              (dir=='F')?"C2BLK":"BLK2C", calp, cbet, N);
   }
   else /* A or B */
   {
      char *bn;
      if (dir == 'F')
      {
         if (cbet == 'N')
            bn = (mtx == 'A') ? "AN2BLK" : "BN2BLK";
         else if (cbet == 'C')
            bn = (mtx == 'A') ? "AC2BLK" : "BC2BLK";
         else if (cbet == 'H')
            bn = (mtx == 'A') ? "AH2BLK" : "BH2BLK";
         else
            bn = (mtx == 'A') ? "AT2BLK" : "BT2BLK";
         fprintf(fpd,"static const cm2am_t ATL_AMM_%s_a%c[%d] =\n{\n",
                 bn, calp, N);
      }
      else
      {
         if (cbet == 'N')
            bn = (mtx == 'A') ? "BLK2AN" : "BLK2BN";
         else if (cbet == 'C')
            bn = (mtx == 'A') ? "BLK2AC" : "BLK2BC";
         else if (cbet == 'H')
            bn = (mtx == 'A') ? "BLK2AH" : "BLK2BH";
         else
            bn = (mtx == 'A') ? "BLK2AT" : "BLK2BT";
         if (mtx == 'A')
            bn = (cbet == 'N') ? "BLK2AN" : "BLK2AT";
         else
            bn = (cbet == 'N') ? "BLK2BN" : "BLK2BT";
         fprintf(fpd,"static const am2cm_t ATL_AMM_%s_a%c[%d] =\n{\n",
                 bn, calp, N);
      }
      fprintf(fpp, "static const float ATL_AMM_%stime_a%c[%d] = \n{\n",
              bn, calp, N);
   }
   for (i=0,cp=cb; cp; i++, cp = cp->next)
   {
      fprintf(fpd, "/* IDX=%3d */ %s%c\n", i, cp->rout, cp->next?',':' ');
      fprintf(fpp, "/* IDX=%3d */ %e%c\n", i, cp->mflop[0], cp->next?',':' ');
   }
   fprintf(fpd, "};\n\n");
   fprintf(fpp, "};\n\n");

   fprintf(fpp, "#endif  /* end ifdef guard */\n");
   fprintf(fpd, "#endif  /* end ifdef guard */\n");
   fclose(fpp);
   fclose(fpd);
}


int *ReadMasterIdx(char pre, int *N)
{
   FILE *fp;
   int *flags;
   char fn[16]={'r','e','s','/',pre,'m','a','s','t','e','r','.','C','P','I',0};
   unsigned int i, n;

   fp = fopen(fn, "r");
   assert(fp);
   assert(fscanf(fp, "%u", &n) == 1);
   flags = malloc(n*sizeof(int));
   assert(flags);
   for (i=0; i < n; i++)
      assert(fscanf(fp, " %x '%*s' ", flags+i) == 1);

   *N = n;
   return(flags);
}

void DoMaster(char pre, FILE *fp, ATL_cpnode_t *cb)
{
   ATL_cpnode_t *cp;
   char *sp;
   int N, i;

   assert(cb);
   fprintf(fp, "#include \"atlas_amm.h\"\n/*\n * Prototypes\n */\n");
   PrintProtos(fp, cb);

   sp = CopyFlag2Str(cb->flag);
   fprintf(fp, "/*\n * Index array\n */\n");
   N = ATL_CountNumberOfCPNodes(cb);
   fprintf(fp, "#define ATL_Cpy%s ATL_%c%s\n", sp, pre, sp);
   fprintf(fp, "#ifndef ATL_DECL_\n");
   fprintf(fp, "   #ifndef ATL_MACRO_ONLY_\n");
   fprintf(fp, "      const extern ATL_cparr_t ATL_Cpy%s[%u];\n", sp, N+N);
   fprintf(fp, "   #endif\n");
   fprintf(fp, "#else\nconst ATL_cparr_t ATL_Cpy%s[%u] =\n{", sp, N+N);
   sp = "\n";
   for (i=0,cp=cb; cp; cp = cp->next, i++)
   {
      ATL_cpflt_t mf, *mfp;
      ATL_iptr_t imf=0;

      mf = 1.0 / (cp->mflop[0]*cp->mb*cp->nb);
      mfp = (ATL_cpflt_t*) &imf;
      *mfp = mf;
      fprintf(fp, "%s   (ATL_cparr_t)%s, ", sp, cp->rout);
      if (sizeof(imf) > sizeof(int))
         fprintf(fp, "0x%lxL", imf);
      else
         fprintf(fp, "0x%x", (int)imf);
      fprintf(fp, "/*%u:%.3e*/", i, (float)mf);
      sp = ",\n";
   }
   fprintf(fp, "\n};\n");
   fprintf(fp, "#endif /* end else of ifndef ATL_DECL */\n");
   fprintf(fp, "\n");
}

void GenMasterFiles(char pre, char *path)
{
   int *flags;
   char *fn;
   int i, pL, L, N;

   pL = strlen(path);
   L = pL + 32;

   fn = malloc(pL+32);
   assert(fn);
   flags = ReadMasterIdx(pre, &N);
   for (i=0; i < N; i++)
   {
      char *sp;
      ATL_cpnode_t *cb;
      FILE *fp;
      const int NC = (pre == 'z' || pre == 'c') ? 2 : 1;
      int k, ic;

      sp = CopyFlag2Str(flags[i]);
      k = sprintf(fn, "res/%ccpyPERF_%s.CPS", pre, sp);
      assert(k < L);
      cb = ReadCPFile(fn);
      assert(cb);
      for (ic=0; ic < NC; ic++)
      {
         if (ic)
         {
            ATL_cpnode_t *cp;
            sp = CopyFlag2Str(flags[i]|(1<<CPF_CONJ));
            for (cp=cb; cp; cp = cp->next)
            {
               cp->flag |= 1<<CPF_CONJ;
               ConjCopyName(cp);
            }
         }
         k = sprintf(fn, "%s/atlas_%c%s.h", path, pre, sp);
         assert(k < L);
         fp = fopen(fn, "w");
         assert(fp);
         sp = fn + pL;
         while (*sp == '/')
            sp++;
         while (*sp != '.')
         {
            *sp = toupper(*sp);
            sp++;
         }
         *sp++ = '_';
         *sp++ = 'H';
         *sp = '\0';
         sp = fn + pL + 1;
         fprintf(fp, "#ifndef %s\n", sp);
         fprintf(fp, "   #define %s 1\n\n", sp);
         DoMaster(pre, fp, cb);
         fprintf(fp, "#endif /* end multiple inclusion guard */\n");
         fclose(fp);
      }
      KillAllCPNodes(cb);
   }
   free(fn);
   free(flags);
}

unsigned int MMIdxCount(ATL_mmnode_t *mb, ATL_mmnode_t *want)
/*
 * RETURNS: index (from 1) of want in mb, 0 if not found
 */
{
   ATL_mmnode_t *mp;
   int i;
   for (i=0, mp=mb; mp; mp = mp->next, i++)
      if (MMKernCompsSame(mp, want))
         return(i+1);
   return(0);
}

double FlopCount(ATL_mmnode_t *mp)
{
   double cnt;
   switch(mp->blask)
   {
   case 1:
      cnt = ((1.0*mp->mbB)*mp->nbB)*mp->kbB;
      break;
   default:
      cnt = ((2.0*mp->mbB)*mp->nbB)*mp->kbB;
   }
   return(cnt);
}
void GenViewFile(char pre, char *path, ATL_view_t *vp, ATL_mmnode_t *mb,
                 ATL_mmnode_t *ML)
{
   ATL_mmnode_t *mp, *MB, *maxPerf=NULL;
   FILE *fpout;
   #define NTHRSH 11
   int THRSH[NTHRSH] = {25, 33, 50, 66, 75, 80, 85, 90, 95, 98, 99};
   ATL_mmnode_t *mpT[NTHRSH];
   char *sp, *fn, *VN;
   char *mnms[10] = {"m", "n", "k", "amm", "A2b", "B2b", "b2C",
                    "b2A", "b2B", "C2b"};
   char *nms[10]={"MB", "NB", "KB", "AMMI",
                "A2BLK", "B2BLK", "BLK2C", "BLK2A", "BLK2B", "C2BLK"};
   int nbits[10], idxT[NTHRSH];
   char **cpnms = nms + 4;
   int *cpidxs[6], *cpnbits=nbits+4;
   double mfMax=0.0;
   int i, d, plen, vlen, len, N, nint, nbit, nkvec, idxMax=0;
   int minMB, minNB, minKB, maxMB, maxNB, maxKB, minKID, maxKID, tbits;
   const unsigned int vflag=vp->flag;
   char upr = pre;
/*
 * Need these because the frontend machine could conceivably have a smaller
 * int type than the machine being compiled for.  Would need to have the
 * tuned machine do the generation to avoid this assertion.  I also assume
 * both machines have same float, but ATLAS typically assumes IEEE storage
 * so this is less risky
 */
   #if (ATL_ISIZE < ATL_FSIZE)
      #error "Integer shorter than float!"
   #endif
   assert(ATL_ISIZE <= sizeof(int));
   assert(mb);
   assert(vp);
   if (pre == 'c')
      upr = 's';
   else if (pre == 'z')
      upr = 'd';
   plen = strlen(path);
   vlen = strlen(vp->nam);
   VN = vp->nam;
   len = plen + 8 + 1 + vlen + 8;
   len = (len < 33) ? 33 : len;  /* need room for cpyPERF_<str>.CPS */
   fn = malloc(len*sizeof(char));
   assert(fn);
   d = sprintf(fn, "%s/atlas_%c%s_view.h", path, pre, vp->nam);
   assert(d < len);
   fpout = fopen(fn, "w");
   assert(fpout);
   for (mp=mb, N=nkvec=0; mp; mp = mp->next, N++)
      if (FLAG_IS_SET(mb->flag, MMF_KVEC))
         nkvec++;
   sp = fn + plen + 1;
   for (i=0; sp[i]; i++)
      sp[i] = toupper(sp[i]);
   sp[i-2] = '_';
   sp[i-1] = 'H';
   fprintf(fpout, "#ifndef %s\n   #define %s 1\n\n", sp, sp);

   for (d=0; d < 2; d++)
   {
      unsigned int m;
      for (m=0; m < 3; m++)
      {
         int FIND=0, i, flag;
         i = d+(d<<1) + m; /* i = 3*d+m */
         flag = (pre == 'd' || pre == 's') ? (1<<CPF_REAL) : 0;
         flag |= (pre == 's' || pre == 'c') ? (1<<CPF_SINGLE) : 0;
         if (m == 2) /* C idx */
         {
            int b;
            flag |= 1<<CPF_CBLK;
            flag |= (d) ? (1<<CPF_TOBLK) : 0;
            for (b=0; b < 4; b++)
            {
               if (vflag & 1<<(CPV_BE1C+b))
               {
                  flag |= 1<<CPF_BE1+b;
                  break;
               }
            }
            for (b=0; b < 3; b++)
            {
               if (vflag & 1<<(CPV_AL1C+b))
               {
                  flag |= 1<<CPF_AL1+b;
                  break;
               }
            }
            if (vflag & (1<<CPV_BE1C))
               flag |= 1<<CPF_BE1;
            else if (vflag & (1<<CPV_BENC))
               flag |= 1<<CPF_BEN;
            else if (vflag & (1<<CPV_BEXC))
               flag |= 1<<CPF_BEX;
            if ((vflag&(CPV_ALLBEC)) && (vflag&(CPV_ALLALC)))
            {
               if (!d && (vflag&(1<<CPV_BLK2C)))
                  FIND = 1;
               else if (d && (vflag&(1<<CPV_C2BLK)))
                  FIND = 1;
            }
         }
         else /* A or B */
         {
            int b;

            flag |= (m == 0) ? (1<<CPF_ABLK) : 0;
            flag |= (d) ? 0 : (1<<CPF_TOBLK);
            for (b=0; b < 3; b++)
            {
               if (vflag & 1<<(CPV_AL1A+b))
               {
                  flag |= 1<<(CPF_AL1+b);
                  break;
               }
            }
            if (vflag & CPV_ALLALA)
            {
               if (!d && (vflag&(1<<CPV_A2BLK)))
                  FIND = 1;
               else if (d && (vflag&(1<<CPV_BLK2A)))
                  FIND = 1;
            }
         }
/*
 *       If view uses this mtx,dir setting, allocate an N+1 iarray, with
 *       1st entry being the minimum value found, rest being offset from min.
 */
         if (FIND)
         {
            int min, k;
            ATL_cpnode_t *cb;
            /*cpflags[i] = flag; */
            sp = CopyFlag2Str(flag);
            k = sprintf(fn, "res/%ccpyPERF_%s.CPS", pre, sp);
            assert(k < len);
            cb = ReadCPFile(fn);
            if (!cb)
            {
               fprintf(stderr, "Missing or empty file '%s'\n", fn);
               assert(cb);
            }
            cpidxs[i] = MMGetCopyIdxsFromList(flag, mb, cb, ML);
            KillAllCPNodes(cb);
            assert(cpidxs[i]);
            assert(cpidxs[i][0] == N);
            min = cpidxs[i][1];
            for (k=1; k < N; k++)
               min = Mmin(min, cpidxs[i][k+1]);
            for (k=1; k <= N; k++)
               cpidxs[i][k] -= min;
            cpidxs[i][0] = min;
            min = cpidxs[i][1];  /* now actually maximum */
            for (k=2; k <= N; k++)
               min = Mmax(min, cpidxs[i][k]);
            cpnbits[i] = ATL_numBitsNeeded(min);
         }
         else
         {
            cpidxs[i] = NULL;
            cpnbits[i] = 0;
         }
      }
   }
   for (tbits=i=0; i < 6; i++)
   {
      fprintf(fpout, "#define ATL_VW%s_MIN_%s      %3u\n", VN, cpnms[i],
              cpidxs[i] ? cpidxs[i][0] : 0);
      fprintf(fpout, "#define ATL_VW%s_NBIT_%s     %3u\n", VN,
              cpnms[i], cpnbits[i]);
      tbits += cpnbits[i];
   }
   maxMB = minMB = mb->mbB;
   maxNB = minNB = mb->nbB;
   maxKB = minKB = mb->kbB;
   if (vflag&(1<<CPV_NOMMI))
   {
      MB = NULL;
      maxKID = minKID = 0;
   }
   else
   {
      MB = ReadMMFileWithPath(upr, "res", "FNLK1.sum");
      assert(MB);
      maxKID = MMIdxCount(MB, mb);
      if (!maxKID)
      {
         fprintf(stderr, "Cannot find ammI for:\n");
         PrintMMLine(stderr, mb);
      }
      assert(maxKID);
      maxKID--;
      minKID = maxKID;
   }
   for (i=0; i < NTHRSH; i++)
      mpT[i] = NULL;
   mfMax = mb->mflop[0];
   idxMax = 0;
   maxPerf = mb;
   for (i=1, mp=mb->next; mp; mp = mp->next, i++)
   {
      int b;
      ATL_mmnode_t p;

      if (mp->mflop[0] > mfMax)
      {
         mfMax = mp->mflop[0];
         idxMax = i;
         maxPerf = mp;
      }
      b = mp->mbB;
      maxMB = Mmax(b, maxMB);
      minMB = Mmin(b, minMB);
      b = mp->nbB;
      maxNB = Mmax(b, maxNB);
      minNB = Mmin(b, minNB);
      b = mp->kbB;
      maxKB = Mmax(b, maxKB);
      minKB = Mmin(b, minKB);
      if (MB)
      {
         b = MMIdxCount(MB, mp);
         assert(b);
         b--;
         maxKID = Mmax(b, maxKID);
         minKID = Mmin(b, minKID);
      }
   }
   nbits[0] = ATL_numBitsNeeded(maxMB-minMB);  /* now set vals to */
   nbits[1] = ATL_numBitsNeeded(maxNB-minNB);  /* number of bits needed */
   nbits[2] = ATL_numBitsNeeded(maxKB-minKB);  /* to store offsets from min */
   nbits[3] = ATL_numBitsNeeded(maxKID-minKID);
   tbits += nbits[0] + nbits[1] + nbits[2] + nbits[3];
   fprintf(fpout, "#define ATL_VW%s_MIN_AMMI  %8u\n", VN, minKID);
   fprintf(fpout, "#define ATL_VW%s_NBIT_AMMI %8u\n", VN, nbits[3]);
   fprintf(fpout, "#define ATL_VW%s_MIN_MB    %8u\n", VN, minMB);
   fprintf(fpout, "#define ATL_VW%s_NBIT_MB   %8u\n", VN, nbits[0]);
   fprintf(fpout, "#define ATL_VW%s_MIN_NB    %8u\n", VN, minNB);
   fprintf(fpout, "#define ATL_VW%s_NBIT_NB   %8u\n", VN, nbits[1]);
   fprintf(fpout, "#define ATL_VW%s_MIN_KB    %8u\n", VN, minKB);
   fprintf(fpout, "#define ATL_VW%s_NBIT_KB   %8u\n", VN, nbits[2]);
   fprintf(fpout, "#define ATL_VW%s_NKVEC     %8u\n", VN, nkvec);
   fprintf(fpout, "\n");
   mp = maxPerf;
   i = idxMax;
   fprintf(fpout, "#define ATL_VW%s_BEST_IDX  %8u\n", VN, idxMax);
   fprintf(fpout, "#define ATL_VW%s_BEST_MB   %8u\n", VN, mp->mbB);
   fprintf(fpout, "#define ATL_VW%s_BEST_NB   %8u\n", VN, mp->nbB);
   fprintf(fpout, "#define ATL_VW%s_BEST_KB   %8u\n", VN, mp->kbB);
   fprintf(fpout, "#define ATL_VW%s_BEST_MU   %8u\n", VN, mp->mu);
   fprintf(fpout, "#define ATL_VW%s_BEST_NU   %8u\n", VN, mp->nu);
   fprintf(fpout, "#define ATL_VW%s_BEST_KU   %8u\n", VN, mp->ku);
   fprintf(fpout, "#define ATL_VW%s_BEST_LCMU %8u\n", VN,
            ATL_iLCM(ATL_iLCM(mp->mu,mp->nu),mp->ku));
   fprintf(fpout, "#define ATL_VW%s_BEST_LCMN %8u\n", VN,
           ATL_iLCM(mp->mu,mp->nu));
   fprintf(fpout, "#define ATL_VW%s_BEST_AMMI %8u\n", VN,
           MMIdxCount(MB?MB:mb, mp));
   fprintf(fpout, "#define ATL_VW%s_BEST_A2BKI %7u\n", VN,
           cpidxs[0]?cpidxs[0][i+1]+cpidxs[0][0]:0);
   fprintf(fpout, "#define ATL_VW%s_BEST_B2BKI %7u\n", VN,
           cpidxs[1]?cpidxs[1][i+1]+cpidxs[1][0]:0);
   fprintf(fpout, "#define ATL_VW%s_BEST_BK2CI %7u\n", VN,
           cpidxs[2]?cpidxs[2][i+1]+cpidxs[2][0]:0);
   fprintf(fpout, "#define ATL_VW%s_BEST_BK2AI %7u\n", VN,
           cpidxs[3]?cpidxs[3][i+1]+cpidxs[3][0]:0);
   fprintf(fpout, "#define ATL_VW%s_BEST_BK2BI %7u\n", VN,
           cpidxs[4]?cpidxs[4][i+1]+cpidxs[4][0]:0);
   fprintf(fpout, "#define ATL_VW%s_BEST_C2BKI %7u\n", VN,
           cpidxs[5]?cpidxs[5][i+1]+cpidxs[5][0]:0);
   fprintf(fpout, "\n");
   for (i=0,mp=mb; mp->next; i++, mp = mp->next);
   fprintf(fpout, "#define ATL_VW%s_LAST_IDX  %8u\n", VN, i);
   fprintf(fpout, "#define ATL_VW%s_LAST_MB   %8u\n", VN, mp->mbB);
   fprintf(fpout, "#define ATL_VW%s_LAST_NB   %8u\n", VN, mp->nbB);
   fprintf(fpout, "#define ATL_VW%s_LAST_KB   %8u\n", VN, mp->kbB);
   fprintf(fpout, "#define ATL_VW%s_LAST_MU   %8u\n", VN, mp->mu);
   fprintf(fpout, "#define ATL_VW%s_LAST_NU   %8u\n", VN, mp->nu);
   fprintf(fpout, "#define ATL_VW%s_LAST_KU   %8u\n", VN, mp->ku);
   fprintf(fpout, "#define ATL_VW%s_LAST_LCMU %8u\n", VN,
            ATL_iLCM(ATL_iLCM(mp->mu,mp->nu),mp->ku));
   fprintf(fpout, "#define ATL_VW%s_LAST_LCMN %8u\n", VN,
           ATL_iLCM(mp->mu,mp->nu));
   fprintf(fpout, "#define ATL_VW%s_LAST_AMMI %8u\n", VN,
           MMIdxCount(MB?MB:mb, mp));
   fprintf(fpout, "#define ATL_VW%s_LAST_A2BKI %7u\n", VN,
           cpidxs[0]?cpidxs[0][i+1]+cpidxs[0][0]:0);
   fprintf(fpout, "#define ATL_VW%s_LAST_B2BKI %7u\n", VN,
           cpidxs[1]?cpidxs[1][i+1]+cpidxs[1][0]:0);
   fprintf(fpout, "#define ATL_VW%s_LAST_BK2CI %7u\n", VN,
           cpidxs[2]?cpidxs[2][i+1]+cpidxs[2][0]:0);
   fprintf(fpout, "#define ATL_VW%s_LAST_BK2AI %7u\n", VN,
           cpidxs[3]?cpidxs[3][i+1]+cpidxs[3][0]:0);
   fprintf(fpout, "#define ATL_VW%s_LAST_BK2BI %7u\n", VN,
           cpidxs[4]?cpidxs[4][i+1]+cpidxs[4][0]:0);
   fprintf(fpout, "#define ATL_VW%s_LAST_C2BKI %7u\n", VN,
           cpidxs[5]?cpidxs[5][i+1]+cpidxs[5][0]:0);
   fprintf(fpout, "\n");

   if (tbits > (ATL_ISIZE<<3))
   {
      for (nbit=nint=i=0; i < 9; i++)
      {
         nbit += nbits[i];
         if (nbit > (ATL_ISIZE<<3))
         {
            nbit = nbits[i];
            nint++;
         }
      }
      if (nbit)
         nint++;
   }
   else
      nint = 1;
   fprintf(fpout, "#define ATL_VW%s_NINT     %u\n", VN, ++nint);
   if (nint == 2)
      fprintf(fpout, "#define ATL_VW%s_IDXMUL(i_) ((i_)+(i_))\n", VN);
   else if (nint == 3)
      fprintf(fpout, "#define ATL_VW%s_IDXMUL(i_) ((i_)+((i_)<<1))\n", VN);
   else if (nint == 4)
      fprintf(fpout, "#define ATL_VW%s_IDXMUL(i_) ((i_)<<2)\n", VN);
   else
      fprintf(fpout, "#define ATL_VW%s_IDXMUL(i_) ((i_)*%u)\n", VN, nint);
   fprintf(fpout, "#define ATL_VW%s_NCASES   %u\n", VN, N);
   fprintf(fpout, "#define ATL_VW%s_MAX_MB   %u\n", VN, maxMB);
   fprintf(fpout, "#define ATL_VW%s_MAX_NB   %u\n", VN, maxNB);
   fprintf(fpout, "#define ATL_VW%s_MAX_KB   %u\n", VN, maxKB);
   fprintf(fpout, "#define ATL_VIEW_%s ATL_%cVIEW_%s\n", VN, pre, VN);

   fprintf(fpout,"\n#ifndef ATL_DECL_\n   const extern unsigned int "
           "ATL_%cVIEW_%s[%u];\n", pre, VN, nint*N);
   fprintf(fpout, "#else\nconst unsigned int ATL_%cVIEW_%s[%u] =\n{",
           pre, VN, nint*N);
   sp = "\n";
   for (i=0, mp=mb; mp; mp = mp->next, i++)
   {
      int k, iv=0, ammI;
      float *timepflop;
/*
 *    Setup threshhold blocking factors, don't go past max perf kernel
 */
      for (k=0; k < NTHRSH; k++)
      {
         if (!mpT[k] && THRSH[k]*0.01*mfMax < mp->mflop[0])
         {
            if (i < idxMax)
            {
               mpT[k] = mp;
               idxT[k] = i;
            }
            else
            {
               mpT[k] = maxPerf;
               idxT[k] = idxMax;
            }
         }
      }
      timepflop = (float*) &k;
      *timepflop = 1.0 / (mp->mflop[0]*FlopCount(mp));
      if (MB)
      {
         ammI = MMIdxCount(MB, mp);
         assert(ammI);
         ammI--;
      }
      else
         ammI = 0;
      fprintf(fpout, "%s/* %u:SPF=%e B=(%u,%u,%u) ammI=%u",
              sp, i, *timepflop, mp->mbB, mp->nbB, mp->kbB, ammI);
      if (cpidxs[0] || cpidxs[1] || cpidxs[2])
         fprintf(fpout, " [A,B]2b=%u,%u b2C=%u",
                 cpidxs[0]?cpidxs[0][i+1]+cpidxs[0][0]:0,
                 cpidxs[1]?cpidxs[1][i+1]+cpidxs[1][0]:0,
                 cpidxs[2]?cpidxs[2][i+1]+cpidxs[2][0]:0);
      if (cpidxs[3] || cpidxs[4] || cpidxs[5])
         fprintf(fpout, " b2[A,B]=%u,%u C2b=%u",
                 cpidxs[3]?cpidxs[3][i+1]+cpidxs[3][0]:0,
                 cpidxs[4]?cpidxs[4][i+1]+cpidxs[4][0]:0,
                 cpidxs[5]?cpidxs[5][i+1]+cpidxs[5][0]:0);
     fprintf(fpout, " */\n");
      fprintf(fpout, "0x%x", k);
/*<ammtime><mb><nb><kb><kernI><a2blkI><b2blkI><blk2cI><blk2aI><blk2bI><c2blkI>*/
      for (nbit=k=0; k < 10; k++)
      {
         const unsigned int nvbits = nbits[k];
         if (nvbits)
         {
            int v;
            switch(k)
            {
            case 0:
               v = mp->mbB - minMB;
               break;
            case 1:
               v = mp->nbB - minNB;
               break;
            case 2:
               v = mp->kbB - minKB;
               break;
            case 3:
               v = ammI - minKID;
               break;
            default:
               v = cpidxs[k-4] ? cpidxs[k-4][i+1] : 0;
            }
            if (nbit+nvbits > (ATL_ISIZE<<3))
            {
               fprintf(fpout, ",0x%x", iv);
               iv = 0;
               nbit = 0;
            }
            iv |= v << nbit;
            nbit += nvbits;
         }
         sp = ",\n";
      }
      if (nbit != 0)
         fprintf(fpout, ",0x%x", iv);
   }
   KillAllMMNodes(MB);
   for (d=0; d < 6; d++)
      free(cpidxs[d]);
   fprintf(fpout, "\n};\n");
   fprintf(fpout, "#endif /* end else of ifndef ATL_DECL_ */\n\n");

   for (i=1,nbit=d=0; d < 10; d++)
   {
      const unsigned int nvbits = nbits[d];

      if (!nvbits)
         fprintf(fpout, "#define ATL_GetVW%s%s(id_) ATL_VW%s_MIN_%s\n",
                 VN, nms[d], VN, nms[d]);
      else
      {
         fprintf(fpout,"#define ATL_GetVW%s%s(id_) \\\n   ", VN, nms[d]);
         if (nbit+nvbits > (ATL_ISIZE<<3))
         {
            i++;
            nbit = 0;
         }
         fprintf(fpout, "((((ATL_VIEW_%s[ATL_VW%s_IDXMUL(id_)+%u])>>%u)&0x%x)+"
                 "ATL_VW%s_MIN_%s)\n", VN, VN, i, nbit, (1<<nvbits)-1,
                 VN, nms[d]);
      }
      nbit += nvbits;
   }
   fprintf(fpout, "#define ATL_GetVW%sInfo(id_, spf_, m_, n_, k_, amm_, "
           "A2b_, B2b_, b2C_, \\\n"
           "                        b2A_, b2B_, C2b_) \\\n{ \\\n", VN);
   fprintf(fpout,"   const unsigned int *ip_=ATL_VIEW_%s+ATL_VW%s_IDXMUL(id_); \\\n",
           VN, VN);
   fprintf(fpout,"   unsigned int v_=ip_[1]; \\\n");
   fprintf(fpout, "   spf_ = *((float*)ip_); \\\n");
   for (i=1,nbit=d=0; d < 10; d++)
   {
      const unsigned int nvbits = nbits[d];
      if (nbit+nvbits > (ATL_ISIZE<<3))
      {
         fprintf(fpout, "   v_ = ip_[%u]; \\\n", ++i);
         nbit = 0;
      }
      fprintf(fpout, "   %s_ = ATL_VW%s_MIN_%s", mnms[d], VN, nms[d]);
      if (!nvbits)
         fprintf(fpout, "; \\\n");
      else
      {
         fprintf(fpout, " + (v_ & 0x%x); \\\n", (1<<nvbits)-1);
         fprintf(fpout, "   v_ >>= %u; \\\n", nvbits);
      }
      nbit += nvbits;
   }
   fprintf(fpout, "} /* end ATL_GetViewInfo */\n");

   fprintf(fpout, "/*\n * Percent of max perf macros\n */\n");
   for (i=0; i <= NTHRSH; i++)
   {
      const int cnt=(i < NTHRSH) ? THRSH[i]:100;
      mp = (i < NTHRSH) ? mpT[i] : maxPerf;
      fprintf(fpout, "#define ATL_VW%s_%dLCMU %d\n", VN, cnt,
              ATL_iLCM(ATL_iLCM(mp->mu,mp->nu),mp->ku));
      fprintf(fpout, "#define ATL_VW%s_%dLCMMN %d\n", VN, cnt,
              ATL_iLCM(mp->mu,mp->nu));
      fprintf(fpout, "#define ATL_VW%s_%dMB %d\n", VN, cnt, mp->mbB);
      fprintf(fpout, "#define ATL_VW%s_%dNB %d\n", VN, cnt, mp->nbB);
      fprintf(fpout, "#define ATL_VW%s_%dKB %d\n", VN, cnt, mp->kbB);
      fprintf(fpout, "#define ATL_VW%s_%dIDX %d\n",VN, cnt,
              (i < NTHRSH) ? idxT[i] : idxMax);
   }

   fprintf(fpout, "/*\n * Generic names for first-included view\n */\n");
   fprintf(fpout, "#ifndef ATL_VIEW_TABLE\n");
   fprintf(fpout, "   #define ATL_GetViewInfo ATL_GetVW%sInfo\n", VN);
   fprintf(fpout, "   #define ATL_GetViewC2BLK ATL_GetVW%sC2BLK\n", VN);
   fprintf(fpout, "   #define ATL_GetViewBLK2C ATL_GetVW%sBLK2C\n", VN);
   fprintf(fpout, "   #define ATL_GetViewBLK2B ATL_GetVW%sBLK2B\n", VN);
   fprintf(fpout, "   #define ATL_GetViewB2BLK ATL_GetVW%sB2BLK\n", VN);
   fprintf(fpout, "   #define ATL_GetViewBLK2A ATL_GetVW%sBLK2A\n", VN);
   fprintf(fpout, "   #define ATL_GetViewA2BLK ATL_GetVW%sA2BLK\n", VN);
   fprintf(fpout, "   #define ATL_GetViewAMMI ATL_GetVW%sAMMI\n", VN);
   fprintf(fpout, "   #define ATL_GetViewKB ATL_GetVW%sKB\n", VN);
   fprintf(fpout, "   #define ATL_GetViewNB ATL_GetVW%sNB\n", VN);
   fprintf(fpout, "   #define ATL_GetViewMB ATL_GetVW%sMB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_TABLE ATL_VIEW_%s\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_C2BKI ATL_VW%s_BEST_C2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_C2BKI ATL_VW%s_LAST_C2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_BK2BI ATL_VW%s_BEST_BK2BI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_BK2BI ATL_VW%s_LAST_BK2BI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_BK2AI ATL_VW%s_BEST_BK2AI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_BK2AI ATL_VW%s_LAST_BK2AI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_BK2CI ATL_VW%s_BEST_BK2CI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_BK2CI ATL_VW%s_LAST_BK2CI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_B2BKI ATL_VW%s_BEST_B2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_B2BKI ATL_VW%s_LAST_B2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_A2BKI ATL_VW%s_BEST_A2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_A2BKI ATL_VW%s_LAST_A2BKI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_AMMI ATL_VW%s_BEST_AMMI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_AMMI ATL_VW%s_LAST_AMMI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_KU ATL_VW%s_BEST_KU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_KU ATL_VW%s_LAST_KU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_NU ATL_VW%s_BEST_NU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_NU ATL_VW%s_LAST_NU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_MU ATL_VW%s_BEST_MU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_MU ATL_VW%s_LAST_MU\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_KB ATL_VW%s_BEST_KB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_KB ATL_VW%s_LAST_KB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_NB ATL_VW%s_BEST_NB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_NB ATL_VW%s_LAST_NB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_MB ATL_VW%s_BEST_MB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_MB ATL_VW%s_LAST_MB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_BEST_IDX ATL_VW%s_BEST_IDX\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_LAST_IDX ATL_VW%s_LAST_IDX\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NCASES ATL_VW%s_NCASES\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NINT ATL_VW%s_NINT\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_IDXMUL ATL_VW%s_IDXMUL\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MAX_MB ATL_VW%s_MAX_MB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MAX_NB ATL_VW%s_MAX_NB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MAX_KB ATL_VW%s_MAX_KB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NKVEC ATL_VW%s_NKVEC\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_MB ATL_VW%s_MIN_MB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_MB ATL_VW%s_NBIT_MB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_NB ATL_VW%s_MIN_NB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_NB ATL_VW%s_NBIT_NB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_KB ATL_VW%s_MIN_KB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_KB ATL_VW%s_NBIT_KB\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_AMMI ATL_VW%s_MIN_AMMI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_AMMI ATL_VW%s_NBIT_AMMI\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_C2BLK ATL_VW%s_MIN_C2BLK\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_C2BLK ATL_VW%s_NBIT_C2BLK\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_BLK2B ATL_VW%s_MIN_BLK2B\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_BLK2B ATL_VW%s_NBIT_BLK2B\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_BLK2A ATL_VW%s_MIN_BLK2A\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_BLK2A ATL_VW%s_NBIT_BLK2A\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_BLK2C ATL_VW%s_MIN_BLK2C\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_BLK2C ATL_VW%s_NBIT_BLK2C\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_B2BLK ATL_VW%s_MIN_B2BLK\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_B2BLK ATL_VW%s_NBIT_B2BLK\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_MIN_A2BLK ATL_VW%s_MIN_A2BLK\n", VN);
   fprintf(fpout, "   #define ATL_VIEW_NBIT_A2BLK ATL_VW%s_NBIT_A2BLK\n", VN);
   fprintf(fpout, "#endif\n /* end generic name redef guard */\n");

   fprintf(fpout, "\n#endif /* end multiple inclusion guard */\n");
   fclose(fpout);
   free(fn);
/*<ammtime><mb><nb><kb><kernI><a2blkI><b2blkI><blk2cI><blk2aI><blk2bI><c2blkI>*/
}


int main(int nargs, char **args)
{
   ATL_view_t *vb, *vp;
   char *path, *nm;
   char pre, dir, upr;
   ATL_cpnode_t *cb;
   ATL_mmnode_t *ML;

   vb = GetFlags(nargs, args, &pre, &path);
   if (pre == 'c')
      upr = 's';
   else if (pre == 'z')
      upr = 'd';
   else
      upr = pre;
   ML = ReadMMFileWithPath(upr, "res", "FNLK1.sum");
   assert(ML);

   GenMasterFiles(pre, path);
   cb = ReadCPFileWithPath(pre, "res", "cpyPERF.CPS");
   for (vp=vb; vp; vp = vp->next)
   {
      ATL_mmnode_t *mb, *mp;
      mb = ReadMMFile(vp->fnam);
      GenViewFile(pre, path, vp, mb, (vp->flag & (1<<CPV_SELFK)) ? NULL : ML);
      KillAllMMNodes(mb);
   }

   KillAllViews(vb);
   KillAllCPNodes(cb);
   KillAllMMNodes(ML);
   free(path);
   return(0);
}
