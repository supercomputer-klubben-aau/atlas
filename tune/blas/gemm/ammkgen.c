/*
 * Takes all amm kernels to be used for anything in ATLAS, and build master
 * index of kernel structures.  This contains no performance info; that is
 * handled by views that reference these structures
 */
#include "atlas_mmgen.h"
#include "atlas_type.h"

void prepSyrkAmm(char pre, ATL_mmnode_t *kb)
/*
 * This routine returns a new list with the SYRK kernels designated with
 * one being called the non-square case, and one the square.  The non-square
 * case will be designated by setting MMF_RIGHT.  The nonsquare case is the
 * syrk kernel with mu != nu, and if all syrk cases are mu=nu, then the
 * second routine of each type will be set to nonsquare.
 */
{
   ATL_mmnode_t *p, *rs=NULL, *cs=NULL, *ru=NULL, *cu=NULL;
   for (p=kb; p; p = p->next)
   {
      if (p->blask == ATL_KSYRK)
      {
         if (!rs)
            rs = p;
         else if (!cs)
            cs = p;
         else if (!ru)
            ru = p;
         else if (!cu)
            cu = p;
         else
            assert(0);  /* expect only 4 kernels max! */
      }
   }
   if (!rs)    /* if no SYRK kernels, just return */
      return;
   assert(cs && ru && cu);   /* expect all 4 kernels otherwise! */
   assert(rs->mu == rs->nu); /* make sure first case is square */
   assert(cs->mu == cs->nu); /* make sure cplx case is square */
   if (FLAG_IS_SET(rs->flag, MMF_COMPLEX))
   {
      ATL_mmnode_t *tp=rs;
      rs = cs;
      cs = tp;
   }
   if (FLAG_IS_SET(ru->flag, MMF_COMPLEX))
   {
      ATL_mmnode_t *tp=ru;
      ru = cu;
      cu = tp;
   }
   assert((cs->flag&(1<<MMF_COMPLEX))&&(cu->flag&(1<<MMF_COMPLEX)));
   ru->flag |= 1<<MMF_RIGHT;
   cu->flag |= 1<<MMF_RIGHT;
}

ATL_cpnode_t *GetKernCopies(char pre, char *outd, ATL_mmnode_t *kb)
/*
 * Creates a copy node for each amm-incompatible copy needed by kerns in kb
 * RETURNS: queue of kern-spec copies to gen & compile
 */
{
   ATL_cpnode_t *cb=NULL;
   ATL_mmnode_t *mp;
   char *fnd, *fn;
   int dlen;
   dlen = strlen(outd) + 1;
   fnd = malloc(dlen+24);
   assert(fnd);
   fn = fnd + dlen;
   sprintf(fnd, "%s/ATL_XssSyrkIntoC_aXbX.c", outd);
   for (mp=kb; mp; mp = mp->next)
   {
      ATL_cpnode_t *cp;
      int cpflag, ib, ia;
      char cpr = pre;
      char bn[4] = {'1', 'N', 'X', '0'};
      cpflag = (pre == 's' || pre == 'c') ? (1<<CPF_SINGLE) : 0;
      if (FLAG_IS_SET(mp->flag, MMF_COMPLEX))
         cpr = (pre == 's') ? 'c' : 'z';
      else
         cpflag |= 1 << CPF_REAL;
      switch(mp->blask)
      {
      case ATL_KSYRK:  /* SYRK needs C copies for all alpha,beta */
         cpflag |= (1<<CPF_SYRK) | (1<<CPF_CBLK);
         cp = GetCopyNodeFromMM(cpflag, mp, NULL);
/*
 *       Node for each (alpha,beta) case, ->rout=kern name
 */
         for (ib=0; ib < 4; ib++)
         {
            for (ia=0; ia < 3; ia++)
            {
               ATL_cpnode_t *p;
               p = CloneCPNode(cp);
               p->flag |= (1<<(CPF_BE1+ib))|(1<<(CPF_AL1+ia));
               p->rout = fnd;
               assert(p->rout);
               fn[4] = cpr;
               if (mp->mu == mp->nu && !FLAG_IS_SET(mp->flag, MMF_RIGHT))
               {
                  fn[5] = 's';
                  fn[6] = 'q';
               }
               else
               {
                  fn[5] = 'u';
                  fn[6] = 'm';
               }
               fn[18] = bn[ia];
               fn[20] = bn[ib];
               p->genstr = GetCopyGenStr(p);
               fn[21] = '\0';
               p->rout = DupString(fn);
               fn[21] = '.';
               p->next = cb;
               cb = p;
            }
         }
         KillCPNode(cp);
         break;
      default:
         fprintf(stderr, "UNKNOWN KERN TYPE %d!\n", mp->blask);
         assert(0);
      }
   }
   free(fnd);
   return(cb);
}

void GenSyrkPerfH(char pre, char *outd, ATL_mmnode_t *sq, ATL_mmnode_t *um)
{
   FILE *fp;

   fp = OpenMMGenHeader(outd, 0, pre, NULL, "amm", "syrkPerf", NULL);
   fprintf(fp, "   #define ATL_sqsyrkMF %e\n", sq->mflop[0]);
   if (um)
   {
      double syMF=*um->mflop, mmMF=um->mflop[1];
      fprintf(fp, "   #define ATL_umsyrkMF %e\n", syMF);
      fprintf(fp, "   #define ATL_umgemmMF %e\n", mmMF);
      fprintf(fp, "   #define ATL_umratio %f\n", syMF/mmMF);
   }
   CloseGenHeader(fp);
}
void GenSyrkH(char pre, char *outd, ATL_mmnode_t *mb, unsigned int bv)
/*
 * Generate syrk-specific header file
 */
{
   ATL_cpnode_t *cb;
   FILE *fp;
   const int RCSAME=bv&1, UM=bv&2;
   const int ntr = (pre == 'd' || pre == 's') ? 2 : 4;
   const char *sh = (UM) ? "um" : "sq";
   const char *SH = (UM) ? "UM" : "SQ";
   const char *ums = (UM) ? "UM":"";
   int ial, ibe, itr, flg;
   const char be[3] = {'0', '1', 'n'};
   char bes[4] = {'1', 'N', 'X', '0'};
   char trs[4] = {'N', 'T', 'C', 'H'};
   char rcs[4] = {'c', 'r', 'c', 'r'};
   char cjs[4] = {' ', ' ', 'C', 'C'};
   char *sfx[2] = {"", "_L2UT"};
   char *hfx[2] = {"", "_L2UH"};

/*
 * um may be same as sq, in which case we gen a file using um's info
 */
   if (!mb)
   {
      if (UM)
      {
         fp = OpenMMGenHeader(outd, 0, pre, NULL, "amm", "umsyrk", NULL);
         fprintf(fp, "   #include atlas_%camm_sqsyrk.h", pre);
         fprintf(fp, "#define ATL_UMSYRKK_NU ATL_SQSYRKK_NU\n");
         fprintf(fp, "#define ATL_UMSYRKK_MU ATL_SQSYRKK_MU\n");
         fprintf(fp, "#define ATL_UMSYRKK_KVEC ATL_SQSYRKK_KVEC\n");
         fprintf(fp, "#define ATL_UMSYRKK_VLEN ATL_SQSYRKK_VLEN\n");
         fprintf(fp, "#define ATL_%cumsyrkK_b0 ATL_%csqsyrkK_b0\n", pre, pre);
         fprintf(fp, "#define ATL_%cumsyrkK_b1 ATL_%csqsyrkK_b1\n", pre, pre);
         fprintf(fp, "#define ATL_%cumsyrkK_bn ATL_%csqsyrkK_bn\n", pre, pre);
         fprintf(fp, "#define ATL_%ca2blk_umsyrkN ATL_%ca2blk_sqsyrkN\n",
                 pre, pre);
/*
 *       redefine C cpy routs: ATL_<pre>[sq,um]SyrkIntoC_a<alp>_b<bet>
 */
         for (ial=0; ial < 3; ial++)
         {
            for (ibe=0; ibe < 4; ibe++)
            {
               int is;
               for (is=0; is < 2; is++)
               {
                  fprintf(fp, "#ifdef Conj_\n");
                  fprintf(fp, "   #define ATL_%cumHerkIntoC_a%cb%c%s "
                          "ATL_%csqHerkIntoC_a%cb%c%s\n",
                           pre, bes[ial], bes[ibe], sfx[is],
                           pre, bes[ial], bes[ibe], hfx[is]);
                  fprintf(fp, "#endif\n");
                  fprintf(fp, "#define ATL_%cumSyrkIntoC_a%cb%c%s "
                          "ATL_%csqSyrkIntoC_a%cb%c%s",
                          pre, bes[ial], bes[ibe], sfx[is],
                          pre, bes[ial], bes[ibe], sfx[is]
                         );
               }
            }
         }
         CloseGenHeader(fp);
      }
      return;
   }
   fp = OpenMMGenHeader(outd, 0, pre, NULL, "amm", UM?"umsyrk":"sqsyrk", NULL);
   fprintf(fp, "#define ATL_%sSYRKK_VLEN %d\n", SH, mb->vlen);
   fprintf(fp, "#define ATL_%sSYRKK_KVEC %d\n", SH,
           FLAG_IS_SET(mb->flag, MMF_KVEC) ? mb->vlen:0);
   if (UM)
      fprintf(fp, "#define ATL_%sSYRKK_MU %d\n", SH, mb->mu);
   fprintf(fp, "#define ATL_%sSYRKK_NU %d\n", SH, mb->nu);
   fprintf(fp, "#define ATL_%sSYRKK_KU %d\n", SH, mb->ku);

   fprintf(fp, "/*\n * Generic names for syrk info\n */\n");
   fprintf(fp, "#ifndef ATL_SYRKK_KU\n  #define ATL_SYRKK_KU "
           "ATL_%sSYRKK_KU\n#endif\n", SH);
   fprintf(fp, "#ifndef ATL_SYRKK_NU\n  #define ATL_SYRKK_NU "
           "ATL_%sSYRKK_NU\n#endif\n", SH);
   fprintf(fp, "#ifndef ATL_SYRKK_KVEC\n  #define ATL_SYRKK_KVEC "
           "ATL_%sSYRKK_KVEC\n#endif\n", SH);
   fprintf(fp, "#ifndef ATL_SYRKK_VLEN\n  #define ATL_SYRKK_VLEN "
           "ATL_%sSYRKK_VLEN\n#endif\n", SH);
   fprintf(fp, "#ifndef ATL_SYRKK_MU\n  #define ATL_SYRKK_MU "
           "ATL_%sSYRKK_%s\n#endif\n", SH, UM?"MU":"NU");
/*
 * Prototype ATL_<pre>_[sq,um]syrkK kernel
 */
   if (RCSAME && (pre == 'c' || pre == 'z'))
   {
      const char upr = (pre == 'z') ? 'd' : 's';
      for (ibe=0; ibe < 3; ibe++)
      {
         fprintf(fp, "#define ATL_%c%ssyrkK_b%c ATL_%c%ssyrkK_b%c\n",
                 pre, sh, be[ibe], upr, sh, be[ibe]);
      }
   }
   for (ibe=0; ibe < 3; ibe++)
   {
      fprintf(fp,
"void ATL_%c%ssyrkK_b%c(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,\n",
              pre, sh, be[ibe]);
      fprintf(fp,
      "                     TYPE*, const TYPE*, const TYPE*, const TYPE*);\n");
      fprintf(fp, "#ifndef ATL_%camsyrkK_b%c\n   "
                  "#define ATL_%camsyrkK_b%c ATL_%c%ssyrkK_b%c\n#endif\n",
               pre, be[ibe], pre, be[ibe], pre, sh, be[ibe]);
   }
/*
 * Prototype/rename all A cpy routs: ATL_<pre>cm2am_syrk<TA>
 * UM case uses gemm's A/B copy, so no need for UM
 */
   if (!UM)
   {
      flg = (1<<CPF_TOBLK)|(1<<CPF_AL1);
      flg |= (pre == 'd' || pre == 's') ? (1<<CPF_REAL):0;
      flg |= (pre == 'c' || pre == 's') ? (1<<CPF_SINGLE):0;
      cb = GetCopyNodeFromMM(flg, mb, NULL);
      cb->rout = GetCopyName(cb, 0);
      fprintf(fp, "#define ATL_%ca2blk_%ssyrkN %s\n",
              pre, sh, cb->rout);
      free(cb->rout);
      cb->flag |= 1<<CPF_TRANS;
      cb->rout = GetCopyName(cb, 0);
      fprintf(fp, "#define ATL_%ca2blk_%ssyrkT %s\n",
              pre, sh, cb->rout);
      KillAllCopyNodes(cb);
      for (itr=0; itr < 2; itr++)
      {
         char tr = itr ? 'T' : 'N';
         fprintf(fp, "void ATL_%ca2blk_%ssyrk%c", pre, sh, tr);
         if (pre == 'd' || pre == 's')
            fprintf(fp,
            "(ATL_CSZT,ATL_CSZT,const TYPE,const TYPE*,ATL_CSZT,TYPE*);\n");
         else
            fprintf(fp, "(ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,"
                        "ATL_CSZT,TYPE*,TYPE*);\n");
         fprintf(fp, "#ifndef ATL_%ca2blk_syrk%c\n   #define "
         "ATL_%ca2blk_syrk%c ATL_%ca2blk_%ssyrk%c\n#endif\n",
            pre, tr,  pre, tr,  pre, sh, tr);
      }
   }
/*
 * Prototype all C cpy routs: ATL_<pre>[sq,um]SyrkIntoC_a<alp>_b<bet>
 */
   for (ial=0; ial < 3; ial++)
   {
      for (ibe=0; ibe < 4; ibe++)
      {
         int is;
         for (is=0; is < 2; is++)
         {
            fprintf(fp, "#ifdef Conj_\n");
            fprintf(fp, "#define ATL_%c%sSyrkIntoC_a%cb%c%s "
                    "ATL_%c%sHerkIntoC_a%cb%c%s\n",
                     pre, sh, bes[ial], bes[ibe], sfx[is],
                     pre, sh, bes[ial], bes[ibe], hfx[is]);
            fprintf(fp, "#endif\n");
            fprintf(fp, "void ATL_%c%sSyrkIntoC_a%cb%c%s\n", pre, sh, bes[ial],
                    bes[ibe], sfx[is]);
            if (pre == 's' || pre == 'd')
               fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR, TYPE *,ATL_CSZT);\n");
            else
               fprintf(fp, "   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,"
                       "const TYPE*,const SCALAR, TYPE *,ATL_CSZT);\n");

            fprintf(fp, "#ifndef ATL_%cSyrkIntoC_a%cb%c%s\n   #define "
            "ATL_%cSyrkIntoC_a%cb%c%s ATL_%c%sSyrkIntoC_a%cb%c%s\n#endif\n",
               pre, bes[ial], bes[ibe], sfx[is],
               pre, bes[ial], bes[ibe], sfx[is],
               pre, sh, bes[ial], bes[ibe], sfx[is]);
         }
      }
   }

   CloseGenHeader(fp);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   fprintf(stderr,
"This routine creates master list of kernel IDs to use, and generates them.\n"
"It expects all input kerns in [s,d]FNLK1.sum, with K-cleaners in ivar.\n");
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr,"USAGE: %s [flags:\n", name);
   fprintf(stderr, "   -b <blasK.sum>: non-gemm kerns file; can be repeated\n");
   fprintf(stderr, "   -o <path>: what directory to output all files to?\n");
   fprintf(stderr, "   -p [s,d]: set type/precision prefix (d) \n"
           "      s/d will generate for complex (c/z) as well\n");

   exit(ierr ? ierr : -1);
}

ATL_mmnode_t *GetFlags
   (int nargs, char **args, char *PRE, char **OUTD, ATL_mmnode_t **BKB)
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
        if (pre == 'z')
           pre = 'd';
        else if (pre == 'c')
           pre = 's';
        break;
      case 'o':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         outd = DupString(args[i]);
         break;
      case 'i':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         if (mb)
            KillAllMMNodes(mb);
         mb = ReadMMFile(args[i]);
         break;
      case 'b': /* -b <blasK.sum> */
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         mp = ReadMMFile(args[i]);
         assert(mp);
         bkb = ATL_JoinMMQs(bkb, mp);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!mb)
   {
      mb = ReadMMFileWithPath(pre, "res", "FNLK1.sum");
      assert(mb);
   }
   if (!outd)
      outd = DupString("tmp");

   *BKB = bkb;
   *OUTD = outd;
   *PRE = pre;
   return(mb);
}

int numBitsNeeded(int N)
{
   int i;
   if (N <= 0)
      return(0);
   for (i=0; (1<<i) <= N; i++);
   return(i);
}
#define NFLGBITS 2  /* KRUNTIME, KMAJOR */
#define NIVAR 8
int *getBitsNeeded(ATL_mmnode_t *mb, int *NINT)
/*
 * 1st NIVAR entries has number of bits required after subtracting MINVALs
 * which are stoted in 2nd NIVAR entries
 */
{
   ATL_mmnode_t *mp;
   unsigned int cnt=1, i, nint, nbits, pbits;
   static unsigned int ib[NIVAR*2];

   assert(mb);
   ib[0] = NFLGBITS;
   ib[NIVAR] = 0;
   ib[1] = ib[1+NIVAR] = mb->mu;
   ib[2] = ib[2+NIVAR] = mb->nu;
   ib[3] = ib[3+NIVAR] = mb->ku;
   ib[4] = ib[4+NIVAR] = mb->vlen;
   ib[5] = ib[5+NIVAR] = mb->kbmin;
   ib[6] = ib[6+NIVAR] = mb->kbmax;
   ib[7] = ib[7+NIVAR] = mb->ivar;
   for (mp=mb->next; mp; mp=mp->next)
   {
      cnt++;
      ib[1] = Mmax(ib[1],mp->mu);
      ib[1+NIVAR] = Mmin(ib[1+NIVAR],mp->mu);
      ib[2] = Mmax(ib[2],mp->nu);
      ib[2+NIVAR] = Mmin(ib[2+NIVAR],mp->nu);
      ib[3] = Mmax(ib[3],mp->ku);
      ib[3+NIVAR] = Mmin(ib[3+NIVAR],mp->ku);
      ib[4] = Mmax(ib[4],mp->vlen);
      ib[4+NIVAR] = Mmin(ib[4+NIVAR],mp->vlen);
      ib[5] = Mmax(ib[5],mp->kbmin);
      ib[5+NIVAR] = Mmin(ib[5+NIVAR],mp->kbmin);
      ib[6] = Mmax(ib[6],mp->kbmax);
      ib[6+NIVAR] = Mmin(ib[6+NIVAR],mp->kbmax);
      ib[7] = Mmax(ib[7],mp->ivar);
      ib[7+NIVAR] = Mmin(ib[7+NIVAR],mp->ivar);
   }
   nint=1;
   nbits = NFLGBITS;
   pbits = ATL_PSIZE<<3;
   for (i=1; i < NIVAR; i++)
   {
      int bits;
      ib[i] = bits = numBitsNeeded(ib[i]-ib[i+NIVAR]);
      nbits += bits;
      if (nbits > pbits) /* not enough room left in int */
      {
         nbits = bits;
         nint++;
      }
   }
   assert(nint <= 3);
   *NINT = nint;
   return(ib);
}

int PrintIntsHex(FILE *fp, ATL_mmnode_t *mp, int nint, int *nbits)
/* RETURNS: flag stored */
{
   ATL_iptr_t ib, ish;
   int flag;

   ib = FLAG_IS_SET(mp->flag, MMF_KRUNTIME);
   ib |= FLAG_IS_SET(mp->flag, MMF_KVEC) ? 2 : 0;
   flag = ib;
   ish = nbits[0];
   if (ish+nbits[1] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->mu-nbits[1+NIVAR])) << ish;
   ish += nbits[1];
   if (ish+nbits[2] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->nu-nbits[2+NIVAR])) << ish;
   ish += nbits[2];
   if (ish+nbits[3] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->ku-nbits[3+NIVAR])) << ish;
   ish += nbits[3];
   if (ish+nbits[4] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->vlen-nbits[4+NIVAR])) << ish;
   ish += nbits[4];
   if (ish+nbits[5] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->kbmin-nbits[5+NIVAR])) << ish;
   ish += nbits[5];
   if (ish+nbits[6] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->kbmax-nbits[6+NIVAR])) << ish;
   ish += nbits[6];
   if (ish+nbits[7] > ATL_PSIZE*8)
   {
      fprintf(fp, ", 0x%lxL", ib);
      ib = ish = 0;
   }
   ib |= ((ATL_iptr_t)(mp->ivar-nbits[7+NIVAR])) << ish;
   ish += nbits[7];
   if (ish)
      fprintf(fp, ", 0x%lxL", ib);
   return(flag);
}

void GenKernArr(char pre, char *outd, ATL_mmnode_t *mb)
{
   FILE *fpout;
   ATL_mmnode_t *mp;
   char *VN[8]={"FLAG", "MU", "NU", "KU", "VLEN", "KBMIN", "KBMAX","K1IDX"};
   unsigned int *nbits;
   unsigned int k, i, nint, N, nkvec;
   ATL_iptr_t ibt;

   assert(sizeof(ATL_iptr_t) == ATL_psize);
   N = ATL_CountNumberOfMMNodes(mb);
   fpout = OpenMMGenHeader(outd, 0, pre, "amm", NULL, "kern", mb);
   fprintf(fpout, "#include \"atlas_%camm_proto.h\"\n\n", pre);

   nbits = getBitsNeeded(mb, &nint);
   k = nint + 3;
   fprintf(fpout, "#define ATL_ENTRYN %u\n", k);
   for (i=0; i < 32 && !((1<<i)&k); i++);
   if (k^(1<<i) == 0)
      fprintf(fpout, "#define ATL_AMM_Idx2Entry(i_) ((i_)<<%u)\n", i);
   else if (k == 5)
      fprintf(fpout, "#define ATL_AMM_Idx2Entry(i_) (((i_)<<2)+(i_))\n");
   else if (k == 6)
      fprintf(fpout, "#define ATL_AMM_Idx2Entry(i_) (((i_)<<2)+((i_)<<1))\n");
   else
      fprintf(fpout, "#define ATL_AMM_Idx2Entry(i_) ((i_)*%u)\n", k);
   fprintf(fpout, "#define ATL_NINTPACK %u\n\n", nint);
   for (i=0; i < NIVAR; i++)
   {
      fprintf(fpout, "#define ATL_AMM_MIN_%s  %u\n", VN[i], nbits[i+NIVAR]);
      fprintf(fpout, "#define ATL_AMM_NBIT_%s %u\n", VN[i], nbits[i]);
   }
   fprintf(fpout, "\n");
   fprintf(fpout, "#define ATL_AMM_ML ATL_%cAMM_ML\n", pre);
   fprintf(fpout,
           "#ifndef ATL_DECL_\n   extern const ATL_iptr_t ATL_AMM_ML[%u];\n",
           N*(3+nint));
   fprintf(fpout, "#else\nconst ATL_iptr_t ATL_AMM_ML[%u]=\n{\n", N*(3+nint));

   mp = mb;
   fprintf(fpout, "(ATL_iptr_t)%s_b1,\n(ATL_iptr_t)%s_b0,\n(ATL_iptr_t)%s_bn\n",
           mp->auth, mp->auth, mp->auth);
   i = PrintIntsHex(fpout, mp, nint, nbits);
   fprintf(fpout, "/* %2u: FL=%x, U=(%u,%u,%u), VL=%u, KB=(%u,%u) KclnI=%d */",
           0, i, mp->mu, mp->nu, mp->ku, mp->vlen, mp->kbmin, mp->kbmax,
           mp->ivar-1);
   nkvec = FLAG_IS_SET(mp->flag, MMF_KVEC);
   for (i=1,mp=mb->next; mp; mp = mp->next, i++)
   {
      fprintf(fpout,
              ",\n(ATL_iptr_t)%s_b1,\n(ATL_iptr_t)%s_b0,\n(ATL_iptr_t)%s_bn\n",
              mp->auth, mp->auth, mp->auth);
      k = PrintIntsHex(fpout, mp, nint, nbits);
      fprintf(fpout,
              "/* %2u: FL=%x, U=(%u,%u,%u), VL=%u, KB=(%u,%u), KclnI=%d */",
              i, k, mp->mu, mp->nu, mp->ku, mp->vlen, mp->kbmin, mp->kbmax,
              mp->ivar-1);
      nkvec += FLAG_IS_SET(mp->flag, MMF_KVEC);
   }
   fprintf(fpout, "\n};\n#endif /* end else of ATL_DECL */\n\n");


   if (nbits[0] == 0)
      fprintf(fpout, "#define ATL_AMM_GetFLAG(idx_) ATL_AMM_MIN_FLAG\n");
   else
   {
      fprintf(fpout, "#define ATL_AMM_GetFLAG(idx_) \\\n");
      fprintf(fpout,
              "   (((ATL_AMM_ML[ATL_AMM_Idx2Entry(idx_)+3])&0x%x)",
              (1<<(nbits[0]))-1);
      if (nbits[0+NIVAR])
         fprintf(fpout, "+ATL_AMM_MIN_FLAG)\n");
      else
         fprintf(fpout, ")\n");
   }
   if (nbits[3] == 0)
      fprintf(fpout, "#define ATL_AMM_GetKU(idx_) ATL_AMM_MIN_KU\n");
   else
   {
      fprintf(fpout, "#define ATL_AMM_GetKU(idx_) \\\n");
      if (nbits[0]+nbits[1]+nbits[2]+nbits[3] <= sizeof(ATL_iptr_t)*8)
         fprintf(fpout,
                 "   ((((ATL_AMM_ML[ATL_AMM_Idx2Entry(idx_)+3])>>%u)&0x%x)",
                 nbits[0]+nbits[1]+nbits[2], (1<<(nbits[3]))-1);
      else
      {
         assert(nbits[0]+nbits[1]+nbits[2] <= sizeof(ATL_iptr_t)*8);
         fprintf(fpout, "   (((ATL_AMM_ML[ATL_AMM_Idx2Entry(idx_)+4])&0x%x)",
                 (1<<(nbits[3]))-1);
      }
      if (nbits[3+NIVAR])
         fprintf(fpout, "+ATL_AMM_MIN_KU)\n");
      else
         fprintf(fpout, ")\n");
   }
   fprintf(fpout, "#define ATL_AMM_GetMNU(idx_, MU_, NU_) \\\n{\\\n");
   fprintf(fpout,
      "   ATL_iptr_t v_ = (ATL_AMM_ML[ATL_AMM_Idx2Entry(idx_)+3])>>%u; \\\n",
            nbits[0]);
   if (nbits[1])
   {
      fprintf(fpout, "   MU_ = ATL_AMM_MIN_MU + (v_ & 0x%x); \\\n",
              (1<<nbits[1])-1);
      fprintf(fpout, "   v_ >>= %u; \\\n", nbits[1]);
   }
   else
      fprintf(fpout, "   MU_ = ATL_AMM_MIN_MU; \\\n");
   if (nbits[2])
      fprintf(fpout, "   NU_ = ATL_AMM_MIN_NU + (v_ & 0x%x); \\\n",
              (1<<nbits[2])-1);
   else
      fprintf(fpout, "   NU_ = ATL_AMM_MIN_NU; \\\n");
   fprintf(fpout, "}\n\n");

   fprintf(fpout, "#define ATL_AMM_iinfo(ia_,FLAG_,MU_,NU_,KU_,VLEN_,"
                  "KBMIN_,KBMAX_,K1IDX_) \\\n");
   fprintf(fpout, "{ \\\n   ATL_iptr_t v_ = *(ia_); \\\n");
   for (nint=k=i=0; i < 8; i++)
   {
      int nb=nbits[i];
      if (nb)
      {
         if (k + nb > ATL_PSIZE*8)
         {
            fprintf(fpout, "   v_ = (ia_)[%u]; \\\n", ++nint);
            k = 0;
         }
         if (nbits[i+NIVAR])
            fprintf(fpout, "   %s_ = ATL_AMM_MIN_%s + (v_&0x%x); \\\n",
                    VN[i], VN[i], (1<<nb)-1);
         else
            fprintf(fpout, "   %s_ = v_ & 0x%x; \\\n", VN[i], (1<<nb)-1);
         if (i != 7)
            fprintf(fpout, "   v_ >>= %u; \\\n", nb);
      }
      else
         fprintf(fpout, "   %s_ = ATL_AMM_MIN_%s; \\\n", VN[i], VN[i]);
   }
   fprintf(fpout, "} /* end ATL_AMM_iinfo definition */\n");

   fprintf(fpout, "\n#define ATL_AMM_kinfo(id_, mm0_, mm1_, mmN_) \\\n");
   fprintf(fpout, "{ \\\n   const ATL_iptr_t *ip_=ATL_AMM_ML + "
                  "ATL_AMM_Idx2Entry(id_); \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mm1_)) = *ip_; \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mm0_)) = ip_[1]; \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mmN_)) = ip_[2]; \\\n");
   fprintf(fpout, "} /* end ATL_AMM_kinfo definition */\n");

   fprintf(fpout, "\n#define ATL_AMM_info(id_,mm0_, mm1_, mmN_,FLAG_,"
                  "MU_,NU_,KU_,VLEN_,KBMIN_,KBMAX_,K1IDX_) \\\n");
   fprintf(fpout, "{ \\\n   const ATL_iptr_t *ip_=ATL_AMM_ML + "
                  "ATL_AMM_Idx2Entry(id_); \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mm1_)) = *ip_; \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mm0_)) = ip_[1]; \\\n");
   fprintf(fpout, "   *((ATL_iptr_t*)(mmN_)) = ip_[2]; \\\n");
   fprintf(fpout, "   ATL_AMM_iinfo(ip_+3,FLAG_,MU_,NU_,KU_,VLEN_,"
                  "KBMIN_,KBMAX_,K1IDX_); \\\n");
   fprintf(fpout, "} /* end ATL_AMM_info definition */\n");

   fprintf(fpout,
           "/*\n * Following macros return 1 if flag set, else 0\n */\n");
   fprintf(fpout, "#define ATL_AMM_KRUNTIME(flg_) ((flg_)&1)\n");
   fprintf(fpout, "#define ATL_AMM_KMAJOR(flg_) (((flg_)>>1)&1)\n");
   CloseGenHeader(fpout);
}
/*
 * generate:
 * ATL_AmmKernbyIdx(idx_, k0_,k1_,kn_, cpidx, flg_, vlen_, mu_, nu_, ku_,
 *                  kbmin_, kbmax_)
 */
void GenKernInf(ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;
   FILE *fpout=NULL;
   unsigned int MAXs[8];
   char *nms[8] = {"flg", "vlen", "mu", "nu", "ku", "kbmin", "kbmax"};
   unsigned int MU, NU, KU, VLEN, KBMAX, KBMIN, KCleanIdx;
   unsigned int nbits, pbits, nint, iv, i;

   assert(mb);
   MAXs[0] = NFLGBITS;
   MAXs[1] = mb->mu;
   MAXs[2] = mb->nu;
   MAXs[3] = mb->ku;
   MAXs[4] = mb->vlen;
   MAXs[5] = mb->kbmin;
   MAXs[6] = mb->kbmax;
   for (mp=mb->next; mp; mp=mp->next)
   {
      MAXs[0] = Mmax(MAXs[0],mp->mu);
      MAXs[1] = Mmax(MAXs[1],mp->nu);
      MAXs[2] = Mmax(MAXs[2],mp->ku);
      MAXs[3] = Mmax(MAXs[3],mp->vlen);
      MAXs[4] = Mmax(MAXs[4],mp->kbmin);
      MAXs[5] = Mmax(MAXs[5],mp->kbmax);
   }
   nint=1;
   nbits = NFLGBITS;
   pbits = ATL_PSIZE<<3;
   for (i=1; i < 7; i++)
   {
      int bits;
      MAXs[i] = bits = numBitsNeeded(MAXs[i]);
      nbits += bits;
      if (nbits > pbits) /* not enough room left in int */
      {
         nbits = bits;
         nint++;
      }
   }
   assert(nint <= 3);
   fprintf(fpout, "ATL_AmmKernbyIdx(idx_, k0_, k1_, kn_, flg_, vlen_, "
           "mu_, nu_, ku_, kbmin_, kbmax_)\\\n{\\\n");
   if (nint = 1)
      fprintf(fpout, "const ATL_iptr_t id_=(idx_)<<1; \\\n");
   else if (nint = 2)
      fprintf(fpout, "const ATL_iptr_t id_=(idx_)+((idx_)<<1);\\\n");
   else if (nint = 3)
      fprintf(fpout, "const ATL_iptr_t id_=(idx_)<<2;\\\n");
   fprintf(fpout, "ATL_iptr_t bits_;\\\n");
   fprintf(fpout, "   (k0_) = (void*) (ATL_AMM_KERN[id_]);\n");
   fprintf(fpout, "   (k1_) = (void*) (ATL_AMM_KERN[id_+1]);\n");
   fprintf(fpout, "   (kn_) = (void*) (ATL_AMM_KERN[id_+2]);\n");
   fprintf(fpout, "   bits_ = ATL_AMM_KERN[id_+3];\\\n");
   fprintf(fpout, "   flg_ = bits_ & 0x%x;\\\n", 1<<(NFLGBITS+1)-1);
   nbits = NFLGBITS;
   nint=3;
   for (i=1; i < 7; i++)
   {
      int bits = MAXs[iv];
      if (nbits + bits > pbits)
      {
         fprintf(fpout, "bits_ = ATL_AMM_KERN[id_+%d];\\\n", ++nint);
         nbits = 0;
      }
      if (!nbits)
         fprintf(fpout, "   %s_ = bits_ & 0x%x;\\\n", nms[i], (1<<(bits+1))-1);
      else
         fprintf(fpout, "   %s_ = ((bits_)>>%u) & 0x%x;\\\n", nms[i],
                 nbits, (1<<(bits+1))-1);
      nbits += bits;
   }
   fprintf(fpout, "} /* end ATL_AmmKernByIdx macro def */\n");
}

void GenKerns(char pre, char *outd, ATL_mmnode_t *mb, ATL_mmnode_t *bkb,
              ATL_cpnode_t *cb)
{
   ATL_mmnode_t *ub, *mp;
   ATL_cpnode_t *cp;
   char *sgen=NULL;
   int i, len, dlen;
   char pr=pre;

   dlen = strlen(outd);
   len = 0;
/*
 * mb has unique compilations, reduce this to unique files, and then generate
 */
   ub = AddUniquePerfMMKernsToList(NULL, mb);
   ub = AddUniquePerfMMKernsToList(ub, bkb);
/*
 * Generate/copy all required kernels.
 */
   printf("\nGENERATING AMM KERNS:\n");
   for (cp=cb; cp; cp=cp->next)
   {
      if ( (i=system(cp->genstr)) )
      {
         fprintf(stderr, "GENSTR RETURNS %d:\n'%s'\n", i, cp->genstr);
         exit(i);
      }
   }
   for (mp=ub; mp; mp = mp->next)
   {
      const int id=mp->ID;
      if (!id)
      {
         printf("   -> %s\n", mp->rout);
         assert(mp->genstr);
         if ( (i=system(mp->genstr)) )
         {
            fprintf(stderr, "GENSTR RETURNS %d:\n'%s'\n", i, mp->genstr);
            exit(i);
         }
      }
      else /* user-supplied kernel */
      {
         ATL_mmnode_t *p;
         printf("   %s -> %s\n", mp->genstr, mp->rout);
         i = strlen(mp->genstr) + strlen(mp->rout) + dlen + 16;
         if (i > len)
         {
            if (sgen)
               free(sgen);
            sgen = malloc(i*sizeof(char));
            assert(sgen);
            len = i;
         }
         sprintf(sgen, "cp AMMCASES/%s %s/%s", mp->genstr, outd, mp->rout);

         if ( (i=system(sgen)) )
         {
            fprintf(stderr, "FAILED CP='%s'\n", sgen);
            exit(i);
         }
      }
   }
   printf("DONE GENERATING AMM KERNS.\n");
   free(sgen);
   KillAllMMNodes(ub);
}

void GenMake(char pre, char *outd, ATL_mmnode_t *mb, ATL_mmnode_t *bkb,
             ATL_cpnode_t *cb)
/*
 * mb files have already been made at least compile-time unique (same source
 * file might occur multiple times due to need to compile with -DKB)
 */
{
   FILE *fp;
   char *fn;
   ATL_mmnode_t *mp, *MB;
   ATL_cpnode_t *cp;
   char *sals[3] = {"1", "N1", "X"};
   char als[3] = {'1', 'n', 'X'};
   char *sbes[4] = {"0", "1", "N1", "X"};
   char bes[4] = {'0', '1', 'n', 'X'};  /* use 1st 3 for mmkerns */
   char ctas[4] = {'N', 'T', 'C', 'H'};
   char dcomp[8] = {'$', '(', 'D', 'M', 'C', ')', '\0'};
   char dflags[12] = {'$','(','D','M','C','F','L','A','G','S',')','\0'};


   MB = AddUniqueMMKernCompList(NULL, mb);
   MB = ATL_JoinMMQs(MB, bkb);
   fn = malloc(strlen(outd)+11);
   assert(fn);
   sprintf(fn, "%s/%cMake_amm", outd, pre);
   fp = fopen(fn, "w");
   assert(fp);
   free(fn);
   fprintf(fp, "include ../Make.inc\n");
   fprintf(fp, "CDEFS2=$(CDEFS)\n\n");
/*
 * Spew out kernels to be compiled
 */
   fprintf(fp, "objs =");
   for (mp=MB; mp; mp = mp->next)
   {
      int ib;
      for (ib=0; ib < 3; ib++)
         fprintf(fp, " \\\n       %s_b%c.o", mp->auth, bes[ib]);
   }
   for (cp=cb; cp; cp = cp->next)
   {
      char *sp;
      fprintf(fp, " \\\n       %s.o", cp->rout);
      sp = strstr(cp->rout, "SyrkIntoC");
      if (sp)
      {
         char pre = sp[-3];
         if (pre == 'z' || pre == 'c')
         {
            sp[0] = '\0';
            fprintf(fp, " \\\n       %sHerk%s.o", cp->rout, sp+4);
            sp[0] = 'S';
         }
      }
   }
   fprintf(fp, "\n");
/*
 * library make targets
 */
   fprintf(fp, "\n\nlib : %clib.grd\nall : %clib.grd\n%clib : %clib.grd\n",
           pre, pre, pre, pre);
   fprintf(fp, "%clib.grd : $(objs)\n", pre);
   fprintf(fp, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\t touch %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n", pre);
   fprintf(fp, "%cclean:\n\t- rm -f $(objs)\n", pre);
   fprintf(fp, "killall : %ckillall\n", pre);
   fprintf(fp, "%ckillall : %cclean\n", pre, pre);
   fprintf(fp, "\t- $(ARCHIVER) d $(ATLASlib) $(objs)\n");
   fprintf(fp, "\t $(RANLIB) $(ATLASlib)\n");
   fprintf(fp, "\t- rm -f ATL_%c*.[S,c]\n\n", pre);
/*
 * Make targets for amm kerns
 */
   fprintf(fp, "#\n#  AMM kernel rules\n#\n");
   fn = (pre == 'd' || pre == 'z') ?  "-DDREAL=1" : "-DSREAL";
   for (mp=MB; mp; mp = mp->next)
   {
      int ib;
      char *comp, *flgs;

      comp = GetMMKernComp(mp, dcomp, dflags, &flgs);
      for (ib=0; ib < 3; ib++)
      {
         char *sp=" ";
         fprintf(fp, "%s_b%c.o : %s\n", mp->auth, bes[ib], mp->rout);
         fprintf(fp, "\t%s $(CDEFS2) %s -DBETA%s=1 \\\n", comp, fn, sbes[ib]);
         if (!FLAG_IS_SET(mp->flag, MMF_KRUNTIME))
            fprintf(fp, "        -DMB=%d -DNB=%d -DKB=%d",
                    mp->mbB, mp->nbB, ((mp->kbB+mp->ku-1)/mp->ku)*mp->ku);
         else
            sp = "        ";
         if (FLAG_IS_SET(mp->flag, MMF_MVA))
         {
            fprintf(fp, "%s-DATL_MOVEA", sp);
            sp = " ";
         }
         if (FLAG_IS_SET(mp->flag, MMF_MVB))
         {
            fprintf(fp, "%s-DATL_MOVEB", sp);
            sp = " ";
         }
         if (FLAG_IS_SET(mp->flag, MMF_MVC))
         {
            fprintf(fp, "%s-DATL_MOVEC", sp);
            sp = " ";
         }
         fprintf(fp, " \\\n        %s \\\n", flgs);
         fprintf(fp,
                 "        -DATL_USERMM=%s_b%c \\\n        -c -o %s_b%c.o \\\n",
                 mp->auth, bes[ib], mp->auth, bes[ib]);
         fprintf(fp, "        %s\n", mp->rout);
      }
   }
   fprintf(fp, "#\n#  non-AMM kernel copy rules\n#\n");
   for (cp=cb; cp; cp = cp->next)
   {
      char *ts;
      if (pre == 'd')
        ts = (cp->flag & (1<<CPF_REAL)) ? "DREAL" : "DCPLX";
      else
        ts = (cp->flag & (1<<CPF_REAL)) ? "SREAL" : "SCPLX";
      fprintf(fp, "%s.o : %s.c\n", cp->rout, cp->rout);
      fprintf(fp, "\t$(KC) $(CDEFS2) -D%s=1 -DATL_USERCPMM=%s $(KCFLAGS) \\\n",
              ts, cp->rout);
      fprintf(fp, "        -c %s.c\n", cp->rout);
      if (!(cp->flag & (1<<CPF_REAL)))
      {
         char *sp;
         sp = strstr(cp->rout, "SyrkIntoC");
         if (sp);
         {
            sp[0]='\0';
            fprintf(fp, "%sHerk%s.o : %sSyrk%s.c\n", cp->rout, sp+4,
                    cp->rout, sp+4);
            if (pre == 'd')
               ts = "-DDCPLX=1 -DConj_=1";
            else
               ts = "-DSCPLX=1 -DConj_=1";
            fprintf(fp,
            "\t$(KC) $(CDEFS2) %s -DATL_USERCPMM=%sHerk%s \\\n",
                    ts, cp->rout, sp+4);
            fprintf(fp, "        $(KCFLAGS) -o $@ -c %sSyrk%s.c\n",
                    cp->rout, sp+4);
            sp[0]='S';
         }
      }
   }
   KillAllMMNodes(MB);
   fclose(fp);
}


void GenTrsmALLT(char pre, char *outd, char sd, ATL_mmnode_t *mb)
{
   char *of;
   FILE *fp;
   ATL_mmnode_t *mp;
   unsigned long long flg=0;
   int k, L, n;
   L = strlen(outd) + 24;
   of = malloc(L);
   assert(of);

   for (n=0, mp=mb; mp; n++, mp=mp->next)
   {
      if (mp->TA == AtlasTrans && mp->TB == AtlasTrans)
         flg |= (1L<<n);
   }
   assert((sizeof(long)<<3) >= n);
   k = sprintf(of, "%s/atlas_%ctrsm%cLN_ALLT.h", outd, pre, sd);
   assert(k<L);
   fp = fopen(of, "w");
   assert(fp);
   fprintf(fp, "#ifndef ATLAS_%cTRSM%cLN_ALLT_H\n",
           toupper(pre), toupper(sd));
   fprintf(fp, "   #define ATLAS_%cTRSM%cLN_ALLT_H 1\n",
           toupper(pre), toupper(sd));
   fprintf(fp, "   #define ATL_trsm_allT(i_) ((0x%lxL>>(i_))&1)\n",
           (unsigned long)flg);
   fprintf(fp, "#endif /* end multiple inclusion guard */\n");
   fclose(fp);
   free(of);
}

ATL_mmnode_t *GenAllSyrkH(char pre, char *outd, ATL_mmnode_t **BKB)
{
   ATL_mmnode_t *rp=NULL, *cp=NULL, *rU=NULL, *cU=NULL, *ret, *bkb=(*BKB);
   unsigned int bv;
   char cpr = (pre == 's') ? 'c' : 'z';

   for (bv=0,ret=bkb; ret && ret->blask == ATL_KSYRK; ret = ret->next, bv++);
   assert(bv == 4);
   cU = bkb;
   rU = cU->next;
   cp = rU->next;
   rp = cp->next;
/*
 * In case reordering has messed up cplx,real ordering expectation, fix
 */
   if (!FLAG_IS_SET(cp->flag, MMF_COMPLEX))
   {
      ATL_mmnode_t *tp=rp;
      rp = cp;
      cp = tp;
      assert(FLAG_IS_SET(cp->flag, MMF_COMPLEX));
   }
   if (!FLAG_IS_SET(cU->flag, MMF_COMPLEX))
   {
      ATL_mmnode_t *tp=rU;
      rU = cU;
      cU = tp;
      assert(FLAG_IS_SET(cU->flag, MMF_COMPLEX));
   }
/*
 * In case reordering messed up (UM pair),(SQ pair) ordering, fix
 */
   if (FLAG_IS_SET(rp->flag, MMF_RIGHT))
   {
      ATL_mmnode_t *tp=rp;
      rp = rU;
      rU = tp;
      tp = cp;
      cp = cU;
      cU = tp;
   }
   assert(!FLAG_IS_SET(rp->flag, MMF_RIGHT)&&!FLAG_IS_SET(cp->flag, MMF_RIGHT));
   assert((rp->mu == rp->nu) && (cp->nu == cp->mu));
   assert(rp->blask == ATL_KSYRK);
   assert(!(rU->mu % rU->nu) || !(rU->nu % rU->mu));
   assert(!(cU->mu % cU->nu) || !(cU->nu % cU->mu));
   bv = MMKernsPerfSame(rp, cp);
   GenSyrkH(pre, outd, rp, bv);
   GenSyrkH(cpr, outd, cp, bv);
   GenSyrkPerfH(pre, outd, rp, rU);
   GenSyrkPerfH(cpr, outd, cp, cU);
   if (bv) /* if real & cplx same, take cplx out to avoid double compile */
   {
      bkb = RemoveMMNodeFromQ(bkb, cp);
      KillMMNode(cp);
   }
/*
 * The code for same kernel presently untested!  HERE HERE RCW.
 */
   if (rU)
   {
      unsigned int bvU;
      bvU = MMKernsPerfSame(rU, cU) | 2;
      if (MMKernsPerfSame(rU, rp)) /* is this real kernel same as mu=nu kern? */
      {
         bkb = RemoveMMNodeFromQ(bkb, rU);
         GenSyrkH(pre, outd, NULL, bvU);
         KillMMNode(rU);
         if (MMKernsPerfSame(cU, cp))
         {
            bkb = RemoveMMNodeFromQ(bkb, cU);
            KillMMNode(cU);
            cU = NULL;
         }
         GenSyrkH(cpr, outd, cU, bvU);
      }
      else if (MMKernsPerfSame(cU, cp)) /* is cplx kernel same as mu=nu kern? */
      {
         GenSyrkH(pre, outd, rU, bvU);
         GenSyrkH(cpr, outd, NULL, bvU);
         bkb = RemoveMMNodeFromQ(bkb, cU);
         KillMMNode(cU);
         cU = NULL;
      }
      else
      {
         bvU = MMKernsPerfSame(rU, cU) | 2;
         GenSyrkH(pre, outd, rU, bvU);
         GenSyrkH(cpr, outd, cU, bvU);
         if (bvU & 1)
         {
            bkb = RemoveMMNodeFromQ(bkb, cU);
            KillMMNode(cU);
         }
      }
/*
 *    Overload RIGHT flag to mean this is mu != nu case, in case both SQ & UM
 *    actually use mu=nu kernels!
 */
   }
   *BKB = bkb;
   return(ret);
}

void GenAllFiles(char pre, char *outd, ATL_mmnode_t *mb, ATL_mmnode_t *bkb,
                 ATL_cpnode_t *cb)
{
   FILE *fpout;
   ATL_mmnode_t *mp;
   int ib, inxt, iaut;
   char bes[3] = {'0', '1', 'n'};
   char fn[24];

   mp = bkb;
   while (mp)
   {
      ATL_mmnode_t *nxt = mp->next;
      char cpr = (pre == 's') ? 'c' : 'z';
      switch(mp->blask)
      {
         ATL_mmnode_t *cp;
         int i;
/*
 *    For SYRK, we expect 4 kernels, with each pair being the real & complex
 *    versions of the same kernel, with the first pair being the mu=nu case
 *    used by pure inner-product and cases ip when diagonal C block's A are
 *    stored in different format than gemm's.  The second pair is for when
 *    ipmenUM is used for gemm.
 */
      case ATL_KSYRK:
         i = mp == bkb;
         nxt = GenAllSyrkH(pre, outd, &mp);
         bkb = (i) ? mp : bkb;
         break;
      default:
         assert(0);
      }
      mp =  nxt;
   }
   fpout = OpenMMGenHeader(outd, 0, pre, "amm", NULL, "proto", mb);
   assert(fpout);
   inxt = GetOffset(&mb->next, mb);
   iaut = GetOffset(&mb->auth, mb);
   for (ib=0; ib < 3; ib++)
      PrintMMProtos(fpout, pre, "KERN", mb, iaut, bes[ib]);
   CloseGenHeader(fpout);
   GenKerns(pre, outd, mb, bkb, cb);
   GenMake(pre, outd, mb, bkb, cb);
   GenKernArr(pre, outd, mb);
}

int main(int nargs, char **args)
{
   char pre;
   char *outd;
   ATL_mmnode_t *mb, *ub, *bkb;
   ATL_cpnode_t *cb=NULL;

   mb = GetFlags(nargs, args, &pre, &outd, &bkb);
/*
 * Prep file for generation.  Free present values, and replace with:
 * ->auth  : kernel name without _b[1,n,0] suffix
 * ->genstr: for ID=0: genstr, else user kernel name (came in ->rout)
 * ->rout  : correct present filename
 */
   PrepMMForGen(pre, outd, "amm", mb);
/*
 * Expect bkb to be real,cplx mu==nu as first two entries, then UM kerns
 * Signal UM kerns to gen by overloading MMF_RIGHT to mean UM even when MU=NU
 */
   if (bkb)
   {
      prepSyrkAmm(pre, bkb);
      PrepMMForGen(pre, outd, "amm", bkb);
      bkb = SortListByIval_G2L(bkb, GetOffset(bkb, &bkb->next),
                               GetOffset(bkb, &bkb->blask));
      cb = GetKernCopies(pre, outd, bkb);
   }
   GenAllFiles(pre, outd, mb, bkb, cb);

   free(outd);
   KillAllCPNodes(cb);
   KillAllMMNodes(mb);
   return(0);
}
