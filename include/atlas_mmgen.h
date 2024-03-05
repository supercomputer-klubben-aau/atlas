#ifndef ATLAS_MMGEN_H
   #define ATLAS_MMGEN_H
/*
 * This file contains helper functions for code generation
 */
#define ATL_GENERATE 1
#include "atlas_cpparse.h"
#include "atlas_mmparse.h"
#include "atlas_mmtesttime.h"

#define CPV_C2BLK  0
#define CPV_BLK2C  1
#define CPV_A2BLK  2
#define CPV_BLK2A  3
#define CPV_BE1C   4
#define CPV_BENC   5
#define CPV_BEXC   6
#define CPV_BE0C   7
#define CPV_AL1C   8
#define CPV_ALNC   9
#define CPV_ALXC  10
#define CPV_AL1A  11
#define CPV_ALNA  12
#define CPV_ALXA  13
#define CPV_NOACP 14
#define CPV_NOCCP 15
#define CPV_NOMMI 16
#define CPV_NOUNR 17
#define CPV_SELFK 18  /* only called with kbB, self-clean if not mul ku */

#define CPV_ALLDIRA ( (1<<CPV_BLK2A)|(1<<CPV_A2BLK) )
#define CPV_ALLDIRC ( (1<<CPV_BLK2C)|(1<<CPV_C2BLK) )
#define CPV_ALLALA ( (1<<CPV_AL1A)|(1<<CPV_ALNA)|(1<<CPV_ALXA) )
#define CPV_ALLALC ( (1<<CPV_AL1C)|(1<<CPV_ALNC)|(1<<CPV_ALXC) )
#define CPV_ALLBEC ( (1<<CPV_BE1C)|(1<<CPV_BENC)|(1<<CPV_BEXC)|(1<<CPV_BE0C) )
#define CPV_ALLSUP ((1<<CPV_NOACP)|(1<<CPV_NOCCP)|(1<<CPV_NOMMI)|(1<<CPV_NOUNR))

typedef struct ViewNode ATL_view_t;
struct ViewNode
{
   char *fnam;
   char *nam;
   int flag;
   ATL_view_t *next;
};

static ATL_view_t *ATL_NewView(int flag, char *nam, char *fnam)
{
   ATL_view_t *p;
   p = malloc(sizeof(ATL_view_t));
   assert(p);
   p->flag = flag;
   p->nam  = nam;
   p->fnam = fnam;
   p->next = NULL;
   return(p);
}

static ATL_view_t *KillView(ATL_view_t *p)
{
   ATL_view_t *next=NULL;
   if (p)
   {
      next = p->next;
      if (p->fnam)
         free(p->fnam);
      if (p->nam)
         free(p->nam);
      free(p);
   }
   return(next);
}

static void KillAllViews(ATL_view_t *p)
{
   while (p)
      p = KillView(p);
}

static char *View2Args(ATL_view_t *p)
{
   const unsigned int flag=p->flag;
   unsigned int i, k;
   char *Ca, *Aa, *Ad, *Cd, *sp, *SUP;
   char scl[4]={'1', 'N', 'X', '0'};
   char Cb[12]={'C', 'b', '=', '\0'};
   char ks[4]={'K','=','1', '\0'};

   ks[2] = (flag&(1<<CPV_SELFK)) ? '0' : '1';
   i = (flag>>CPV_C2BLK)&3;
   if (i == 3)
      Cd ="Cd=I,F";
   else if (i == 0)
      Cd = "";
   else
      Cd = (i == 1) ? "Cd=F" : "Cd=I";

   i = (flag>>CPV_A2BLK)&3;
   if (i == 3)
      Ad ="Ad=I,F";
   else if (i == 0)
      Ad = "";
   else
      Ad = (i == 1) ? "Ad=F" : "Ad=I";

   if (Cd[0] != '\0')
   {
      i = (flag>>CPV_BE1C)&0xF;
      if (i == 0xF || !i)
         strcat(Cb, "0,1,N,X");
      else if (i)
      {
         for (k=0; (i&1) == 0 && k <= 4; k++, i>>=1);
         Cb[3] = scl[k];
         sp = Cb+4;
         for (i>>=1; i; k++, i>>=1)
         {
            if (i&1)
            {
               *sp = ',';
               sp[1] = scl[k];
               sp += 2;
            }
         }
         *sp = '\0';
      }
   }
   else
      Cb[0] = '\0';
   if (Ad[0] != '\0')
   {
      i = (flag>>CPV_AL1A)&7;
      switch(i)
      {
      case 0:
         Aa = "";
         break;
      case 1:
         Aa = "Aa=1";
         break;
      case 2:
         Aa = "Aa=N";
         break;
      case 3:
         Aa = "Aa=1,N";
         break;
      case 4:
         Aa = "Aa=X";
         break;
      case 5:
         Aa = "Aa=1,X";
         break;
      case 6:
         Aa = "Aa=N,X";
         break;
      default:
         Aa = "Aa=1,N,X";
      }
   }
   else
      Aa = "";
   if (Cd[0] != '\0')
   {
      i = (flag>>CPV_AL1C)&7;
      switch(i)
      {
      case 0:
         Ca = "";
         break;
      case 1:
         Ca = "Ca=1";
         break;
      case 2:
         Ca = "Ca=N";
         break;
      case 3:
         Ca = "Ca=1,N";
         break;
      case 4:
         Ca = "Ca=X";
         break;
      case 5:
         Ca = "Ca=1,X";
         break;
      case 6:
         Ca = "Ca=N,X";
         break;
      default:
         Ca = "Ca=1,N,X";
      }
   }
   else
      Ca = "";
   i = (flag&CPV_ALLSUP)>>CPV_NOACP;
   switch(i)
   {
   case 0x0: /* 0b0000 */
      SUP = "";
      break;
   case 0x1: /* 0b0001 */
      SUP = " S=A";
      break;
   case 0x2: /* 0b0010 */
      SUP = " S=C";
      break;
   case 0x3: /* 0b0011 */
      SUP = " S=C,A";
      break;
   case 0x4: /* 0b0100 */
      SUP = " S=M";
      break;
   case 0x5: /* 0b0101 */
      SUP = " S=M,A";
      break;
   case 0x6: /* 0b0110 */
      SUP = " S=M,C";
      break;
   case 0x7: /* 0b0111 */
      SUP = " S=M,C,A";
      break;
   case 0x8: /* 0b1000 */
      SUP = " S=U";
      break;
   case 0x9: /* 0b1001 */
      SUP = " S=U,A";
      break;
   case 0xA: /* 0b1010 */
      SUP = " S=U,C";
      break;
   case 0xB: /* 0b1011 */
      SUP = " S=U,C,A";
      break;
   case 0xC: /* 0b1100 */
      SUP = " S=U,M";
      break;
   case 0xD: /* 0b1101 */
      SUP = " S=U,M,A";
      break;
   case 0xE: /* 0b1110 */
      SUP = " S=U,M,C";
      break;
   case 0xF: /* 0b1111 */
      SUP = " S=U,C,M,U";
      break;
   default:
      assert(i == 0);
   }

   i = strlen(p->nam) + strlen(p->fnam) + strlen(Cb) + 26 + 2*5 + 2*4
     + strlen(SUP) + strlen(ks);
   sp = malloc(i);
   assert(sp);
   k = sprintf(sp, "-V %s %s %s %s %s %s %s %s %s", Ca, Aa, Cb, Ad, Cd,
               ks, SUP, p->nam, p->fnam);
   assert(k < i);
   return(sp);
}

/* procedure 1 */
int CPV_ScalStr2bits(char *st, int ibit)
{
   int flag=0;
   GET_SCAL:
   {
      switch(*st)
      {
      case '1':
         flag |= (1L<<ibit);
         break;
      case 'N':
         flag |= (1L<<(ibit+1));
         break;
      case 'X':
         flag |= (1L<<(ibit+2));
         break;
      case '0':
         flag |= (1L<<(ibit+3));
         break;
      default:
         assert(0);
      }
   }
   if (*(++st) == ',')
   {
      st++;
      goto GET_SCAL;
   }
   return(flag);
}

/* procedure 2 */
int CPV_DirStr2bits(char *st, int ibit)
{
   int flag;
   char ch = *st++;
   flag = (ch == 'F' || ch == 'f') ? (1<<ibit) : (1<<(ibit+1));
   if (*st++ == ',')
   {
      char ch = *st;
      flag |= (ch == 'F' || ch == 'f') ? (1<<ibit) : (1<<(ibit+1));
   }
   return(flag);
}

/* procedure 3 */
void PrintBegIfdef(FILE *fpout, char *nm)
{
   int i, k;
   assert(nm && fpout);
   for (k=0; k < 2; k++)
   {
      if (!k)
         fprintf(fpout, "#ifndef ");
      else
         fprintf(fpout, "   #define ");
      for (i=0; nm[i]; i++)
      {
         char ch = nm[i];
         ch = toupper(ch);
         ch = (ch == '.') ? '_' : ch;
         fputc(ch, fpout);
      }
      if (!k)
         fprintf(fpout, "\n");
      else
         fprintf(fpout, " 1\n");
   }
}

/* procedure 4 */
ATL_cpnode_t *GetCopyNodeFromMM(int flag, ATL_mmnode_t *mp, ATL_mmnode_t *ML)
/*
 * ML NULL means rout only called with K=kbB, so self-clean if kbB%ku != 0. Else
 * ML is master list of mm kernels, where ivar is set to the entry #+1 of the
 * K-clean node (ivar=0 means self-cleaning or no K-clean needed).
 */
{
   ATL_cpnode_t *cp;

   cp = GetCPNode();
   cp->flag = flag;
   if (FLAG_IS_SET(mp->flag, MMF_KVEC))
      cp->kvec = (mp->vlen > 1) ? mp->vlen : 0;
   if (mp->blask == ATL_KSYRK && (flag&(1<<CPF_CBLK)))
      cp->flag |= (1<<CPF_SYRK);
   else if (mp->blask == ATL_KSYMM)
      cp->flag |= (1<<CPF_SYMM);
   if (flag&(1<<CPF_CBLK))
   {
      cp->STGID = mp->stgC;
      cp->mb = mp->mbB;
      cp->nb = mp->nbB;
      cp->mu = mp->mu;
      cp->nu = mp->nu;
   }
   else /* A or B */
   {
      if (flag&(1<<CPF_ABLK))  /* A */
      {
         cp->STGID = mp->stgA ;
         cp->mb = mp->mbB;
         cp->nb = mp->kbB;
         cp->nu = mp->mu;
      }
      else                      /* B */
      {
         cp->STGID = mp->stgB;
         cp->mb = mp->kbB;
         cp->nb = mp->nbB;
         cp->nu = mp->nu;
      }
      cp->mu = 1;
      if (FLAG_IS_SET(mp->flag, MMF_KVEC))  /* KVEC only cleans to vlen */
         cp->mu = mp->vlen;
      else if (mp->blask) /* non-gemm kernels always pad to ku */
         cp->mu = mp->ku; /* so can handle all K in face of unrolling */
      else if (ML)/* to find copy padding, must examine cleaning node in ML */
      {
         ATL_mmnode_t *ml;
         ml = MMKernCompIsPresent(ML, mp);
         assert(ml);
         if (ml->ivar)
         {
            const int n=ml->ivar-1;
            int i;
            for(i=0, ml=ML; i < n && ml; i++, ml = ml->next);
            assert(ml);
            assert(MMKernsCompat(mp, ml));
         }
         else ml = mp;
         cp->mu = ml->ku;
      }
      else if (mp->kbB % mp->ku)  /* self-cleaning with non-mult kbB */
         cp->mu = mp->ku;
      else /* no cleanup needed because kbB is multiple of ku */
         cp->mu = 1;
   }
   return(cp);
}

/* procedure 5 */
ATL_cpnode_t *GetAllCopyNodesFromMM(int flag, ATL_mmnode_t *mb,
                                    ATL_mmnode_t *ML)
{
   ATL_cpnode_t *cb=NULL, *cp;
   ATL_mmnode_t *mp;

   for (mp=mb; mp; mp = mp->next)
   {
      cp = GetCopyNodeFromMM(flag, mp, ML);
      cp->next = cb;
      cb = cp;
   }
   return(ReverseCPQ(cb));
}

/* procedure 6 */
int *MMGetCopyIdxsFromList(int flag, ATL_mmnode_t *mb, ATL_cpnode_t *cb,
                           ATL_mmnode_t *ML)
/*
 * RETURNS: integer array of indices (starting from 0) in cb, 1st elt is len
 */
{
   ATL_mmnode_t *mp;
   int i, N, *idxs;
   if (!mb)
      return(NULL);
   N = ATL_CountNumberOfMMNodes(mb);
   idxs = malloc((N+1)*sizeof(int));
   assert(idxs);
   idxs[0] = N;
   for (i=0,mp=mb; mp; i++, mp = mp->next)
   {
      ATL_cpnode_t *cm, *cp;
      int id;

      cm = GetCopyNodeFromMM(flag, mp, ML);
      assert(cm);
      for (id=0,cp=cb; cp; cp = cp->next, id++)
         if (!CopyAreDiff(cp, cm))
            break;
      if (!cp)
      {
         fprintf(stderr, "CANNOT FIND COPY FOR:\n");
         PrintCPLine(stderr, cm);
         PrintMMLine(stderr, mp);
         fprintf(stderr, "\nFAILING LIST: %p\n", cb);
         PrintCPNodes(stderr, cb);
         assert(cp);
      }
      KillCPNode(cm);
      idxs[i+1] = id;
   }
   return(idxs);
}

/* procedure 7 */
ATL_cpnode_t *MMGetCopiesFromList(int flag, ATL_mmnode_t *mb, ATL_cpnode_t *cb,
                                  ATL_mmnode_t *ML)
{
   ATL_cpnode_t *mcb=NULL;
   ATL_mmnode_t *mp;
   for (mp=mb; mp; mp = mp->next)
   {
      ATL_cpnode_t *cm, *cp;
      cm = GetCopyNodeFromMM(flag, mp, ML);
      cp = FindEquivCopy(cb, cm);
      if (!cp)
      {
         fprintf(stderr, "NO COPY: ID=%d '%s', mu=%d, nu=%d, flg=%x\n",
                 cm->ID, cm->rout?cm->rout:"gen", cm->mu, cm->nu, cm->flag);
         PrintCPLine(stderr, cm);
      }
      assert(cp);
      KillCPNode(cm);
      cm = CloneCPNode(cp);
      cm->next = mcb;
      mcb = cm;
   }
   return(ReverseCPQ(mcb));
}


/* procedure 8 */
static int MMCopyGetKernID(int flag)
{
   int kernID;
   if (flag&(1<<CPF_TOBLK)) /* col-maj to block */
   {
      if (flag&(1<<CPF_CBLK)) /* copy C matrix */
      {
         if (flag&(1<<CPF_SYRK))
            kernID = ATL_KSKCPFC;
         else
            kernID = ATL_KGECPFC;
      }
      else /* involves A or B matrix */
         kernID = ATL_KGECPFA;
   }
   else  /* block to col-maj */
   {
      if (flag&(1<<CPF_CBLK)) /* copy C matrix */
      {
         if (flag&(1<<CPF_SYRK))
            kernID = ATL_KSKCP2C;
         else
            kernID = ATL_KGECP2C;
      }
      else /* involves A or B matrix */
         kernID = ATL_KGECP2A;
   }
   return(kernID);
}

/* procedure 9 */
void MMCopyTimePrep1(char pre, int kernID, int ialp, int ibet, ATL_mmnode_t *mp)
{
   mp->flag = (mp->flag & (~MMF_MVSET))|(1<<MMF_MVC);
   mp->blask = kernID;
   mp->ID = 0;
   mp->TB = mp->TA = AtlasTrans; /* generally, worst case perf-wise */
   if (mp->rout)
     free(mp->rout);
   mp->rout = DupString("ATL_tmp.c");
   if (mp->comp)
      free(mp->comp);
   if (mp->cflags)
      free(mp->cflags);
   if (mp->auth)
      free(mp->auth);
   mp->auth = mp->comp = mp->cflags = NULL;
   if (mp->genstr)
     free(mp->genstr);
   mp->genstr = MMGetCpGenString(pre, mp, ialp, ibet);
   mp->ivar = 1;  /* effective VLEN for copy kern (not amm) */
}

/* procedure 10 */
void MMCopyTimePrep(char pre, int flag, ATL_mmnode_t *mb)
{
   ATL_mmnode_t *mp;
   int ialp, ibet=0;
   int kernID;

   kernID = MMCopyGetKernID(flag);
   ialp = CopyGetAlphaI(flag);
   if (flag&(1<<CPF_CBLK))
      ibet = CopyGetBetaI(flag);
   for (mp=mb; mp; mp = mp->next)
     MMCopyTimePrep1(pre, kernID, ialp, ibet, mp);
}


/* procedure 11 */
static ATL_cpnode_t *GetMMCopyNode(int ID, int mu, int nu, int kvec, int flag)
{
   ATL_cpnode_t *cp;
   cp = calloc(1, sizeof(ATL_cpnode_t));
   cp->ID = ID;
   cp->mu = mu;
   cp->nu = nu;
   cp->kvec = kvec > 1 ? kvec : 0;
   cp->flag = flag;
   assert(cp);
   return(cp);
}

/* procedure 12 */
static ATL_cpnode_t *CloneCopyNode(ATL_cpnode_t *b)
{
   ATL_cpnode_t *p;
   p = malloc(sizeof(ATL_cpnode_t));
   memcpy(p, b, sizeof(ATL_cpnode_t));
   p->next = NULL;
/*
 * Now get our own copies of strings
 */
   if (p->genstr)
      p->genstr = DupString(p->genstr);
   if (p->rout)
      p->rout = DupString(p->rout);
   return(p);
}

/* procedure 13 */
static ATL_cpnode_t *KillCopyNode(ATL_cpnode_t *cp)
{
   ATL_cpnode_t *next=NULL;
   if (cp)
   {
      if (cp->genstr)
         free(cp->genstr);
      if (cp->rout)
         free(cp->rout);
      next = cp->next;
      free(cp);
   }
   return(next);
}

/* procedure 14 */
static ATL_cpnode_t *KillAllCopyNodes(ATL_cpnode_t *cp)
{
   while (cp)
      cp = KillCopyNode(cp);
   return(NULL);
}

/* procedure 15: RETURNS ptr to entry in b if n already there, else NULL */
ATL_cpnode_t *FindCopy(ATL_cpnode_t *b, ATL_cpnode_t *n)
{
   ATL_cpnode_t *p;
   for (p=b; p; p = p->next)
      if (!CopyAreDiff(p, n))
         return(p);
   return(NULL);
}

/* procedure 16 */
ATL_cpnode_t *FindLastCopyNode(ATL_cpnode_t *b)
{
   if (b)
   {
      while(b->next)
         b = b->next;
   }
   return(b);
}

/* procedure 17 */
static ATL_cpnode_t *AddUniqueCopyNode(ATL_cpnode_t *cb, ATL_cpnode_t *cnb)
/*
 * Adds all copy entries in cnb to cb, if cb doesn't already have a
 * functionally equivalent entry. cnb is unchanged.
 */
{
   ATL_cpnode_t *p;
   for (p=cnb; p; p = p->next)
   {
      if (!FindCopy(cb, p))
      {
         ATL_cpnode_t *np;
         np = CloneCopyNode(p);
         np->next = cb;
         cb = np;
      }
   }
   return(cb);
}

/* procedure 18 */
static void PrepMMForGen(char pre, char *outd, char *nm, ATL_mmnode_t *mb)
/*
 * Prep mb for generation.  Free present values, and replace with:
 * ->auth  : kernel name without _b[1,n,0] suffix
 * ->genstr: for ID=0: genstr, else user kernel name (came in ->rout)
 * ->rout  : correct present filename (used in generation & compilation)
 */
{
   ATL_mmnode_t *mp;
   char suff[3]={'.', 'c', '\0'};
   char *od;

   od = NewMergedString(outd, "/");
   for (mp=mb; mp; mp = mp->next)
   {
      char cpr = pre;
      if (mp->flag & (1<<MMF_COMPLEX))
         cpr = (pre == 's') ? 'c' : 'z';
      if (mp->flag & (1<<MMF_KUISKB))
         assert((mp->flag & (1<<MMF_KRUNTIME)) == 0);
      if (mp->auth)
         free(mp->auth);
      mp->auth = GetMMKernName(cpr, nm, mp);
      if (mp->ID) /* user-supplied kernel */
      {
         int h;
         char ch = 'c', och;

         assert(!mp->genstr);
         assert(mp->rout);
         h = strlen(mp->rout);
         mp->genstr = mp->rout;
         mp->rout = GetMMFilename(pre, nm, mp);
      }
      else  /* generated kernel */
      {
         char *sp, *sp1;
         if (mp->rout)
            free(mp->rout);
         if (mp->genstr)
            free(mp->genstr);
         sp = GetMMFilename(cpr, nm, mp);
         mp->rout = NewMergedString(od, sp);  /* want to gen in outd/ */
         sp1 = MMGetGenString(cpr, mp);
         mp->genstr = NewMergedString(sp1, " ");
         free(sp1);
         free(mp->rout);
         mp->rout = sp; /* want filename w/o path for GenMake */
      }
   }
   free(od);
}

/* procedure 19 */
static FILE *OpenMMGenHeader
(
   char *outd,/* directory to write to must be non-NULL */
   int UID,   /* 0 for ATLAS files, >0 for uammsearch */
   char pre,  /* [d,s,z,c] */
   char *ip,  /* "[ip,op] */
   char *cn,  /* "[ge,sq,mn,dm,dn,nk,mk]" */
   char *nm,  /* "[blk,kern,flag,perf,sum,]" */
   ATL_mmnode_t *mb
)
{
   char *fn, bs[6];
   int len, i;
   FILE *fp;

   /* NAME WILL BE: <outd>/atlas_<pre><ip><cn>[ID]_<nm>.h */
   assert(nm && outd);
   len = strlen(outd) + 1; /* outd+/ */
   len += 6 + 1 + 1 + 2;  /* atlas+pre+_+.h */
   if (ip)
      len += strlen(ip);
   else
      ip = "";
   if (cn)
      len += strlen(cn);
   else
      cn = "";
   if (UID)
      len += NumDecDigits(UID);
   len += strlen(nm);
   bs[0] = 0;

   fn = malloc(len+1);
   assert(fn);
   if (UID)
      i = sprintf(fn, "%s/atlas_%c%s%s%u_%s.h", outd, pre, ip, cn, UID, nm);
   else
      i = sprintf(fn, "%s/atlas_%c%s%s_%s.h", outd, pre, ip, cn, nm);
   assert(i <= len);
   fp = fopen(fn, "w");
   assert(fp);
/*
 * Now change file part of string to guard text: uppercase & ".h" -> "_H";
 */
   fn[i-1] = 'H';
   fn[i-2] = '_';
   for (i -= 3; fn[i] != '/'; i--)
      fn[i] = toupper(fn[i]);
   i++;
   fprintf(fp, "#ifndef %s\n", fn+i);
   fprintf(fp, "   #define %s 1\n", fn+i);
   fprintf(fp, "#include \"atlas_amm.h\"\n");
   if (mb && strcmp(nm, "perf"))
   {
      i = CountListEntries(mb, GetOffset(&mb->next, mb));
      if (cn[0] == 'g' && cn[1] == 'e')
      {
         fprintf(fp, "#ifdef ATL_AMM_NCASES\n   #if ATL_AMM_NCASES != %d\n", i);
         fprintf(fp, "      #error \"NCASES MISMATCH!\"\n   #endif\n");
         fprintf(fp, "#else\n   #define ATL_AMM_NCASES %d\n#endif\n", i);
      }
      else
      {
         fprintf(fp, "#ifdef ATL_%sAMM_NCASES\n   #if ATL_%sAMM_NCASES != %d\n",
                 cn, cn, i);
         fprintf(fp, "      #error \"NCASES MISMATCH!\"\n   #endif\n");
         fprintf(fp, "#else\n   #define ATL_%sAMM_NCASES %d\n#endif\n", cn, i);
      }
   }

   free(fn);
   return(fp);
}

/* procedure 20 */
static void CloseGenHeader(FILE *fp)
{
   fprintf(fp, "\n#endif\n");
   fclose(fp);
}

/* procedure 21 */
static char *GetSafeIntType(void *vb, int nxtoff, int off, int *CNT)
{
   int max, min, n;
   *CNT = GetIntMaxMinAtOff(vb, nxtoff, off, &max, &min);
   if (min < 0)
   {
      if (max < (1<<7))
         return("char");
      else if (max < (1<<15))
         return("short");
      return("int");
   }
   else if (max < (1<<8))
      return("unsigned char");
   else if (max < (1<<16))
      return("unsigned short");
   return("unsigned int");
}


char *FillInMMScalarSuff(char ta, char alp, char bet)
{
   static char ab[8];
   if (!(alp|bet))
      ab[0] = '\0';
   else
   {
      int i;
      i = 1;
      ab[0] = '_';
      if (!bet)
      {
         ab[1] = ta;
         ab[2] = 'a';
         ab[3] = alp;
         ab[4] = '\0';
      }
      else if (!alp)
      {
         ab[1] = 'b';
         ab[2] = bet;
         ab[3] = '\0';
      }
      else
      {
         ab[1] = 'a';
         ab[2] = alp;
         ab[3] = '_';
         ab[4] = 'b';
         ab[5] = bet;
         ab[6] = '\0';
      }
   }
   return(ab);
}

/* procedure 22 */
void PrintMMProtos(FILE *fp, char pre, char *nm, ATL_mmnode_t *mb, int off,
                   char bet)
{
   ATL_mmnode_t *mp;

   if (bet == '0')
      fprintf(fp, "#if !defined(NO%s_b%c) && !defined(NO%s_K1)\n", nm,bet, nm);
   else
      fprintf(fp, "#if !defined(NO%s_b%c) && !defined(NO%s_K1_b%c)\n",
              nm,bet, nm,bet);
   for (mp=mb; mp; mp = mp->next)
      fprintf(fp, "void %s_b%c\n%s\n%s\n", GetStrAtOff(mp, off), bet,
              "   (ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,TYPE*,",
              "    const TYPE*,const TYPE*,const TYPE*);");
   fprintf(fp, "#endif\n\n");
}

/* procedure 23 */
void PrintMMCpProtosA(FILE *fp, char *nm, ATL_cpnode_t *cb, char dir, char alp)
{
   ATL_cpnode_t *cp;
   char *pro;
   char pre;

   pre = CopyGetPre(cb->flag);
   if (pre == 'd' || pre == 's')
   {
      pro = (dir == 'F') ?
           "(ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*)"
         : "(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*)";
   }
   else
   {
      pro = (dir == 'F') ?
      "(ATL_CSZT,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT,const TYPE*,const TYPE*);"
      : "(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,TYPE*,TYPE*);";
   }
   fprintf(fp, "#ifndef NO%s_a%c\n", nm, alp);
   for (cp=cb; cp; cp = cp->next)
   {
      char ta;
      ta = CopyGetTrans(cp->flag);
      fprintf(fp, "void %s_%ca%c\n   %s;\n", cp->rout, ta, alp, pro);
   }
   fprintf(fp, "#endif\n\n");
}

/* procedure 24 */
void PrintMMCpProtosC(FILE *fp, char *nm, ATL_cpnode_t *cb,
                      char dir, char alp, char bet)
{
   ATL_cpnode_t *cp;
   char *ab;
   char *pro;
   char pre;

   pre = CopyGetPre(cb->flag);
   ab = FillInMMScalarSuff(0, alp, bet);
   if (pre == 'd' || pre == 's')
   {
      pro = (dir == 'F') ?
"(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,TYPE*,ATL_CSZT);"
:
"(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,TYPE*);";
   }
   else
   {
      pro = (dir == 'F') ?
"(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,const SCALAR,\n     TYPE *,ATL_CSZT);"
            :
"   (ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,ATL_CSZT,const SCALAR,\n     TYPE*,TYPE*);";
   }
   fprintf(fp, "#ifndef NO%s%s\n", nm,ab);
   for (cp=cb; cp; cp = cp->next)
      fprintf(fp, "void %s%s\n   %s\n", cp->rout, ab, pro);
   fprintf(fp, "#endif\n\n");
}

/* procedure 25 */
void PrintIntArrAtOff
(
   FILE *fp,
   char pre,   /* [d,s,c,z] */
   char *nm,   /* name of the array, optionally decorated by _aX_bX */
   void *sb,
   int nxtoff, /* offset (bytes) in struct to next pointer */
   int ioff,   /* offset in struct to integer val we are printing */
   char alp,
   char bet
)
/*
 * Dumps integer stored at ioff in sb's struct to file fp using name nm.

 */
{
   char *cp = sb;
   char *type;
   char *ab;
   int i, n;

   ab = FillInMMScalarSuff(0, alp, bet);

   type = GetSafeIntType(sb, nxtoff, ioff, &n);
   fprintf(fp, "#ifndef NO%s\n", nm);
   fprintf(fp, "static const %s ATL_AMM_%s%s[%d] =\n{\n", type, nm, ab, n);
   for (i=0; cp; i++)
   {
      char *nxt;
      nxt = GetStrAtOff(cp, nxtoff);
      fprintf(fp, "%8d%c  /* CNT=%d */\n", *((int*)(cp+ioff)), nxt?',':' ', i);
      cp = nxt;
   }
   fprintf(fp, "};\n#endif\n\n");
}
/* procedure 26 */
void PrintStrArrAtOff(FILE *fp, char pre, char *nm, void *vb,
                      int nxtoff, int off, char *type,
                      char ta, char alp, char bet)
/*
 * Dumps integer stored at byte offset off mb's MMNodes to file fp using name nm
 */
{
   char *cp;
   int i, n;
   char *ab;

   n = CountListEntries(vb, nxtoff);
   ab = FillInMMScalarSuff(ta, alp, bet);
   fprintf(fp, "#ifndef NO%s%s\n", nm,ab);
   if (!bet) /* A/B copy funcs */
      fprintf(fp, "static const %s ATL_AMM_%s_a%c[%d] =\n{\n",type,nm,alp,n);
   else
      fprintf(fp, "static const %s ATL_AMM_%s%s[%d] =\n{\n", type, nm, ab, n);
   for (i=0,cp=vb; cp; i++)
   {
      char *next;
      next = GetStrAtOff(cp, nxtoff);
      fprintf(fp, "/* IDX=%3d */ %s%s%c\n",
              i, GetStrAtOff(cp, off), ab, next ? ',':' ');
      cp = next;
   }
   fprintf(fp, "};\n#endif\n\n");
}

/* procedure 27 */
char *GetMMKernComp(ATL_mmnode_t *mmp, char *dcomp, char *dflags, char **flgs)
{
   char *comp = dcomp;
   if (mmp->comp)
   {
      comp = (mmp->comp[0] == 'g' && mmp->comp[1] == 'c' &&
              mmp->comp[2] == 'c' &&
             (mmp->comp[3] == '\0' || mmp->comp[3] == ' '))
             ? "$(GOODGCC)" : mmp->comp;
      *flgs = mmp->cflags;
   }
   else
      *flgs = dflags;
   return(comp);
}

void PrintMakeTargs(FILE *fp, char pre)
{
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
}
#endif  /* end guard around atlas_mmgen.h */
