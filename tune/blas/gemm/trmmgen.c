/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 Majedul Sujon
 */
#include "atlas_mmgen.h"
#include "atlas_type.h"

void GenTrCopyHeader(FILE *fp, char pre, char sd, char up0, char ta0, char up,
      char ta, ATL_mmnode_t *mp)
{
/*
 * Acopy always TRMM copy
 * LN: TT2blk-Lo
 * UT: TN2blk-Up
 * two additional cases: unit-diagonal, conjugate
 *
 * NOTE:
 *
 * format:
 *    ATL_da2blk_trmm_LLN_Lo_a1               ATL_dcpFromATt_Lo_a1_1x24_0
 */
/*
 * selecting correct copies from table:
 *    S  U  T  kernels   cpyA    cpB
 *    ---------------------------------
 *    L  L  N    TRMM=1   tAT,L   AN
 *    L  L  T    TRMM=2   tAN,L   AN
 *    L  U  N    TRMM=2   tAT,U   AN
 *    L  U  T    TRMM=1   tAN,U   AN
 *    R  L  N    TRMM=3    AT    tAN,L
 *    R  L  T    TRMM=4    AT    tAT,L
 *    R  U  N    TRMM=4    AT    tAN,U
 *    R  U  T    TRMM=3    AT    tAT,U
 *
 */
   int i, j;
   char *uplo;
   char *uld; /* _DiaU */
   const char al[3] = {'1', 'N', 'X'};
   unsigned int UR, kvec;
   const char uls[2] = {'L', 'U'};
   const char tas[2] = {'N', 'T'};
   char at, tta;

   kvec = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
   uplo = (up == 'L')? "Lo" : "Up";
   if (sd == 'L')
   {
      at = (ta == 'N')? 'T' : 'N';
      UR = mp->mu;
   }
   else  /* R */
   {
      at = ta;
      UR = mp->nu;
   }
   for (j=0 ; j < 2; j++)
   {
      if (j == 1)
      {
         uld = "diagU_";
         fprintf(fp, "/*\n* Unit Diagonal \n*/\n");
      }
      else
         uld = "";

      for (i=0; i < 3; i++)
      {
/*
 *       copy name format
 *       ATL_<pre>cp[Into,From]A[N,T][g,k,y,r,s]_aX_<ku>x<nu>_<kvec>
 */
         if (pre == 'c' || pre == 'z')
         {
            tta = (at == 'N')? 'C' : 'H';
            fprintf(fp, "#ifdef Conj_\n");
            fprintf(fp,
                    "#define ATL_%ca2blk_trmm_%c%c%c_%s_%sa%c "
                    "ATL_%ccpFromA%ct_%s_%sa%c_%ux%u_%u\n",
                    pre, sd, up0, ta0, uplo, uld, al[i],
                    pre, tta, uplo, uld, al[i], mp->ku, UR, kvec);
            fprintf(fp, "#else\n");
         }
         fprintf(fp,
                 "#define ATL_%ca2blk_trmm_%c%c%c_%s_%sa%c "
                 "ATL_%ccpFromA%ct_%s_%sa%c_%ux%u_%u\n",
                 pre, sd, up0, ta0, uplo, uld, al[i],
                 pre, at, uplo, uld, al[i], mp->ku, UR, kvec);
         if (pre == 'c' || pre == 'z')
            fprintf(fp, "#endif\n");

         fprintf(fp, "void ATL_%ca2blk_trmm_%c%c%c_%s_%sa%c", pre, sd, up0, ta0,
                 uplo, uld, al[i]);
         if (pre == 'd' || pre == 's')
            fprintf(fp, "(ATL_CSZT,const TYPE,const TYPE*,ATL_CSZT,TYPE*);\n");
         else
            fprintf(fp,
                 "(ATL_CSZT,const TYPE*,const TYPE*, ATL_CSZT,TYPE*,TYPE*);\n");

         fprintf(fp, "#ifndef ATL_%ctrmm_a2blk_%s_%sa%c\n   #define "
                 "ATL_%ctrmm_a2blk_%s_%sa%c "
                 "ATL_%ca2blk_trmm_%c%c%c_%s_%sa%c\n#endif\n",
                 pre, uplo, uld, al[i],  pre, uplo, uld, al[i],  pre, sd, up0,
                 ta0, uplo, uld, al[i]);
      }
   }

}

void GenTrmmHead(char pre, char *outd, char sd, char up, char ta,
      ATL_mmnode_t *mp)
{
   char *of;
   FILE *fp;
   int i, flg, k, L;
   char *uplo;
   const char *st[5] = {"VLEN", "KVEC", "MU", "NU", "KU"};
   char uld [6]; /* _DiaU */
   const char be[3] = {'0', '1', 'n'};
   const char al[3] = {'1', 'N', 'X'};
   unsigned int UR, kvec;
   const char uls[2] = {'L', 'U'};
   const char tas[2] = {'N', 'T'};
   char at;
   ATL_cpnode_t *cb;

   L = strlen(outd) + 24;
   of = malloc(L);
   assert(of);
/*
 * LN represents both LN and UT case
 */
   if (ta == 'N')
      k = sprintf(of, "%s/atlas_%cutrmm%c_LN.h", outd, pre, sd);
   else
      k = sprintf(of, "%s/atlas_%cutrmm%c_LT.h", outd, pre, sd);
   assert(k<L);
   fp = fopen(of, "w");
   assert(fp);
   free(of);
/*
 * inclusion guard
 */
   fprintf(fp, "#ifndef ATLAS_%cTRMM%c_%c%c_H\n",
           toupper(pre), toupper(sd), toupper(up), toupper(ta));
   fprintf(fp, "   #define ATLAS_%cTRMM%c_%c%c_H 1\n\n",
           toupper(pre), toupper(sd), toupper(up), toupper(ta));

   fprintf(fp, "#include \"atlas_amm.h\"\n");
/*
 * data
 */
   assert(mp);
   fprintf(fp, "#define ATL_TRMM_%c%c%c_KB %d\n", sd, up, ta, mp->kbB);
   fprintf(fp, "#define ATL_TRMMK_%c%c%c_VLEN %d\n", sd, up, ta, mp->vlen);
   kvec = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
   fprintf(fp, "#define ATL_TRMMK_%c%c%c_KVEC %u\n", sd, up, ta, kvec);
   fprintf(fp, "#define ATL_TRMMK_%c%c%c_MU %d\n", sd, up, ta, mp->mu);
   fprintf(fp, "#define ATL_TRMMK_%c%c%c_NU %d\n", sd, up, ta, mp->nu);
   fprintf(fp, "#define ATL_TRMMK_%c%c%c_KU %d\n", sd, up, ta, mp->ku);

   for (i=0; i<5; i++)
   {
      fprintf(fp, "#ifndef ATL_TRMMK_%s\n", st[i]);
      fprintf(fp, "   #define ATL_TRMMK_%s ATL_TRMMK_%c%c%c_%s\n", st[i], sd,
              up, ta, st[i]);
      fprintf(fp, "#endif\n");
   }
   fprintf(fp, "#ifndef ATL_TRMM_KB\n");
   fprintf(fp, "   #define ATL_TRMM_KB ATL_TRMM_%c%c%c_KB\n", sd, up, ta);
   fprintf(fp, "#endif\n");
/*
 * kernel info
 */
   fprintf(fp, "/*\n* TRMM KERNEL INFO\n*/\n");
   for (i=0; i < 3; i++)
   {
      fprintf(fp,
              "void ATL_%ctrmmK_%c%c%c_b%c("
              "ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*, \n",
              pre, sd, up, ta, be[i]);
      fprintf(fp,
              "                            TYPE*, "
              "const TYPE*, const TYPE*, const TYPE*);\n");
      fprintf(fp, "#ifndef ATL_%ctrmmK_b%c\n", pre, be[i]);
      fprintf(fp, "   #define ATL_%ctrmmK_b%c ATL_%ctrmmK_%c%c%c_b%c\n", pre,
              be[i], pre, sd, up, ta, be[i]);
      fprintf(fp, "#endif\n");
   }
/*
 * blk2C
 * NOTE: we need copy routines for all alpha
 */
   fprintf(fp, "/*\n* BLK2C COPY\n*/\n");
   flg = (1<<CPF_CBLK);
   flg |= (pre == 'd' || pre == 's') ? (1<<CPF_REAL):0;
   flg |= (pre == 'c' || pre == 's') ? (1<<CPF_SINGLE):0;
   cb = GetCopyNodeFromMM(flg, mp, NULL);
   assert(cb);
/*
 * blk2c copy is same as gemm... so we need gemm style blk2c copy
 * so, disable trmm flag for this specific copy routine
 */
   cb->flag &= ~(1<<CPF_TRMM);
   cb->flag |= (1<<CPF_BE0);
   for (i=0; i < 3; i++)
   {
      if (al[i] == '1')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_ALX));
         cb->flag |= (1<<CPF_AL1);
      }
      else if (al[i] == 'N')
      {
         cb->flag &= ~((1<<CPF_ALX)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALN);
      }
      else if (al[i] == 'X')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALX);
      }
      if(i) free(cb->rout);
      cb->rout = GetCopyName(cb, 0);
      fprintf(fp, "#define ATL_%cblk2c_trmm_%c%c%c_a%cb0 %s\n",
              pre, sd, up, ta, al[i], cb->rout);
      fprintf(fp, "void ATL_%cblk2c_trmm_%c%c%c_a%cb0", pre, sd, up, ta, al[i]);
      if (pre == 'd' || pre == 's')
         fprintf(fp,
                 "(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const SCALAR,"
                 "TYPE*, ATL_CSZT);\n");
      else
         fprintf(fp,
                 "(ATL_CSZT,ATL_CSZT,const SCALAR,const TYPE*,const TYPE*,"
                 "const SCALAR,TYPE*, ATL_CSZT);\n");
      fprintf(fp, "#ifndef ATL_%ctrmm_blk2c_a%cb0\n   #define "
              "ATL_%ctrmm_blk2c_a%cb0 ATL_%cblk2c_trmm_%c%c%c_a%cb0\n#endif\n",
              pre, al[i], pre, al[i], pre, sd, up, ta, al[i]);
   }
   KillAllCopyNodes(cb);
/*
 * NOTE: depending on the side, Acopy and Bcopy changes
 * side = left: Acopy is trmm copy and Bcopy is gemm copy
 * side = right:Acopy is gemm copy and Bcopy is trmm copy although trmm copy is
 * a2blk and gemm copy is b2blk in header file
 */
/*
 * Bcopy: always gemm copy
 */
   fprintf(fp, "/*\n* B2BLK COPY\n*/\n");
   flg = 0;
   if (sd == 'R')
   {
      flg |= (1<<CPF_ABLK);
      flg |= (1<<CPF_TRANS);
   }

   flg |= (1<<CPF_TOBLK); /* b2blk*/
   flg |= (pre == 'd' || pre == 's') ? (1<<CPF_REAL):0;
   flg |= (pre == 'c' || pre == 's') ? (1<<CPF_SINGLE):0;
   cb = GetCopyNodeFromMM(flg, mp, NULL);
   assert(cb);
   cb->flag &= ~(1<<CPF_TRMM);
   for (i=0; i < 3; i++)
   {
      if (al[i] == '1')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_ALX));
         cb->flag |= (1<<CPF_AL1);
      }
      else if (al[i] == 'N')
      {
         cb->flag &= ~((1<<CPF_ALX)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALN);
      }
      else if (al[i] == 'X')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALX);
      }
      if(i) free(cb->rout);
      cb->rout = GetCopyName(cb, 0);
      fprintf(fp, "#define ATL_%cb2blk_trmm_%c%c%c_a%c %s\n",
              pre, sd, up, ta, al[i], cb->rout);
      fprintf(fp, "void ATL_%cb2blk_trmm_%c%c%c_a%c", pre, sd, up, ta, al[i]);
         if (pre == 'd' || pre == 's')
            fprintf(fp,
            "(ATL_CSZT,ATL_CSZT,const TYPE,const TYPE*,ATL_CSZT,TYPE*);\n");
         else
            fprintf(fp, "(ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,"
                        "ATL_CSZT,TYPE*,TYPE*);\n");
         fprintf(fp, "#ifndef ATL_%ctrmm_b2blk_a%c\n   #define "
         "ATL_%ctrmm_b2blk_a%c ATL_%cb2blk_trmm_%c%c%c_a%c\n#endif\n",
            pre, al[i],  pre, al[i],  pre, sd, up, ta, al[i]);
   }
   KillAllCopyNodes(cb);
/*
 * In eaqch header file, we have two options.
 * LN/UT and LT/UN
 */
   fprintf(fp, "/*\n* A2BLK COPY\n*/\n");
   if (ta == 'N')
   {
      GenTrCopyHeader(fp, pre, sd, 'L', 'N', 'L', 'N', mp);
      GenTrCopyHeader(fp, pre, sd, 'L', 'N', 'U', 'T', mp);
   }
   else /* LT UN case */
   {
      GenTrCopyHeader(fp, pre, sd, 'L', 'T', 'L', 'T', mp);
      GenTrCopyHeader(fp, pre, sd, 'L', 'T', 'U', 'N', mp);
   }
   fprintf(fp, "#endif\n");
}


void GenMicroTrmm(char pre, char *outd, char sd, char up, char ta,
                           ATL_mmnode_t *mb)
{
   int i, L;
   int trmmid=0;
   int mu=mb->mu, nu=mb->nu;
   char *ln, *mk;
   char rt[256];

   L = strlen(outd) + 256;
   ln = malloc(L);
   assert(ln);
/*
 * TRMM=1 (left lower), 2 (left upper), 3 (right lower), 4 (right upper).
 * Transform:
 *    LLN : 1           RLN : 3
 *    LLT : 2           RLT : 4
 *    LUN : 2           RUN : 4
 *    LUT : 1           RUT:  3
 *
 *    transpose reverses the upper/lower
 */
   if (sd == 'R')
      trmmid = 2;
   if (ta == 'N')
      trmmid += (up=='L')? 1 : 2;
   else
      trmmid += (up=='L')? 2 : 1;

   i = sprintf(rt, "%s/ATL_%ctrmmK_%c%c%c_%ux%ux%u.c", outd, pre, sd, up, ta,
               mb->mu, mb->nu, mb->ku);
   assert(i<256);
   if (mb->rout) free(mb->rout);
   mb->rout = DupString(rt);

   mk = MMGetGenString(pre, mb);
   i = sprintf(ln, "%s TRMM=%d", mk, trmmid);
   /*printf("genstr = %s\n", ln);*/
   free(mk);
   assert(i < L);
   assert(!Sys2File(ln, NULL));
   free(ln);
}

void AddCopyDirective(FILE *fp, char *nam, char pre, char alc, unsigned int UR)
{
   int i, k;
   int sz;
   char *cpn, *ccpn;

   sz = strlen(nam);
   sz += 10; /* might be added: diagU_, [up,lo]_*/
   cpn = malloc(sz);
   assert(cpn);
/*
 * split the string at a1/aX/aN
 */
   assert(nam);
   k = strlen(nam);
   for (i=0; i < k; i++, nam++)
   {
      if (*nam == 'a')
      {
         assert(*(nam+1)!='\0'); /* assumption: must consists of a[1,N,X] */
         if ( *(nam+1) == '1' || *(nam+1) == 'N' || *(nam+1) == 'X')
            break;
      }
      cpn[i] = *nam;
   }
   cpn[i] = '\0';
/*
 * for conj: convert ATt to AHt, ANt to ACt
 */
   if (pre == 'z' || pre == 'c')
   {
      ccpn = malloc(sz);
      assert(ccpn);
      strcpy(ccpn, cpn);
      k = strlen(ccpn);
      for (i=0; i < k-2; i++)
      {
         if (ccpn[i] == 'A' && ccpn[i+1] == 'N' && ccpn[i+2] == 't')
            ccpn[i+1] = 'C';
         else if (ccpn[i] == 'A' && ccpn[i+1] == 'T' && ccpn[i+2] == 't')
            ccpn[i+1] = 'H';
      }
   }

   fprintf(fp, "#ifdef SCALAR\n  #undef SCALAR\n#endif\n");
   fprintf(fp, "#ifdef TYPE\n  #undef TYPE\n#endif\n");
   fprintf(fp, "#ifdef SREAL\n  #undef SREAL\n#endif\n");
   fprintf(fp, "#ifdef DREAL\n  #undef DREAL\n#endif\n");
   fprintf(fp, "#ifdef SCPLX\n  #undef SCPLX\n#endif\n");
   fprintf(fp, "#ifdef DCPLX\n  #undef DCPLX\n#endif\n");
   fprintf(fp, "#ifdef ALPHA1\n  #undef ALPHA1\n#endif\n");
   fprintf(fp, "#ifdef ALPHAN\n  #undef ALPHAN\n#endif\n");
   fprintf(fp, "#ifdef ALPHAN1\n  #undef ALPHAN1\n#endif\n");
   fprintf(fp, "#ifdef ALPHAX\n  #undef ALPHAX\n#endif\n");
   fprintf(fp, "#ifdef BETA1\n  #undef BETA1\n#endif\n");
   fprintf(fp, "#ifdef BETAN\n  #undef BETAN\n#endif\n");
   fprintf(fp, "#ifdef BETAN1\n  #undef BETAN1\n#endif\n");
   fprintf(fp, "#ifdef BETAX\n  #undef BETAX\n#endif\n");
/*
 * directives for type
 */
   if (pre == 's')
   {
      fprintf(fp, "#define SREAL 1\n");
      fprintf(fp, "#define TYPE float\n");
      fprintf(fp, "#define SCALAR TYPE\n");
   }
   else if (pre == 'd')
   {
      fprintf(fp, "#define DREAL 1\n");
      fprintf(fp, "#define TYPE double\n");
      fprintf(fp, "#define SCALAR TYPE\n");
   }
   else if (pre == 'c')
   {
      fprintf(fp, "#define SCPLX 1\n");
      fprintf(fp, "#define TYPE float\n");
      fprintf(fp, "#define SCALAR TYPE*\n");
   }
   else if (pre == 'z')
   {
      fprintf(fp, "#define DCPLX 1\n");
      fprintf(fp, "#define TYPE double\n");
      fprintf(fp, "#define SCALAR TYPE*\n");
   }
   else
   {
      printf("UNKNOWN TYPE\n");
      assert(0);
   }
/*
 * directives alpha
 */
   if (alc == '1')
      fprintf(fp, "#define ALPHA1 1\n");
   else if (alc == 'X')
      fprintf(fp, "#define ALPHAX 1\n");
   else if (alc == 'N')
   {
      fprintf(fp, "#define ALPHAN 1\n");
      fprintf(fp, "#define ALPHAN1 1\n");
   }
/*
 * UR &
 */
   fprintf(fp, "#define ATL_NU %u\n", UR);
   fprintf(fp, "\n#ifdef UnitDiag_\n");
   fprintf(fp, "   #ifdef Upper_\n");
   if (pre == 'z' || pre == 'c')
   {
      fprintf(fp, "      #ifdef Conj_\n");
      fprintf(fp, "         #define ATL_USERCPMM %sUp_diagU_%s\n", ccpn, nam);
      fprintf(fp, "      #else\n");
   }
   fprintf(fp, "         #define ATL_USERCPMM %sUp_diagU_%s\n", cpn, nam);
   if (pre == 'z' || pre == 'c') fprintf(fp, "      #endif\n");
   fprintf(fp, "   #else\n");
   if (pre == 'z' || pre == 'c')
   {
      fprintf(fp, "      #ifdef Conj_\n");
      fprintf(fp, "         #define ATL_USERCPMM %sLo_diagU_%s\n", ccpn, nam);
      fprintf(fp, "      #else\n");
   }
   fprintf(fp, "         #define ATL_USERCPMM %sLo_diagU_%s\n", cpn, nam);
   if (pre == 'z' || pre == 'c')
      fprintf(fp, "      #endif\n");
   fprintf(fp, "   #endif\n");
   fprintf(fp, "#else\n");
   fprintf(fp, "   #ifdef Upper_\n");
   if (pre == 'z' || pre == 'c')
   {
      fprintf(fp, "      #ifdef Conj_\n");
      fprintf(fp, "         #define ATL_USERCPMM %sUp_%s\n", ccpn, nam);
      fprintf(fp, "      #else\n");
   }
   fprintf(fp, "         #define ATL_USERCPMM %sUp_%s\n", cpn, nam);
   if (pre == 'z' || pre == 'c')
      fprintf(fp, "      #endif\n");
   fprintf(fp, "   #else\n");
   if (pre == 'z' || pre == 'c')
   {
      fprintf(fp, "      #ifdef Conj_\n");
      fprintf(fp, "         #define ATL_USERCPMM %sLo_%s\n", ccpn, nam);
      fprintf(fp, "      #else\n");
   }
   fprintf(fp, "         #define ATL_USERCPMM %sLo_%s\n", cpn, nam);
   if (pre == 'z' || pre == 'c')
      fprintf(fp, "      #endif\n");
   fprintf(fp, "   #endif\n");
   fprintf(fp, "#endif\n\n");

   if (pre == 'z' || pre == 'c')
      free(ccpn);
   free(cpn);
}

void GenTrCopy(ATL_mmnode_t *mp, ATL_cpnode_t *cb, char *outd, char pre,
      char sd, char up, char ta)
{
   int i, k, L;
   int UR, kvec;
   char *ln, *cpg;
   char rt[256];
   const char al[3] = {'1', 'N', 'X'};
   char tta;
   FILE *fp;

   L = strlen(outd) + 256;
   ln = malloc(L);
   assert(ln);

   cb->flag &= ~(1<<CPF_TRANS);
/*
 * for Left TRMM: we need tAT when ta=N and tAN when ta=T
 * since we generate the copyname by our own, tta needs to reflect that
 */
   #if 0
   if (sd == 'L')
      cb->flag |= (ta=='N') ? (1<<CPF_TRANS) : 0;
   else
      cb->flag |= (ta=='T') ? (1<<CPF_TRANS) : 0;
   #else
   if (sd == 'L')
      tta = (ta=='N')? 'T' : 'N';
   else
      tta = ta;
      cb->flag |= (tta=='T') ? (1<<CPF_TRANS) : 0;
   #endif
/*
 * NOTE: we don't need up/lo as the parameter for TRMM-copy anymore
 */
   /*if (up == 'L')
      cb->flag |= (1<<CPF_LOWER);
   else
      cb->flag &= ~(1<<CPF_LOWER);
   */
   UR = (sd=='L')? mp->mu : mp->nu;
   kvec = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;

   for (i=0; i < 3; i++)
   {
      if (al[i] == '1')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_ALX));
         cb->flag |= (1<<CPF_AL1);
      }
      else if (al[i] == 'N')
      {
         cb->flag &= ~((1<<CPF_ALX)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALN);
      }
      else if (al[i] == 'X')
      {
         cb->flag &= ~((1<<CPF_ALN)|(1<<CPF_AL1));
         cb->flag |= (1<<CPF_ALX);
      }
#if 0
      cpg = GetCopyName(cb, 0);
#else
      cpg = malloc(64*sizeof(char));
      assert(cpg);
      sprintf(cpg, "ATL_%ccpFromA%ct_a%c_%ux%u_%u", pre, tta, al[i], mp->ku,
             UR, kvec);
#endif
      k = sprintf(rt, "%s/%s.c", outd, cpg);
      assert(k<256);
/*
 *    add copy directives at the start of the source file
 */
      fp = fopen(rt, "w");
      assert(fp);
      AddCopyDirective(fp, cpg, pre, al[i], cb->nu);/*always nu ifnot CPF_CBLK*/
      fclose(fp);
      free(cpg);
      if (cb->rout) free(cb->rout);
      cb->rout = DupString("ATL_tmp.c");
      cpg = GetCopyGenStr(cb);
      k = sprintf(ln, "%s", cpg);
      free(cpg);
      assert(k < L);
      assert(!Sys2File(ln, NULL));
/*
 *    append files
 */
      k = sprintf(ln, "cat %s  >> %s", cb->rout, rt);
      assert(k < L);
      if (system(ln))
      {
         fprintf(stderr, "ERROR IN CMD: '%s'\n", ln);
         assert(0);
      }
   }
   free(ln);
}

void GenCopyTrmm(char pre, char *outd, char sd, char up, char ta,
                           ATL_mmnode_t *mb)
/*
 * LN/UT: tAT-L, tAN-U
 * LT/UN: tAN-L, tAT-U
 */
{
   int i, L, flg;
   char *ln, *rt;
   ATL_cpnode_t *cb;

/*
 * copy name format
 *    ATL_<pre>cp[Into,From]A[N,T][g,k,y,r,s,t]_aX_<ku>x<nu>_<kvec>
 */
   flg = (1<<CPF_TOBLK);
   flg |= (pre == 'd' || pre == 's') ? (1<<CPF_REAL):0;
   flg |= (pre == 'c' || pre == 's') ? (1<<CPF_SINGLE):0;
   if (sd == 'L')
      flg |= (1<<CPF_ABLK);

   cb = GetCopyNodeFromMM(flg, mb, NULL);
   assert(cb);
   cb->flag |= (1<<CPF_TRMM); /* set TRMM flag */
/*
 * up is always 'L'... ta can be 'N' or 'T'
 *
 * LN/UT: tAT-L, tAN-U
 * LT/UN: tAN-L, tAT-U
 */
   /*printf("Case: %c%c%c\n", sd, up, ta);*/
   if (ta == 'N')
   {
      GenTrCopy(mb, cb, outd, pre, sd, 'L', 'N');
      GenTrCopy(mb, cb, outd, pre, sd, 'U', 'T');
   }
   else
   {
      GenTrCopy(mb, cb, outd, pre, sd, 'L', 'T');
      GenTrCopy(mb, cb, outd, pre, sd, 'U', 'N');
   }
   KillAllCopyNodes(cb);
}

void GenMakeTrmm(char pre, char *outd, char up, char *sds, char *tas,
                 ATL_mmnode_t **MBs)
{
   int k, i, UR;
   FILE *fp;
   char *fn, *typ, *comp, *sp, *flgs;
   const char be[3] = {'0', '1', 'n'};
   const char al[3] = {'1', 'N', 'X'};
   int bi, ai, kvec;
   char tta;

   k = strlen(outd) + 24;
   fn = malloc(k);
   assert(fn);
   i = sprintf(fn, "%s/Makefile", outd);
   fp = fopen(fn, "w");
   assert(fp);
   free(fn);
   fprintf(fp, "include ../Make.inc\n");
   fprintf(fp, "CDEFS2=$(CDEFS)\n\n");

   fprintf(fp, "lib : %clib.grd\n", pre);
   fprintf(fp, "all : %clib.grd\n", pre);
   fprintf(fp, "clean : %cclean\n\n", pre);

   fprintf(fp, "obj = ");
   if (pre == 'c' || pre == 's')
   {
      comp = "$(SKC)";
      flgs = "$(SKCFLAGS)";
   }
   else
   {
      comp = "$(DKC)";
      flgs = "$(DKCFLAGS)";
   }
   if (pre == 'z')
      typ = "DCPLX";
   else if (pre == 'c')
      typ = "SCPLX";
   else if (pre == 's')
      typ = "SREAL";
   else
      typ = "DREAL";

   for (k=0; k < 2; k++)  /* 0: objname, 1: target*/
   {
/*
 *    kernels
 */
      int s, di;
      for (s=i=0; s < 2; s++)
      {
         int it;
         char sd=sds[s];
         for (it=0; it < 2; it++)
         {
            ATL_mmnode_t *mp;
            char ta = tas[it];
            for (mp=MBs[i++]; mp; mp=mp->next)
            {
               for (bi=0; bi < 3; bi++)
               {
                  if (!k)
                     fprintf(fp, "\\\n      ATL_%ctrmmK_%c%c%c_b%c.o ",
                             pre, sd, up, ta, be[bi]);
                  else /* target phase */
                  {
                     fprintf(fp, "ATL_%ctrmmK_%c%c%c_b%c.o :"
                             " ATL_%ctrmmK_%c%c%c_%ux%ux%u.c $(deps)\n",
                              pre, sd, up, ta, be[bi],
                              pre, sd, up, ta, mp->mu, mp->nu, mp->ku);
                     fprintf(fp, "\t%s $(CDEFS2) -D%s=1 -DBETA%c=1 \\\n", comp,
                           typ, (be[bi]=='n')?'N':be[bi]);
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
                     "        -DATL_USERMM=ATL_%ctrmmK_%c%c%c_b%c \\\n"
                     "        -c -o ATL_%ctrmmK_%c%c%c_b%c.o \\\n",
                             pre, sd, up, ta, be[bi],
                             pre, sd, up, ta, be[bi]);
                     fprintf(fp, "        ATL_%ctrmmK_%c%c%c_%ux%ux%u.c \n",
                             pre, sd, up, ta, mp->mu, mp->nu, mp->ku);
                  }
               }
            }
         }
      }
/*
 *    copies
 *    e.g., ATL_zcpFromACt_Up_diagU_a1_1x4_0
 */
      for (di=0; di < 2; di++)
      {
         char *dia = di ? "_diagU": "";
         for (s=i=0; s < 2; s++)
         {
            int it;
            char sd=sds[s];
            for (it=0; it < 2; it++)
            {
               ATL_mmnode_t *mp;
               char ta = tas[it];
               for (mp=MBs[i++]; mp; mp=mp->next)
               {
                  UR = (sd=='L')? mp->mu : mp->nu;
                  kvec = FLAG_IS_SET(mp->flag, MMF_KVEC) ? mp->vlen:0;
                  for (ai=0; ai < 3; ai++)
                  {
                     int ul;
                     for (ul=0; ul < 2; ul++)
                     {
                        int c;
                        const int NCNJ = ((pre == 'c' || pre == 'z')) ? 2:1;
                        for (c=0; c < NCNJ; c++)
                        {
                           if (c)
                              tta = (ta == 'N')? 'C' : 'H';
                           else
                              tta = ta;
                           if (!k)
                           {
                              fprintf(fp,
                              "\\\n      ATL_%ccpFromA%ct_%s%s_a%c_%ux%u_%u.o ",
                                 pre, tta, (ul?"Up":"Lo"), dia, al[ai], mp->ku,
                                 UR, kvec);
                           }
                           else
                           {
                              fprintf(fp,
                                      "ATL_%ccpFromA%ct_%s%s_a%c_%ux%u_%u.o :"
                                      " ATL_%ccpFromA%ct_a%c_%ux%u_%u.c\n",
                                      pre, tta, (ul?"Up":"Lo"), dia, al[ai],
                                      mp->ku, UR, kvec,
                                      pre, ta, al[ai], mp->ku, UR, kvec);

                              fprintf(fp, "\t%s $(CDEFS2) %s", comp, flgs);
                              if (di)
                                 fprintf(fp, " -DUnitDiag_=1");
                              if (ul)
                                 fprintf(fp, " -DUpper_=1");
                              else
                                 fprintf(fp, " -DLower_=1");
                              if (c)
                                 fprintf(fp, " -DConj_=1");
                              fprintf(fp,
                              " \\\n\t-c ATL_%ccpFromA%ct_a%c_%ux%u_%u.c \\\n",
                                      pre, ta, al[ai], mp->ku, UR, kvec);
                              fprintf(fp,
                              "\t-o ATL_%ccpFromA%ct_%s%s_a%c_%ux%u_%u.o\n",
                                      pre, tta, (ul?"Up":"Lo"), dia, al[ai],
                                      mp->ku, UR, kvec);
                           }
                        }
                     }
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

void GenAllTrmm(char pre, char *outd)
{
   int L, i, s;
   char sds[2] = {'L', 'R'};
   char tas[2] = {'N', 'T'};
   char ups[2] = {'L', 'U'};
   char *od;
   char fn[16];
   ATL_mmnode_t *MBs[4];
#if 1
/*
 * Delete and the re-create <pre>UTRSM subdir
 */
   L = strlen(outd) + 8 + 7;
   od = malloc(L);
   assert(od);
   i = sprintf(od, "rm -rf %s/%cUTRMM", outd, pre);
   assert(i < L);
   Sys2File(od, NULL);
   od[0]='m'; od[1]='k'; od[2]='d'; od[3]='i'; od[4]='r'; od[5]=' ';
   assert(!Sys2File(od, NULL));
#endif
   strcpy(fn, "trmmKSU.sum");

   for (i=s=0; s < 2; s++)
   {
      char sd=sds[s];
      int it;
      for (it=0; it < 2; it++)
      {
         char ta = tas[it];
         char up = ups[it];
         ATL_mmnode_t *mb;
         fn[5] = sd;
         fn[6] = up;
         mb = ReadMMFileWithPath(pre, "res", fn);
         if (!mb)
            fprintf(stderr, "Problem with file='res/%c%s'!\n", pre, fn);
         assert(mb);
         mb->blask = ATL_KTRMM; /* forced to TRMM */
         MBs[i] = mb;
         GenMicroTrmm(pre, od+7, sd, 'L', ta, mb);
         GenCopyTrmm(pre, od+7, sd, 'L', ta, mb);
         GenTrmmHead(pre, od+7, sd, 'L', ta, mb);
         i++;
      }
   }
   GenMakeTrmm(pre, od+7, 'L', sds, tas, MBs);
   for (i=0; i < 4; i++)
      KillAllMMNodes(MBs[i]);
   free(od);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   /*fprintf(stderr,
"This will create the <pre>UTRMM subdir in <outdir>, and populate it with\n"
"all utrmm kernels and their Makefile.\n");*/
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
   GenAllTrmm(pre, outd);
   free(outd);
   return(0);
}
