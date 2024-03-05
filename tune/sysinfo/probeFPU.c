#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   #define GAS_X86 1
   #define SSE    0
   #define AVX    1
   #define AVXFMA 2
   #define NISA   3
   char *isas[NISA] = {"SSE", "AVX", "AVXFMA"};
#else
   #define GAS_X86 0
   #if defined(ATL_GAS_ARM64)
      char *isas[1] = {"ARM64"};
      #define ARM64 0
      #define NISA 1
   #elif defined(ATL_GAS_ARM)
      char *isas[1] = {"ARM32"};
      #define ARM 0
      #define NISA 1
   #elif defined(ATL_GAS_PPC)
      #define NISA   2
      #define VSX 0
      #define PWR 1
      char *isas[NISA] = {"VSX", "POWER"};
   #else
      #error "ISA not presented supported!"
   #endif
#endif
#define LNLEN 128
double getMflop(char pre, int isa, int vlenb, int nacc, int nmul, double mf)
{
   char ln[LNLEN];
   int i;
   FILE *fpin;
   i = sprintf(ln, "make fpucase isa=%s nacc=%d nmul=%d vlen=%d mf=%f "
               "typ=\"%s\"",
               isas[isa], nacc, nmul, vlenb, mf > 0.0?mf:0.0,
               pre=='s'?"-DSREAL=1":"");
   assert(i < LNLEN);
   if (system(ln))
      return(-1.0);
   fpin = fopen("res/fpuS.out", "r");
   if (!fpin)
      return(-2.0);
   if (fscanf(fpin, "%le", &mf) != 1)
      return(-3.0);
   fclose(fpin);
   return(mf);
}

int findISA(char pre, int *VLENB)
{
   #if defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   int i, vl;
   printf("\nPROBING FOR ISA:\n");
   for (i=NISA-1; i >= 0; i--)
   {
      if (getMflop(pre, i, 16, 4, i==AVXFMA?0:1, 0.0) > 0.0)
      {
         printf("   %s : yes!\n", isas[i]);
         break;
      }
      else
         printf("   %s : no.\n", isas[i]);
   }
   if (*VLENB == 0)
   {
      printf("\nPROBING FOR VLEN (bytes):\n");
      for (vl=64; vl >= 32; vl >>= 1)
      {
         if (getMflop(pre, i, vl, 4, 1, 0.0) > 0.0)
         {
            printf("   VLEN=%u: yes!\n", vl);
            break;
         }
         else
            printf("   VLEN=%u: no.\n", vl);
      }
      *VLENB = vl;
   }
   return(i);
   #elif defined(ATL_AVXZ)
      *VLENB = 64;
      return(AVXFMA);
   #elif defined(ATL_AVXMAC)
      *VLENB = 32;
      return(AVXFMA);
   #elif defined(ATL_AVX)
      *VLENB = 32;
      return(AVX);
   #elif defined(ATL_SSE1)
      *VLENB = 16;
      return(SSE);
   #elif defined(ATL_GAS_PPC)
      *VLENB = 16;
      return(VSX);
   #elif defined(ATL_GAS_ARM64)
      *VLENB = 16;
      return(ARM64);
   #elif defined(ATL_GAS_ARM)
      *VLENB = 0;
      return(ARM);
   #endif
   *VLENB = 0;
   return(-1);
}

int findNREG(char pre, int isa, int vlenb)
{
   #if GAS_X86
      const int nmul = (isa == AVXFMA) ? 0:1, nextra = nmul+1;;
   #else
      const int nmul=0, nextra=1;;
   #endif
   int nreg=128;
   printf("PROBING FOR NUMBER OF REGISTERS WITH VLENB=%u:\n", vlenb);
   printf("   TRYING NREG=%d: ", nreg);
   while (getMflop(pre, isa, vlenb, nreg-nextra, 1, 0.0) <= 0.0)
   {
      nreg >>= 1;
      printf("NO.\n   TRYING NREG=%d: ", nreg);
      assert(nreg > 0);
   }
   printf("YES!\nNUMBER OF REGISTERS DETECTED AS: %d\n\n", nreg);
   return(nreg);
}

int findLAT(char pre, int isa, int vlenb, int nreg, int *NMUL, double *PEAK)
/*
 * Finds: peak mflop and effective latency of FMA of selected isa & vlenb.
 * If this ISA has no FMA (separate mul & add) also finds how many regs must
 * be dedicated to multiply registers.  If machine supports register renaming
 * in hardware (modern x86), nmul almost always 1, but w/o it will usually
 * the stages in the fmul pipe.  This code will detect FMA latency as max of
 * add and mul latency, usually, so not optimal solution for separate mul&add,
 * but these are not common in modern architecture anyway.
 * RETURNS: 0 if selected code will not compile, else detected FMA latency.
 */
{
   #if GAS_X86
      int NOFMA=(isa == AVXFMA) ? 0:1;
   #else
      int NOFMA=0;
   #endif
   int nmulB, naccB, nacc, nmul;
   double mfB=0.0, mf, mfF=(*PEAK);

   printf("PROBING FOR EFFECTIVE FMA LATENCY, VLEN=%u:\n", vlenb);
/*
 * If no FMA, then try maximum NMUL setting that will fit in nreg:
 *    nmul+nacc+1 <= nreg, nmul=nacc = (nreg-1)/2;
 */
   if (NOFMA)
   {
      nacc = (nreg-2)>>1;
      nacc = (nacc < 1) ? 1:nacc;
      mfB = getMflop(pre, isa, vlenb, nacc, nacc, mfF);
      naccB = nmulB = nacc;
      printf("   nacc=%d, nmul=%d, mf=%.2f\n", nacc, nacc, mfB);
   }
   else
      *NMUL = nmul = 0;
   nacc = nreg-2;
   mf = getMflop(pre, isa, vlenb, nacc, NOFMA, mfF);
   printf("   nacc=%d, nmul=%d, mf=%.2f\n", nacc, nmul, mf);
   if (mf > mfB)
   {
      mfB = mf;
      naccB = nacc;
      nmulB = NOFMA;
   }
   if (mfB <= 0.0)
   {
      *PEAK = 0.0;
      *NMUL = 0;
      return(0);
   }
   nacc = naccB;
   do
   {
      int nacc0=nacc;
      nacc = (nacc > 16) ? (nacc>>1):nacc-1;
      nmul = NOFMA;
      mf = getMflop(pre, isa, vlenb, nacc, 1, mfF);
      printf("   nacc=%d, nmul=%d, mf=%.2f\n", nacc, nmul, mf);
      if (NOFMA && nacc + nacc + 1 <= nreg)
      {
         double mfM;
         mfM = getMflop(pre, isa, vlenb, nacc, nacc, mfF);
         printf("   nacc=%d, nmul=%d, mf=%.2f\n", nacc, nacc, mfM);
         if (mfM > mf)
         {
            nmul = nacc;
            mf = mfM;
         }
      }
      if (mf*1.03 <= mfB)  /* definite drop-off in performance */
      {
         naccB = nacc0;
         break;
      }
      else
      {
         mfB = (mf > mfB) ? mf : mfB;
         naccB = nacc;
         nmulB = nmul;
      }

   }
   while (naccB > 1);
/*
 * If best latency comes from halving size, need to make estimate exact
 */
   if (naccB >= 16)
   {
      double mf, mfL=mfB;
      unsigned int big=naccB, sml=naccB>>1;
      while (big-sml > 1)
      {
         unsigned int mid = sml + ((big-sml)>>1);
         mf = getMflop(pre, isa, vlenb, mid, mid, mfF);
         printf("   nacc=%d, nmul=%d, mf=%.2f\n", mid, mid, mf);
         if (mf*1.03 <= mfB)  /* midpt still too slow, adjust upwards */
            sml = mid;
         else /* need to adjust big downwards */
            big = mid;
      }
   }
   printf("EFFECTIVE FMA LATENCY DETECTED AS: %d, NMUL=%d\n\n", naccB, nmulB);
   *PEAK = mfB;
   *NMUL = nmulB;
   return(naccB);
}

void findAllPeaks(char pre, int isa, int nreg, int nmul, double peak)
{
   double mfF;
   int i, vlens[3] = {64, 32, 16};
   char fn[16];
   FILE *fp;

   mfF = .05*peak;
   sprintf(fn, "res/%cSIMD.idx", pre);
   fp = fopen(fn, "w");
   assert(fp);
   for (i=isa; i >= 0; i--)
   {
      int v;
      for (v=0; v < 3; v++)
      {
         int lat, vl=vlens[v], nmul;
         double mfpk;
         lat = findLAT(pre, i, vl, nreg, &nmul, &mfpk);
         fprintf(fp, "'%s' %d %.2f %d %d\n",  isas[i], vl, mfpk, lat, nmul);
         printf("******************************\n");
         printf("'%s' %d %.2f %d %d\n",  isas[i], vl, mfpk, lat, nmul);
         printf("******************************\n\n");
      }
   }
   fclose(fp);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p <pre> [d] : Date type prefix, one of [s,d]\n");
   fprintf(stderr, "   -v <vlen> [probe]: vector len in bytes to force\n");
   exit(ierr ? ierr : -1);
}

char getFlags(int nargs, char **args, int *VLEN)
{
   char pre='d';
   int i;

   *VLEN = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'p':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         pre = args[i][0];
         break;
      case 'v':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *VLEN = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (pre == 'S' || pre == 's' || pre == 'C' || pre == 'c')
       pre = 's';
   else
      pre = 'd';
   return(pre);
}
/*
 * This probe finds the minimal vector operation necessary to get peak
 * It presently only works for AMD64.
 */
int main(int nargs, char **args)
{
   double peak=0.0;
   FILE *fp;
   #if GAS_X86
      int FMA;
   #else
      int FMA=1;
   #endif
   int isa, nreg, vlenb, nlat, nmul;
   char pre='d';
   char fn[32];

   pre = getFlags(nargs, args, &vlenb);

   isa = findISA(pre, &vlenb);
   #if GAS_X86
      FMA = (isa == AVXFMA);
   #endif
   nreg = findNREG(pre, isa, vlenb);
   assert(nreg >= 0);
   nlat = findLAT(pre, isa, vlenb, nreg, &nmul, &peak);
   assert(nlat > 0 && nlat < nreg && peak > 0.0);
/*
 * For backwards compatibility, write <pre>MULADD file
 */
   sprintf(fn, "res/%cMULADD", pre);
   fp = fopen(fn, "w");
   assert(fp);
   fprintf(fp, "%d\n%d\n%.2f\n%d\n", FMA, nlat, peak, nreg);
   fprintf(fp, "#'<simd>' <mflop> <nreg> <vlen bytes> <FMA latency> <nmul>\n");
   fprintf(fp, "'%s' %.2f %u %u %u %u\n", isas[isa], peak, nreg, vlenb,
           nlat, nmul);
   fclose(fp);
/*
 * Confirm peak performance of best case using a forced mflop high enough
 * we think timing is repeatable.
 * getMflop call ensures last generated fpu.S is an optimal one.
 */
   peak = getMflop(pre, isa, vlenb, nlat, nmul, .1*peak);
   fp = fopen("fnlfpu.S", "w");
   assert(fp);
   fprintf(fp, "#define VL%u 1\n", vlenb);
   fclose(fp);
   assert(!system("cat fpu.S >> fnlfpu.S"));
/*
 * Store results in res/<pre>SIMD
 */
   fn[5] = 'S'; fn[6] = 'I'; fn[7] = 'M'; fn[8] = 'D'; fn[9] = '\0';
   fp = fopen(fn, "w");
   assert(fp);
   fprintf(fp, "'%s' %.2f %u %u %u %u\n", isas[isa], peak, nreg, vlenb,
           nlat, nmul);
   fprintf(fp, "#'<simd>' <mflop> <nreg> <vlen bytes> <FMA latency> <nmul>\n");
   fclose(fp);
/*
 * Ask for peak code to take .05 seconds in <pre>mflop.frc
 */
   sprintf(fn, "../blas/gemm/%cmflops.frc", pre);
   fp = fopen(fn, "w");
   assert(fp);
   fprintf(fp, "%le", .05*peak);
   fclose(fp);

   #if 0
      findAllPeaks(pre, isa, nreg, nmul, peak);
   #endif
   return(0);
}
