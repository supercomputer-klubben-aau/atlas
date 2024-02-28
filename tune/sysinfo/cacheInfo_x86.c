/*
 * Wrote this from gcc manual + http://www.sandpile.org/x86/cpuid.htm
 * as well as AMD's 25481.pdf (CPUID Specification), and some Intel docs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define KSH 10 /* kilo-byte (1024) shift */
#define MSH 20 /* mega-byte (104857) shift */

#define do_cpuid(eax_, ebx_, ecx_, edx_)  __asm__("cpuid" \
           :"=a"(eax_), "=b"(ebx_), "=c"(ecx_), "=d"(edx_) \
           :"a"(eax_), "c"(ecx_) \
           : )

#define uint unsigned int
#define NONE 0xFFFFFFFF
unsigned int nways_amd(unsigned int t)
/*
 * AMD trolling hard by hardwiring in partial lookup table for associativity
 * and doing it differently for every level.  This is the L3 lookup.
 * RETURNS: nways, or 0 if fully associative
 */
{
   unsigned int nway;
   switch (t)     /* AMD asks: why not hardwire assoc in lookup table? */
   {
   case 0:        /* no such level */
      return(NONE);  /* return 0 from caller, don't set any ptrs */
   case 0x1:
   case 0x2:
   case 0x4:
      nway = t;
      break;
   case 0x6:
      nway = 8;
      break;
   case 0x8:
      nway = 16;
      break;
   case 0xA:
      nway = 32;
      break;
   case 0xB:
      nway = 48;
      break;
   case 0xC:
      nway = 64;
      break;
   case 0xD:
      nway = 96;
      break;
   case 0xE:
      nway = 128;
      break;
   case 0xF:
      nway = 0;
      break;
   default:
      assert(0);
   }
   return(nway);
}
int IsPow2(uint N)
/*
 * returns 1<<i for power of two, or 0 zero or non-power-of-two
 */
{
   int i;
   if (N == 0)
      return(0);
   for (i=0; i & (1<<i); i++); /* find first nonzero */
   return((((1<<i)^N) == 0) ? i : 0);
}

unsigned int GetNumPhysCore_amd(void)
/*
 * Returns the number of physical cores on AMD machines supporting 0x80000008.
 * Will be double the number of FPUs on later archs, unfortunately.
 * NOTE: already tested 0x80000008 supported before calling.
 */
{
   unsigned int eax, ebx, ecx, edx;
   eax = 0x80000008;
   do_cpuid(eax, ebx, ecx, edx);
   return((ecx&0XFF)+1);
}

unsigned int SafeGetNumPhysCore_amd(void)
{
   unsigned int eax, ebx, ecx, edx;
   eax = 0x80000000;
   do_cpuid(eax, ebx, ecx, edx);
   if (eax >= 0x80000008)
   {
      eax = 0x80000008;
      do_cpuid(eax, ebx, ecx, edx);
      return((ecx&0XFF)+1);
   }
   return(0);
}

uint cacheInfo_amd(uint lvl, uint *NSET, uint *NWAY, uint *INCL,
                   uint *NTHR, uint *NCOR)
{
   unsigned int eax, ebx, ecx, edx, LS, nway, nset, t, LPT;
   unsigned long cs;
/*
 * AFAIK, L1 & L2 always private, L3 shared for AMD
 */
   if (lvl > 2)  /* AMD CPUID presently only allows up to 3 lvls cache */
      return(0);
   ecx = 0;
   eax = 0x80000005 + (lvl != 0);
   do_cpuid(eax, ebx, ecx, edx);
   if (!lvl)  /* hard to overstate the madness of AMD's cpuid: */
   {          /* they've hardwired the # of cache lvls in! */
      *INCL = 1;
      LS = ecx & 0x7F;
      LPT = (ecx>>8)&0xFF;
      t = (edx>>16)&0xFF;  /* way encoding */
      nway = (t == 0x01) ? 1 : (t>>1); /* will fix fully assoc later */
      cs = (ecx>>24)*1024;
      if (t == 0xFF)   /* fully associative */
      {
         nway = cs / LS;
         nset = 1;
      }
      else
         nset = cs / (LS*nway);
   }
   else /* L2:ecx, L3:edx */
   {
      t = (lvl == 1) ? ecx : edx;
      nway = nways_amd((t>>12)&0xF);
      if (nway == NONE)  /* if this level of cache not there/disabled */
         return(0);      /* let caller know */
      LS = t & 0xFF;
      LPT = (t>>8)&0xF;
      if (lvl == 1)           /* L2 deterministic unlike L3 on cache size */
         cs = (ecx>>16)*1024;
/*
 *    For L3, we get i*512KB <= CS <= (i+1)*512KB.  For fully-associative,
 *    just declare it as the smaller size.  For a given associativity, we
 *    know nset is a power of two, which can further constrain CS.
 */
      else                    /* L3 doesn't fully specify cache size */
      {                       /* guess low bound unless it doesn't work */
         cs = (edx>>18)*512*1024; /* minCS is i*512KB, max (i+1)*512KB */
         if (nway)  /* we can determine nsets */
         {
            nset = cs / nway;
            if (nset*nway != cs) /* if nway not mul of cs, cs not right! */
               cs += 512*1024;
            else if (!IsPow2(nset)) /* if nsets not pow2, CS not right! */
               cs += 512*1024;
         }
      }
      if (!nway)              /* full assoc */
      {
         nset = 1;
         nway = cs / LS;
      }
      else                    /* cache with nway-associativity */
         nset = cs / (LS*nway);
   }
/*
 * AMD has "lines per tag" field, and they never explain what it is beyond
 * the phrase itself.  AFAIK, lines per tag is 1 in every cache ever made.
 * only way >1 makes sense to me is if you treat cache like a invalidating
 * write buffer, where you can load only part of the cacheline, and indicate
 * some bytes in the line are invalid.  Then, this would represent the number
 * of valid bits you use.  However, I've never heard of a cache working like
 * that, and all the AMD docs I see say have this as 0 (invalid cache) or 1.
 * so, I'm going to assert these cases, and force a probe if this isn't true,
 * since I'm not sure I understand what LPT really is.
 */
   assert(LPT == 0 || LPT == 1);
   *NWAY = nway;
   *NSET = nset;
/*
 * AFAIK, L1&L2 always private, L3 public for AMD, with all caches except L1
 * essentially being victim caches (i.e., not inclusive, though not truly
 * exclusive either).
 */
   *INCL = (lvl == 0) ? 1:0;
   *NTHR = (lvl == 2) ? 2 : 1; /* AFAIK, L1 & L2 private, L3 shared for AMD */
   *NCOR = SafeGetNumPhysCore_amd();
   return(LS);
}

int do_cpuid2Lvl_amd(uint lvl, uint *EAX, uint *EBX, uint *ECX, uint *EDX)
/*
 * Loops over AMD's levels, until lvl=lvl in data cache is achieved (L1->lvl=0)
 * RETURNS: -1 for no such level, else the actual lvl needed by cpuid
 */
{
   int i;
   int id, tlvl=lvl, glvl, eax, ebx, ecx, edx;
   do
   {
      eax = 0x8000001D;
      ecx = tlvl++;                  /* keep increasing try lvl til one works */
      do_cpuid(eax, ebx, ecx, edx);
      id = eax&0x1F;                 /* cache type info returned about */
      glvl = (eax>>5)&7;             /* what cache lvl did I get? */
printf("TRY:id=%u, lvlret=%u, lvlwnt=%u\n", id, glvl, lvl+1);
      if (!id)                       /* if we are out of cache levels */
         return(-1);                 /* tell caller no match */
   }
   while((id != 1 && id != 3) || glvl != lvl+1);
   *EAX = eax;
   *EBX = ebx;
   *ECX = ecx;
   *EDX = edx;
   printf("OK! :id=%u, lvlret=%u, lvlwnt=%u\n", id, glvl, lvl+1);
   return(glvl);
}

uint cacheInfo_amdD(uint lvl, uint *NSET, uint *NWAY, uint *INCL,
                    uint *NTHR, uint *NCOR)
/*
 * This is AMD's more intel-like cache interface, works on newer procs.
 * RETURNS: 0 for no such level, else line size
 */
{
   unsigned int eax, ebx, ecx, edx, LS, lv;
   int i;

/*
 * If we are past the number of cache levels, return 0
 */
   if (do_cpuid2Lvl_amd(lvl, &eax, &ebx, &ecx, &edx) == -1)
      return(0);

   *NCOR = GetNumPhysCore_amd();
   *NTHR = ((eax>>14)&0xFFF) + 1;
   *NWAY = (ebx>>22) + 1;
   *NSET = ecx + 1;
   *INCL = (edx>>1)&1;
   LS = (ebx & 0xFFF) + 1;
   return(LS);
}
#define do_cpuid_int(lvl_, o1_, o2_, o3_, o4_)  __asm__("cpuid" \
           :"=a"(o1_), "=b"(o2_), "=c"(o3_), "=d"(o4_) \
           :"a"(0x00000004), "c"(lvl_) \
           : )
int tlbLookup(unsigned char ch, int *N, long *SZ)
{
   unsigned lv=0;  /* 0 means not about TLB, else lvl of data TLB described */
   switch(ch)
   {
   case 0x3:
      *N = 64;
      *SZ = 4<<KSH;
      return(1);
   case 0x4:
      *N = 8;
      *SZ = 4<<MSH;
      return(1);
   case 0x5:
      *N = 32;
      *SZ = 4<<MSH;
      return(2);
   case 0x56:
      *N=16;
      *SZ = 4<<MSH;
      return(1);
   case 0x57:
   case 0x59:
      *N=16;
      *SZ = 4<<KSH;
      return(1);
   case 0x5A:
      *N=32;
      *SZ = 4<<MSH;
      return(1);
   case 0x5B:
      *N=64;
      *SZ = 4<<MSH;
      return(1);
   case 0x5C:
      *N=128;
      *SZ = 4<<MSH;
      return(1);
   case 0x5D:
      *N=256;
      *SZ = 4<<MSH;
      return(1);
   case 0x63:
      *N=32+4;
      *SZ = 4<<MSH;
      return(1);
   case 0x64:
      *N=512;
      *SZ = 4<<KSH;
      return(1);
   case 0xa0:
      *N=32;
      *SZ = 4<<KSH;
      return(1);
   case 0xb3:
      *N=128;
      *SZ = 4<<KSH;
      return(1);
   case 0xb4:
      *N=256;
      *SZ = 4<<KSH;
      return(2);
   case 0xba:
      *N=64;
      *SZ = 4<<KSH;
      return(2);
   case 0xc0:
      *N=8;
      *SZ = 4<<MSH;
      return(1);
   case 0xc1:
      *N=1024;
      *SZ = 2<<MSH;
      return(-2);  /* shared L2 TLB */
   case 0xc2:
      *N=16;
      *SZ = 2<<MSH;
      return(0);
   case 0xc3:
      *N=1536+16;
      *SZ = 2<<MSH;
      return(-2);  /* shared L2 TLB */
   case 0xc4:
      *N=32;
      *SZ = 2<<MSH;
      return(1);
   case 0xca:
      *N=512;
      *SZ = 4<<KSH;
      return(-3);  /* shared L3 TLB */
   default:
      return(0);
   }
   return(0);
}
uint tlbInfo_int(uint *NPRV, uint **NTLBs)
/*
 * *NPRIV: number of private levels: *NPRIV <= return
 * *NTLB : ptr to array of TLB entry counts for all levels
 * RETURNS: ntlb levels, and thus number of entries in *NTLB
 */
{
   unsigned int i=0, nprv=0, nlvl=0;
   static int Ns[3] = {0,0,0};
   unsigned int regs[4];
   unsigned char *byts=(unsigned char*)regs;

   do
   {
      unsigned int k;

      regs[0] = 2;
      do_cpuid(regs[0], regs[1], regs[2], regs[3]);
      for (k=1; k < 16; k++)
      {
         int lv, ntlb;
         long pgsz;   /* value unused for now */
         lv = tlbLookup(byts[k], &ntlb, &pgsz);
         if (!lv)
            continue;
         if (lv > 0)
            nprv = (lv > nprv) ? lv:nprv;
         else
            lv = -lv;
         assert(lv < 4);
         nlvl = (lv > nlvl) ? lv:nlvl;
         Ns[lv-1] += ntlb;
      }
   }
   while (++i < byts[0]);

   *NTLBs = Ns;
   *NPRV = nprv;
   return(nlvl);
}

uint cacheInfo(uint lvl, uint *NSET, uint *NWAY, uint *INCL,
               uint *NTHR, uint *NCOR)
/*
 * RETURNS: 0 for no such level, else line size
 */
{
   unsigned int eax, ebx, ecx, edx, LS, lv;

   do_cpuid_int(lvl, eax, ebx, ecx, edx);
/*
 * If we are past the number of cache levels, return 0
 */
   if ((eax&0xF) == 0)
      return(0);
/*
 * On some systems, lvl=1 is L1 instruction cache, so check if the level
 * returned (starting from 0) doesn't match the query, and if so query next
 * level, and then assert it matches
 */
   lv = (eax>>5)&7;
   if (lv != lvl+1)  /* lvl=1 is i-cache, meaning we must ask for next lvl */
   {
      do_cpuid_int(lvl+1, eax, ebx, ecx, edx);
      if ((eax&0xF) == 0)
         return(0);
      lv = (eax>>5)&7;
   }
   assert(lv == lvl+1);

   *NCOR = ((eax>>26)&0x3F) + 1;
   *NTHR = ((eax>>14)&0xFFF) + 1;
   *NSET = ecx + 1;
   *NWAY = (ebx>>22) + 1;
   *INCL = (edx>>1)&1;
   LS = (ebx & 0xFFF) + 1;
   return(LS);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n", ierr,
              flag?flag:"OUT-OF_ARGUMENTS");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr,"USAGE: %s [flags]:\n", name);
   fprintf(stderr,"   -o </path/cacheheader.h> : (res/atlas_cache.h)\n");
   exit(ierr ? ierr : -1);
}

FILE *GetFlags(int nargs, char **args)
{
   FILE *fpout=NULL;
   int i;

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);

      switch(args[i][1])
      {
      case 'o':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        fpout = fopen(args[i], "w");
        assert(fpout);
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   if (!fpout)
   {
      fpout = fopen("res/atlas_cache.h", "w");
      assert(fpout);
   }
   return(fpout);
}

void *getInf(void)
/*
 * Determines whether to use AMD or Intel method to determine caches
 */
{
   int eax, ebx, ecx, edx;
   eax=0;
   do_cpuid(eax, ebx, ecx, edx);
   printf("cpu0: %x %x %x %x\n", eax, ebx,edx,ecx);

   if (ebx == 0x756e6547 && edx == 0x49656e69 && ecx == 0x6c65746e &&
       eax >= 4)
      return(cacheInfo);  /* Spells out GenuineIntel in little endian hex */
   if (ebx == 0x68747541 && edx == 0x69746e65 && ecx == 0x444d4163)
   {  /* Spells out AuthenticAMD in little endian hex */
      eax = 0x80000000;
      do_cpuid(eax, ebx, ecx, edx);
      printf("eax=%x want=%x, gap=%d\n", eax, 0x8000001D, eax-0x8000001D);
      if (eax >= 0x8000001D)
         return(cacheInfo_amdD);
      else if (eax >= 0x80000006)
         return(cacheInfo_amd);
   }
   return(NULL);
}
/*
 * This file either uses CPUID to determine the caches, or dies with
 * assertion failed.  If this routine fails to run or build, we instead
 * run cachesearch to empirically estimate main cache features.
 */
int main(int nargs, char **args)
{
   unsigned int nset, nway, LS, INCL, nthr, ncor=0, nlvl=1;
   unsigned int nhyp=1, safeLS, LLP, i;
   unsigned long cs, llp_cs, llp_prev_cs;
   FILE *fp;
   char *cn[2]={"LLPC", "LLC"};
   int lvls[2];
   uint (*cacheInf)(uint lvl, uint *NSET, uint *NWAY, uint *INCL,
                   uint *NTHR, uint *NCOR);

/*
 * See if we are using Intel or AMD-style cacheInfo, or giving up
 */
   cacheInf = getInf();
   assert(cacheInf != NULL);
/*
 * OK, on system where this method should work, create atlas_cache.h
 */
   fp = GetFlags(nargs, args);
   fprintf(fp, "/* generated by ATLAS/tune/sysinfo/cacheInfo_x86.c */\n");
   fprintf(fp, "#ifndef ATLAS_CACHE_H\n   #define ATLAS_CACHE_H 1\n");
/*
 * For L1Size, figure out number of cores that are really hyperthreading
 */
   LS = cacheInf(0, &nset, &nway, &INCL, &nhyp, &ncor);
   assert(LS);
   llp_prev_cs = llp_cs = cs = nset*nway*LS;
   printf("L1: INCL=%d nsets=%u, nway=%u LS=%u sz=%u"
          ", hypthr=%u, ncores=%u\n",
          INCL, nset, nway, LS, nset*nway*LS, nhyp, ncor);
   fprintf(fp, "#define x86_NHYPTHR %u /* hyperthreading number */\n", nhyp);
   fprintf(fp, "#define x86_NCORPAK %u /* ncores per package */\n\n", ncor);
   fprintf(fp, "#define L1C_SZ %luL\n", cs);
   fprintf(fp, "#define L1C_LS %u\n", LS);
   fprintf(fp, "#define L1C_NSET %u\n", nset);
   fprintf(fp, "#define L1C_NWAY %u\n", nway);
   fprintf(fp, "#define L1C_INCLUSIVE 1\n");
   fprintf(fp, "#define L1C_NSHARE 1\n");
   fprintf(fp, "#define L1C_SURE 1\n");
   fprintf(fp, "#ifdef SREAL\n   #define L1C_ELTS %lu\n", cs>>2);
   fprintf(fp, "   #define L1C_LSELTS %u\n", LS>>2);
   fprintf(fp, "#elif defined(DREAL) || defined(SCPLX)\n   "
               "#define L1C_ELTS %lu\n", cs>>3);
   fprintf(fp, "   #define L1C_LSELTS %u\n", LS>>3);
   fprintf(fp, "#elif defined(DCPLX)\n   #define L1C_ELTS %lu\n", cs>>4);
   fprintf(fp, "   #define L1C_LSELTS %u\n", LS>>4);
   fprintf(fp, "#endif\n\n");
   safeLS = LS;
   LLP = 0;

/*
 * First find the number of levels, and whether there is a hyperthreading
 * multiplier on the ncores sharing cache or not
 */
   nlvl=0;
   while ( (LS = cacheInf(nlvl, &nset, &nway, &INCL, &nthr, &ncor)) )
   {
      unsigned int nshar;

      nlvl++;
      if (nhyp > 1)
      {
         nshar = nthr / nhyp;
         assert(nshar*nhyp == nthr);
      }
      else
         nshar = nthr;
      cs = nset*nway*LS;
      printf("L%u: INCL=%d nsets=%u, nway=%u LS=%u sz=%lu (%lu) nshar=%u\n",
             nlvl, INCL, nset, nway, LS, cs, cs/nshar, nshar);
      if (nshar > 1)  /* we want per-node size */
         cs /= nshar;
      fprintf(fp, "#define L%uC_SZ %luL\n", nlvl, cs);
      fprintf(fp, "#define L%uC_LS %u\n", nlvl, LS);
      fprintf(fp, "#define L%uC_NSET %u\n", nlvl, nset);
      fprintf(fp, "#define L%uC_NWAY %u\n", nlvl, nway);
      if (nlvl > 1)
         fprintf(fp, "#define L%uC_INCLUSIVE %u\n", nlvl, INCL);
      fprintf(fp, "#define L%uC_NSHARE %u\n", nlvl, nshar);
      fprintf(fp, "#define L%uC_SURE 1\n", nlvl);
      fprintf(fp, "#ifdef SREAL\n   #define L%uC_ELTS %luL\n", nlvl, cs>>2);
      fprintf(fp, "   #define L%uC_LSELTS %u\n", nlvl, LS>>2);
      fprintf(fp, "#elif defined(DREAL) || defined(SCPLX)\n   "
                  "#define L%uC_ELTS %luL\n", nlvl, cs>>3);
      fprintf(fp, "   #define L%uC_LSELTS %u\n", nlvl, LS>>3);
      fprintf(fp, "#elif defined(DCPLX)\n   #define L%uC_ELTS %luL\n",
              nlvl, cs>>4);
      fprintf(fp, "   #define L%uC_LSELTS %u\n", nlvl, LS>>4);
      fprintf(fp, "#endif\n\n");
      safeLS = (LS > safeLS) ? LS : safeLS;
      if (nthr <= nhyp )
      {
         llp_prev_cs = llp_cs;
         LLP = nlvl;
         llp_cs = cs;
      }
   }
   printf("\nLAST PRIVATE LEVEL = %u\n", LLP);
   fprintf(fp, "#define LLPC_LVL %u\n", LLP);
   fprintf(fp, "#define LLPC_SZ  L%uC_SZ\n", LLP);
   fprintf(fp, "#define LLPC_LSELTS L%uC_LSELTS\n", LLP);
   fprintf(fp, "#define LLC_LVL %u\n", nlvl);
   fprintf(fp, "#define LLC_LSELTS L%uC_LSELTS\n", nlvl);
   if (INCL)
   {
      fprintf(fp, "#define LLC_SZ   L%uC_SZ\n", nlvl);
      fprintf(fp, "#define LLC_ELTS L%uC_ELTS\n", nlvl);
   }
   else
   {
      unsigned long sz;
      sz = (LLP == nlvl) ? cs+llp_prev_cs : cs+llp_cs;
      fprintf(fp, "#define LLC_SZ %luL\n", cs+llp_prev_cs);
      fprintf(fp, "#ifdef SREAL\n   #define LLC_ELTS %luL\n", sz>>2);
      fprintf(fp, "#elif defined(DREAL) || defined(SCPLX)\n   "
                  "#define LLC_ELTS %luL\n", sz>>3);
      fprintf(fp, "#elif defined(DCPLX)\n   #define LLC_ELTS %luL\n",
              sz>>4);
      fprintf(fp, "#endif\n\n");
   }
   fprintf(fp, "#define ATL_SAFELS %u\n", safeLS);
   {
      unsigned int *tlbs, ntlbL, nprv, ttlb=0, tptlb=0;

      ntlbL = tlbInfo_int(&nprv, &tlbs);
      fprintf(fp, "\n#define TLB_PLVLS %u\n", nprv);
      fprintf(fp, "#define TLB_LVLS %u\n", ntlbL);
      for (i=0; i < ntlbL; i++)
      {
         fprintf(fp, "#define NTLB%u %u\n", i, tlbs[i]);
         tptlb += (i < nprv) ? tlbs[i] : 0;
         ttlb += tlbs[i];
      }
      fprintf(fp, "#define TOTPTLB %u\n", tptlb);
      fprintf(fp, "#define TOTTLB %u\n", ttlb);
   }

   fprintf(fp, "\n#endif /* done multiple inclusion guard */\n");
   fclose(fp);
   return(0);
}
