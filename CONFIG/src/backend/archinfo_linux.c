#include "atlconf.h"

void PrintUsage(char *name, int i)
{
   fprintf(stderr, "USAGE: %s -v (verb) -b (@ bits) -a (arch) -n (ncpu) -c <ncache> -C <lvl> (cache size) -m (Mhz) -t (cpu throttling)\n", name);
   exit(i);
}

int GetFlags(int nargs, char **args, int *CacheLevel)
{
   int i, flag = 0;

   *CacheLevel = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0], i);
      switch(args[i][1])
      {
      case 'n':
         flag |= Pncpu;
         break;
      case 'c':
         flag |= Pncache;
         break;
      case 'C':
         if (++i > nargs)
            PrintUsage(args[0], i);
         *CacheLevel = atoi(args[i]);
         break;
      case 'v':
         flag |= Pverb;
         break;
      case 'm':
         flag |= PMhz;
         break;
      case 'a':
         flag |= Parch;
         break;
      case 'b':
         flag |= P64;
         break;
      case 't':
         flag |= Pthrottle;
         break;
      default:
         PrintUsage(args[0], i);
      }
   }
   if (!flag)
     flag = Parch | P64;
   return(flag);
}

/*
 * RETURNS: CPU part
 * -1 means value is not found
 */
#ifndef uint
   #define uint unsigned int
#endif
uint getArmInfo(uint *IMPL, uint *ARCH, uint *VAR, uint *REV)
{
   uint part;
   char *res;
   *IMPL = *ARCH = *VAR = *REV = part = -1;
   res = atlsys_1L(NULL, "grep -F 'CPU part' /proc/cpuinfo", 0, 0);
   if (res)
   {
      char *sp;
      uint i;
      sp = strstr(res, "0x");
      if (sp)
         if (sscanf(sp, "%x", &i) == 1)
            part = i;
      free(res);
   }
   res = atlsys_1L(NULL, "grep -F 'CPU variant' /proc/cpuinfo", 0, 0);
   if (res)
   {
      char *sp;
      uint i;
      sp = strstr(res, "0x");
      if (sp)
         if (sscanf(sp, "%x", &i) == 1)
            *VAR = i;
      free(res);
   }
   res = atlsys_1L(NULL, "grep -F 'CPU implement' /proc/cpuinfo", 0, 0);
   if (res)
   {
      char *sp;
      uint i;
      sp = strstr(res, "0x");
      if (sp)
         if (sscanf(sp, "%x", &i) == 1)
            *IMPL = i;
      free(res);
   }

   res = atlsys_1L(NULL, "grep -F 'CPU architect' /proc/cpuinfo", 0, 0);
   if (res)
   {
      char *sp;
      uint i;
      sp = strstr(res, ": ");
      if (sp)
         if (sscanf(sp+1, "%u", &i) == 1)
            *ARCH = i;
      free(res);
   }
   res = atlsys_1L(NULL, "grep -F 'CPU revision' /proc/cpuinfo", 0, 0);
   if (res)
   {
      char *sp;
      uint i;
      sp = strstr(res, ": ");
      if (sp)
         if (sscanf(sp+1, "%u", &i) == 1)
            *REV = i;
      free(res);
   }
   return(part);
}
enum MACHTYPE ProbeArch()
{
   enum ARCHFAM fam;
   enum MACHTYPE mach=MACHOther;
   int ierr, i;
   char *res;

   fam = ProbeArchFam(NULL);
   switch(fam)
   {
   case AFPPC:
      res = atlsys_1L(NULL, "cat /proc/cpuinfo | grep -F cpu", 0, 0);
      if (res)
      {
#if 0
         if (strstr(res, "604e")) mach = PPC604e;
         else if (strstr(res, "604")) mach = PPC604;
         else
#endif
         if (strstr(res, "G4")) mach = PPCG4;
         else if (strstr(res, "7400")) mach = PPCG4;
         else if (strstr(res, "7410")) mach = PPCG4;
         else if (strstr(res, "7447")) mach = PPCG4;
         else if (strstr(res, "7455")) mach = PPCG4;
         else if (strstr(res, "PPC970FX")) mach = PPCG5;
         else if (strstr(res, "PPC970MP")) mach = PPCG5;
         else if (strstr(res, "POWER9")) mach = IbmPwr9;
         else if (strstr(res, "POWER8")) mach = IbmPwr8;
         else if (strstr(res, "POWER7")) mach = IbmPwr7;
         else if (strstr(res, "POWER6")) mach = IbmPwr6;
         else if (strstr(res, "POWER5")) mach = IbmPwr5;
         else if (strstr(res, "POWER4")) mach = IbmPwr4;
         else if (strstr(res, "e6500")) mach = Pwre6500;
         free(res);
      }
      break;
   case AFMIPS:
      res = atlsys_1L(NULL, "grep -F 'cpu model' /proc/cpuinfo", 0, 0);
      if (res)
      {
         if (res[0] != '\0')
         {
            if (strstr(res, "ICE9"))
               mach = MIPSICE9;
/*
 *          I have no access to what cpuinfo on Linux does for this procs,
 *          so this is a WAG as to what it would say
 */
            else if (strstr(res, "R10000") || strstr(res, "R12000") ||
                     strstr(res, "R12000") || strstr(res, "R14000"))
               mach = MIPSR1xK;
         }
         free(res);
      }
      break;
   case AFARM:
      res = atlsys_1L(NULL, "grep -F 'Processor' /proc/cpuinfo", 0, 0);
      if (!res)
         res = atlsys_1L(NULL, "grep -F cpu /proc/cpuinfo", 0, 0);
      if (res)
      {
         if (strstr(res, "ARMv7") || strstr(res,"v7l"))
         {
            free(res);
            res = atlsys_1L(NULL, "grep -F 'CPU part' /proc/cpuinfo", 0, 0);
            if (res)
            {
               char *sp;
               sp = strstr(res, "0x");
               if (sp)
               {
                  int i;
                  if (sscanf(sp, "%x", &i) == 1)
                  {
                     i &= 0xFF;
                     switch(i)
                     {
                     case 7:
                        mach = ARM7;
                        break;
                     case 9:
                        mach = ARM9;
                        break;
                     case 15:
                        mach = ARM15;
                        break;
                     case 17:
                        mach = ARM17;
                        break;
                     default:;
                     }
                  }
               }
            }
         }
         else if (strstr(res, "AArch64") || strstr(res, "aarch64"))
         {
            free(res);
            res = atlsys_1L(NULL, "grep -F 'Hardware' /proc/cpuinfo", 0, 0);
            if (res && strstr(res, "X-Gene "))
               mach = ARM64xg;
            else
            {
               if(res)
                  free(res);
               res = atlsys_1L(NULL, "grep -F 'CPU part' /proc/cpuinfo", 0, 0);
               if (res)
               {
                  char *sp;
                  sp = strstr(res, "0x");
                  if (sp)
                  {
                     int i;
                     if (sscanf(sp, "%x", &i) == 1)
                     {
                        i &= 0xFF;
                        switch(i)
                        {
                        case 3:
                           mach = ARM64a53;
                           break;
                        case 7:
                           mach = ARM64a57;
                           break;
                        default:;
                        }
                     }
                  }
               }
            }
         }
         free(res);
      }
/*
 *    If on ARM, and haven't found arch yet, try a list of exact matches from
 *    /proc/cpuinfo that I've blindly matched from machines I have access to
 */
      if (mach == MACHOther)
      {
         uint part, impl, arch, var, rev;
         part = getArmInfo(&impl, &arch, &var, &rev);
/*
 *       Unfortunately, later linuxes seem to have no idea what my early
 *       ARM64xgene1 board is, as most info comes back 0.  So, say anytime
 *       implementer is 0x50 and arch=8 (only filled in values) that that is
 *       what we've got
 */
         if (impl == 0x50)
         {
            if (arch == 8 && !var && !part & !rev)
               mach = ARM64xg;
         }
/*
 *       My Cavium ThunderX system
 */
         else if (impl == 0x43 && part==0x0a1 && var == 0x1 && arch == 8)
            mach = ARM64thund;
         else if (impl == 0x41 && part==0xd07 && var == 0x1 && arch == 8)
            mach = ARM64A1100;
      }
      break;
   case AFIA64:
      res = atlsys_1L(NULL, "grep -F 'Itanium' /proc/cpuinfo", 0, 0);
      if (res && res[0] == '\0')
      {
         free(res);
         res = NULL;
      }
      if (!res)
         res = atlsys_1L(NULL, "grep -F \"model name\" /proc/cpuinfo", 0, 0);
      if (res)
      {
         if (res[0] != '\0')
         {
            if (strstr(res, "Itanium 2") || strstr(res, "McKinley"))
               mach = IA64Itan2;
            else if (strstr(res, "Itanium")) mach = IA64Itan;
         }
         free(res);
      }
      break;
   case AFX86:
      res = atlsys_1L(NULL, "grep -F 'model name' /proc/cpuinfo", 0, 0);
      if (res && res[0] == '\0')
      {
         free(res);
         res = NULL;
      }
      if (!res)
         res = atlsys_1L(NULL, "grep -F model /proc/cpuinfo", 0, 0);
      if (res && res[0] == '\0')
      {
         free(res);
         res = NULL;
      }
      if (res)
      {
         if (strstr(res, "Pentium"))
         { /* Pentium of some flavor */
            if (strstr(res, " III ")) mach = IntPIII;
            else if (strstr(res, " II ")) mach = IntPII;
            else if (strstr(res, "Pro")) mach = IntPPRO;
            else if (strstr(res, "MMX")) mach = IntP5MMX;
            else if (strstr(res, " 4 "))
            {
               char *rs2;
               rs2 = atlsys_1L(NULL,
                      "grep -F 'model' /proc/cpuinfo | grep -F -v 'name'", 0, 0);
               if (rs2)
               {
                  i = GetLastInt(res);
                  if (i < 3) mach = IntP4;
                  else if (i == 3) mach = IntP4E;
                  free(rs2);
               }
            }
         }
         else if (strstr(res, "Atom"))
            mach = IntAtom;
         else if (strstr(res, "Core"))
         {
            if (strstr(res, "i7"))
            {
               if (strstr(res, "4770") || strstr(res, "4765") ||
                   strstr(res, "8700")) /* 8700: coffee-lake */
                  mach = IntCorei3;
               else if (strstr(res, "2600") || strstr(res, "3770"))
                  mach = IntCorei2;
               else if (strstr(res, "6700"))
                  mach = IntCorei4;
               else
                  mach = IntCorei1;
            }
            if (strstr(res, "i5"))
            {
               if (strstr(res, "4670") || strstr(res, "4570") ||
                   strstr(res, "4430") || strstr(res, "8600") ||
                   strstr(res, "8400")) /* 8600,8400: coffee-lake */
                  mach = IntCorei3;
               else if (strstr(res, "i5-2500") || strstr(res, "i5-2400") ||
	           strstr(res, "i5-2390") || strstr(res, "i5-2300"))
                  mach = IntCorei2;
               else
                  mach = IntCorei1;
            }
         }
         else if (strstr(res, "Xeon")) /* dreaded Xeon-is-anything */
         {
            if (strstr(res, "E5420")) mach = IntCore2;
            else if (strstr(res, "X7560")) mach = IntCorei1;
         }
         else if (strstr(res, "A8-3850"))
            mach = AmdLlano;
         else if (strstr(res, "Efficeon")) mach = TMEff;
         else if (strstr(res, "Athlon HX")) mach = AmdHammer;
         else if (strstr(res, "Opteron") || strstr(res, "Hammer") ||
                  strstr(res, "Athlon(tm) 64"))
            mach = AmdHammer;
         else if (strstr(res, "Athlon")) mach = AmdAthlon;
         else if (strstr(res, "AMD-K7")) mach = AmdAthlon;
         free(res);
      }
      break;
/*
 *    Add these back if we get machine access and can test
 */
   case AFSPARC:  /* don't know here anymore */
      res = atlsys_1L(NULL, "grep -F cpu /proc/cpuinfo", 0, 0);
      if (res)
      {
         if (strstr(res, "UltraSparc T2")) mach = SunUST2;
         else if (strstr(res, "UltraSparc T1")) mach = SunUST1;
         else if (strstr(res, "UltraSparc IV")) mach = SunUSIV;
         else if (strstr(res, "UltraSparc III")) mach = SunUSIII;
         else if (strstr(res, "UltraSparc II")) mach = SunUSII;
         else if (strstr(res, "UltraSparc I")) mach = SunUSI;
         else if (strstr(res, "UltraSparc")) mach = SunUSX;
         free(res);
      }
      break;
   case AFALPHA:
      #if 0
      res[0] = '\0';
      ierr = CmndOneLine(NULL, "grep -F 'model name' /proc/cpuinfo", res);
      if (ierr || res[0] == '\0')
         ierr = CmndOneLine(NULL, "grep -F model /proc/cpuinfo", res);
      if (!ierr && res[0] != '\0')
      {
         if (strstr(res, "EV5")) mach = Dec21164;
         else if (strstr(res, "EV4")) mach = Dec21064;
         else if (strstr(res, "EV6")) mach = Dec21264;
      }
      #endif
      break;
   case AFS390:
      res = atlsys_1L(NULL, "cat /proc/cpuinfo | grep -F \"processor \"", 0, 0);
      if (res)
      {
         if (strstr(res, "2094") || strstr(res, "2096")) mach = IbmZ9;
         else if (strstr(res, "2097") || strstr(res, "2098")) mach = IbmZ10;
         else if (strstr(res, "2817") || strstr(res, "2818")) mach = IbmZ196;
         else if (strstr(res, "2827") || strstr(res, "2828")) mach = IbmZ12;
         else if (strstr(res, "2964") || strstr(res, "2965")) mach = IbmZ13;
         else mach = IbmZ13;  /* looks risky to me, but IBM folks did it */
         free(res);
      }
      break;
   default:
#if 0
      if (!CmndOneLine(NULL, "grep -F 'cpu family' /proc/cpuinfo", res))
         if (strstr(res, "PA-RISC 2.0")) mach = HPPA20;
#else
     ;
#endif
   }
   return(mach);
}
int ProbeNCPU()
{
   int ncpu = 0;
   char *res;

   #if 0
   if (mach == Dec21264 || mach == Dec21164 || mach == Dec21064)
   {
      if ( !CmndOneLine(NULL, "grep -F 'cpus detected' /proc/cpuinfo", res) )
         ncpu = GetLastInt(res);
   }
   #endif
   if (!ncpu)
   {
      FILE *fpres;
      int len=0;
      fpres = atlsys(NULL, "grep '^processor' /proc/cpuinfo", 0, 0);
      if (fpres)
      {
         for (ncpu=0; res = ATL_fgets(res, &len, fpres); ncpu++);
         fclose(fpres);
      }
   }
   if (!ncpu)
   {
      res = atlsys_1L(NULL, "grep 'ncpus active' /proc/cpuinfo", 0, 0);
      if (res)
      {
         ncpu = GetFirstInt(res);
         free(res);
      }
   }
   return(ncpu);
}

int ProbePointerBits(int *sure)
{
   int i;
   char *uname;
   char *cmnd, *res;

   *sure = 0;
/*
 * This probe should be running on backend; if its ptr length is 8, we've
 * definitely got a 64 bit machine
 * NOTE: getting 4 could be a result of compiler flags on a 64-bit arch,
 * so reverse is not dispositive
 */
   if (sizeof(void*) == 8)
   {
      *sure = 1;
      return(64);
   }

/*
 * Note this is a weak probe, archinfo_x86 much better . . .
 */
   uname = FindUname(NULL);
   i = strlen(uname) + 4;
   cmnd = malloc(sizeof(char)*i);
   assert(cmnd);
   sprintf(cmnd, "%s -a", uname);

   res = atlsys_1L(NULL, cmnd, 0, 0);
   free(cmnd);
   if (res)
   {
/*
 *    If uname is a known 64-bit platform, we're sure we've got OS support
 *    for 64bits (may not have compiler support, but that's not our fault)
 */
      if (strstr(res, "x86_64") || strstr(res, "ppc64") || strstr(res, "ia64"))
      {
         *sure = 1;
         free(res);
         return(64);
      }
      free(res);
   }
   return(32);
}

int ProbeMhz()
{
   int mhz=0;
   char *res;
   res = atlsys_1L(NULL, "grep -F 'cpu MHz' /proc/cpuinfo", 0, 0);
   if (res)
   {
      mhz = GetFirstDouble(res) + 0.5;
      free(res);
   }
   if (!mhz)
   {
      res = atlsys_1L(NULL, "cat /proc/cpuinfo | grep -F clock | grep -F MHz",0,0);
      if (res)
      {
         mhz = GetLastLongWithRound(res);
         free(res);
      }
   }
/*
 * Try Linux/SPARC lookup
 */
   if (!mhz)
   {
      res = atlsys_1L(NULL, "grep -F 'ClkTck' /proc/cpuinfo", 0, 0);
      if (res)
      {
         mhz = GetLastHex(res) / 1000000;
         free(res);
      }
   }
/*
 * If nothing, try cpufreq-info of various flavors
 */
   if (!mhz)
   {
      res = atlsys_1L(NULL, "cpufreq-info -f -m", 0, 0);
      if (res)
      {
         char *sp;
         sp = strstr(res, "Hz");
         if (sp)
         {
           double dmhz;
            dmhz = GetFirstDouble(res);
            if (dmhz != 0.0)
            {
               switch(sp[-1])
               {
               case 'G':
                  mhz = dmhz / 1000;
                  break;
               case ' ':
                  mhz = dmhz / 1000000;
                  break;
               case 'K':
                  mhz = dmhz / 1000;
                  break;
               case 'M':
               default:
                  mhz = dmhz;
                  break;
               }
            }
         }
         free(res);
      }
   }
   if (!mhz)
   {
      res = atlsys_1L(NULL, "cpufreq-info -f", 0, 0);
      if (res)
      {
         mhz = GetLastInt(res);  /* assumes clock speed given in KHz */
         free(res);
         if (mhz)
            mhz = mhz / 1000;
      }
   }
/*
 * RCW: Patch supplied by IBM.  I have no idea why it is written this way.
 * S390 CPU speed is in bogomips rather than mhz
 */
   if (!mhz)
   {
      res = atlsys_1L(NULL, "cat /proc/cpuinfo | grep -F bogomips",0,0);
      if (res)
      {
         double result = GetFirstDouble(res);
         result = (result*100.0)/88.0;  /* RCW: WTF? */
         /*if (result > 5000) result = 5000;*/
         mhz = result;
         free(res);
      }
   }
   return(mhz);
}

int ProbeThrottle()
/*
 * RETURNS: 1 if cpu throttling is detected, 0 otherwise
 */
{
   int iret=0;
   int imax=0, imin=0, icur=0;
   char *res;

/*
 * If cpufreq directory doesn't exist, guess no throttling.  If
 * cpufreq exists, and cur Mhz < max, throttling is enabled,
 * throttling also enabled if governer is not "performance", and min freq
 * is less than max
 */
   res = atlsys_1L(NULL,
                   "cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq",
                   0, 0);
   if (res)
   {
      imax = GetFirstInt(res);
      free(res);
      res = atlsys_1L(NULL,
            "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_min_freq", 0, 0);
      assert(res);
      imin = GetFirstInt(res);
      free(res);
      res = atlsys_1L(NULL,
            "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq", 0, 0);
      assert(res);
      icur = GetFirstInt(res);
      free(res);
      res = atlsys_1L(NULL,
            "cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor", 0, 0);
      assert(res);
      if (icur < imax)
         iret = 1;
      else if (!strstr(res, "performance") && imin < imax)
         iret = 1;
      free(res);
   }
   return(iret);
}

int main(int nargs, char **args)
{
   int flags, CacheLevel, ncpu, mhz, bits, sure;
   enum MACHTYPE arch=MACHOther;

   flags = GetFlags(nargs, args, &CacheLevel);
   if (flags & Parch)
   {
      arch = ProbeArch();
      if (flags & Pverb)
         printf("Architecture detected as %s.\n", machnam[arch]);
      printf("MACHTYPE=%d\n", arch);
   }
   if (flags & Pncpu)
      printf("NCPU=%d\n", ProbeNCPU());
   if (flags & PMhz)
      printf("CPU MHZ=%d\n", ProbeMhz());
   if (flags & Pthrottle)
      printf("CPU THROTTLE=%d\n", ProbeThrottle());
   if (flags & P64)
   {
      bits = ProbePointerBits(&sure);
      printf("PTR BITS=%d, SURE=%d\n", bits, sure);
   }

/*
 * Here for future, presently unsupported
 */
   if (flags & Pncache)
      printf("NCACHES=0\n");
   if (flags & PCacheSize)
      printf("%d Cache size (kb) = 0\n", CacheLevel);
   return(0);
}
