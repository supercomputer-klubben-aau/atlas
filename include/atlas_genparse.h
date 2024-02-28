#ifndef ATLAS_GENPARSE_H
   #define ATLAS_GENPARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#define NASMD 11
enum ASMDIA
   {ASM_None=0, gas_x86_32, gas_x86_64, gas_sparc, gas_ppc, gas_parisc,
    gas_mips, gas_arm, gas_arm64, gas_wow64, gas_s390};
static char *ASMNAM[NASMD] =
   {"",     "GAS_x8632", "GAS_x8664", "GAS_SPARC", "GAS_PPC", "GAS_PARISC",
    "GAS_MIPS", "GAS_ARM", "GAS_ARM64", "GAS_WOW64", "GAS_S390"};
/*
 * Basic data structure for forming queues with some minimal info
 */
typedef struct SIDNode ATL_sidnode_t;
struct SIDNode  /* holds string, integer, and double */
{
   double d;
   char *str;
   int i;
   ATL_sidnode_t *next;
};

#define SET_FLAG(bits_, MMR_, val_) \
{\
   if (val_) (bits_) |= (1<<(MMR_)); \
   else (bits_) &= ~(1<<(MMR_)); \
}
#define FLAG_IS_SET(field_, bit_) ( ((field_) & (1<<(bit_))) != 0 )

/* procedure 1 */
int Sys2File(char *cmnd, char *file)
/*
 * Executes cmnd, and redirects stdout to file
 * RETURNS: 0 for success, else error code from system call
 */
{
   int lF, lC, i;
   char *sp, *sys;
   assert(cmnd);
   if (!file)
      file = "/dev/null";
   else if (!strcmp(file, "stdout"))
      file = NULL;
   lC = strlen(cmnd);
   if (file)
      lF = strlen(file);
   else
      lF = 0;
   sys = malloc((lC+lF+10)*sizeof(char));
   assert(sys);
   for (i=0; i < lC; i++)
      sys[i] = cmnd[i];
   if (file)
   {
      sp = sys + lC;
      *sp = ' ';
      sp[1] = '>';
      sp[2] = ' ';
      sp += 3;
      for (i=0; i < lF; i++)
         sp[i] = file[i];
      sp[lF] = ' ';
      sp[lF+1] = '2';
      sp[lF+2] = '>';
      sp[lF+3] = '&';
      sp[lF+4] = '1';
      sp[lF+5] = '\0';
   }
   else
      sys[lC] = '\0';
   i = system(sys);
   if (i)
      fprintf(stderr, "\nERROR %d for '%s'!\n\n", i, sys);
   free(sys);
   return(i);
}

/* procedure 2 : allocates ATL_sidnode_t */
static ATL_sidnode_t *ATL_NewSIDNode(void)
{
   ATL_sidnode_t *sp;
   sp = calloc(1, sizeof(ATL_sidnode_t));
   assert(sp);
   return(sp);
}

/* procedure 3 : allocates ATL_sidnode_t */
static ATL_sidnode_t *ATL_FreeSIDNode(ATL_sidnode_t *die)
{
   ATL_sidnode_t *sp=NULL;
   if (die)
   {
      sp = die->next;
      if (die->str)
         free(die->str);
      free(die);
   }
   return(sp);
}

/* procedure 4 */
static int GetL1CacheElts(char pre)
{
   FILE *L1f;
   int L1Size, i;

   L1f = fopen("res/L1CacheSize", "r");
   if (!L1f)
   {
      assert(system("make res/L1CacheSize\n") == 0);
      L1f = fopen("res/L1CacheSize", "r");
   }
   assert(L1f);
   assert(fscanf(L1f, "%d", &L1Size) == 1);
   fclose(L1f);
   assert(L1Size > 0);
   if (pre == 'c' || pre == 'd')
      i = 1024/8;
   else if (pre == 's')
      i = 1024/4;
   else if (pre == 'z')
      i = 1024/16;
   else
      assert(0);
   return(i*L1Size);
}

/* procedure 5 */
static char *ExtendString(char *str, int len)
/*
 * Given already-allocated str, allocate a new string of total length len,
 * copy str to it (strlen(str)<= len)
 */
{
   char *sp;
   sp = realloc(str, len);
   assert(sp);
   return(sp);
}

/* procedure 6 */
static char *GetStrForSprintf
(
   char *form,  /* format string that will be passed to printf */
   int extra,   /* extra chars over format length to allocate */
   char *old    /* original ptr to pass to realloc */
)
{
   old = realloc(old, strlen(form)+extra);
   assert(old);
   return(old);
}

/* procedure 7 */
static int NumDecDigits(int num)
/*
 * RETURNS: number of decimal digits required to hold num, wt sign of neg #s
 *          counted as a digit
 */
{
   int bound, nd;
   if (num < 0)
   {
      nd = 2;
      num = -num;
   }
   else
      nd = 1;
   for (bound=9; num > bound; nd++)
      bound = bound*10 + 9;
   return(nd);
}

/* procedure 8 */
static int NumHexDigits(unsigned int num)
{
   int dig;

   for (dig=1; num; num >>= 4)
      dig++;
   return(dig);
}

/* procedure 9 */
static char *DupString(char *str)
{
   int i,n;
   char *s;

   if (!str)
      return(NULL);
   n = strlen(str)+1;
   s = malloc(sizeof(char)*n);
   assert(s);
   for (i=0; i < n; i++)
      s[i] = str[i];
   return(s);
}

/* procedure 10 */
static char *NewMergedString(char *st0, char *st1)
/*
 * RETURNS: new string with st1 concatonated to st0
 */
{
   int i, n0, n1;
   char *s;

   if (!st0)
   {
      if (!st1)
         return(NULL);
      else
         return(DupString(st1));
   }
   n0 = strlen(st0);
   n1 = strlen(st1) + 1;
   s = malloc(sizeof(char)*(n0+n1));
   assert(s);
   for (i=0; i < n0; i++)
      s[i] = st0[i];
   s += n0;
   for (i=0; i < n1; i++)
      s[i] = st1[i];
   return(s-n0);
}

/* procedure 11 */
static char *GetSingleQuoteString(char *str)
{
   char *sp;
   int i, n;

   assert(str[0] == '\'');
   for (i=1; str[i] && str[i] != '\''; i++);
   assert(str[i]);
   sp = malloc(i*sizeof(char));
   for (n=i,i=1; i < n; i++)
      sp[i-1] = str[i];
   sp[n-1] = '\0';
   return(sp);
}

/* procedure 12 */
static FILE *OpenFileWithPath(char *path, char *file, char *perm)
{
   char *fn;
   FILE *fp;
   assert(file && perm);
   if (path)
   {
      int plen, flen;
      char *fn;

      plen = strlen(path);
      flen = strlen(file);
      fn = malloc(plen+flen+2);
      assert(fn);
      memcpy(fn, path, plen);
      fn[plen] = '/';
      memcpy(fn+plen+1, file, flen);
      fn[plen+flen+1] = '\0';
      fp = fopen(fn, perm);
      free(fn);
   }
   else
      fp = fopen(file, perm);

   return(fp);
}

/* procedure 13 */
static int asmNames2bitfield(char *str)
/*
 * Takes str containing an assembly name list.  The list is ended by the first
 * white space or end of string.  List items are seperated by ',', and there
 * can be no whitespace in list.
 * RETURNS: bitfield with bits set corresponding to assemblies, 0 on error.
 */
{
   char asmname[64];
   int i, KeepOn, bits=0;

   do
   {
      for (i=0; !isspace(str[i]) && str[i] != ',' && str[i] && i < 64; i++)
         asmname[i] = str[i];
      asmname[i] = '\0';
      KeepOn = str[i] == ',';
      str += i+1;
      if (i >= 64)
         return(0);  /* no asm name > 63 in length */
      for (i=0; i < NASMD; i++)
      {
         if (!strcmp(ASMNAM[i], asmname))
         {
            bits |= (1<<i);
            break;
         }
      }
   }
   while(KeepOn);
   return(bits);
}

/* procedure 14 */
static int GetDoubleArr(char *str, int N, double *d)
/*
 * Reads in a list with form "%le,%le...,%le"; N-length d recieves doubles.
 * RETURNS: the number of doubles found, or N, whichever is less
 */
{
   int i=1;
   assert(sscanf(str, "%le", d) == 1);
   while (i < N)
   {
      str = strstr(str, ",");
      if (!str)
         break;
      str++;
      assert(sscanf(str, "%le", d+i) == 1);
      i++;
   }
   return(i);
}

/* procedure 15 */
static char *GetLongerString(char *shrt, int newlen)
/*
 * Allocates new string of size newlen, copies shrt into it, and frees shrt.
 */
{
   char *sp;

   sp = malloc(sizeof(char)*newlen);
   assert(sp);
   if (shrt)
   {
      strcpy(sp, shrt);
      free(shrt);
   }
   else if (newlen >= 0)
      sp[0] = '\0';
   return(sp);
}

/* procedure 16 */
static char *GetOneLine(FILE *fpin)
/*
 * RETURNS: string of one line from stream fpin,  NULL if stream exhausted.
 */
{
   const int inc=256;
   static int len=0;
   static char *ln, *sp;
   int i, j, KeepOn;

   if (!len)
   {
      ln = malloc(inc*sizeof(char));
      assert(ln);
      len = inc;
   }
   if (!fgets(ln, len, fpin))
      return(NULL);

   for (i=0; ln[i]; i++);  /* find end of string */
   if (!i) return(ln);
   while (ln[i-1] != '\n')    /* if last char not \n, read rest of line */
   {
      len += inc;
      ln = GetLongerString(ln, len);
      if (!fgets(ln+i, inc, fpin))
         return(ln);
       for (; ln[i]; i++);  /* find end of string */
   }
   return(ln);
}

/* procedure 17 */
static char *GetJoinedLines(FILE *fpin)
/*
 * Gets lines from file fpin; if last non-whitespace char is '\', joins lines
 * RETURNS: line from file including joining, NULL if fpin exhausted
 */
{
   char *ln, *sp;
   static char *join=NULL;
   static int jlen=0;
   int i, j, k;

   sp = ln = GetOneLine(fpin);
   if (!sp)
      return(NULL);
   j = 0;   /* current length of join string */
   if (ln)
   {
      for (i=0; ln[i]; i++);  /* find end of string */
      if (!i) return(NULL);
      for (i--; isspace(ln[i]) && i > 0; i--);  /* find last non-wspace char */
      while (ln[i] == '\\')
      {
         if (jlen < j+i+3)
         {
            jlen = j+i+i+3;
            join = GetLongerString(join, jlen);
         }
         for (k=0; k < i; k++)
            join[j+k] = ln[k];
         j += k;
         join[j++] = ' ';
         join[j] = '\0';
         ln = GetOneLine(fpin);   /* get new line that should be joined */
         assert(ln);              /* can't end file with continue */
         for (i=0; ln[i]; i++);   /* find end of new line */
         for (i--; isspace(ln[i]) && i > 0; i--); /* find last non-wspc char */
         sp = join;
      }
      if (sp == join)
      {
         if (jlen < j+i+3)
         {
            jlen = j+i+i+3;
            join = GetLongerString(join, jlen);
         }
         for (k=0; k <= i; k++)
            join[j+k] = ln[k];
         j += k;
         join[j] = '\n';
         join[j+1] = '\0';
         sp = join;
      }
   }
   return(sp);
}

/* procedure 18 */
static void ATL_isortG2L(int N, int *X)
/*
 * Sort X from greatest-to-least via selection sort (O(N^2))
 */
{
   int i;
   if (N < 2)
      return;
   for (i=0; i < N-1; i++)
   {
      int imax=X[i], idx=i, j;
      for (j=i+1; j < N; j++)
      {
         const int xj=X[j];
         if (xj > imax)
         {
            imax = xj;
            idx = j;
         }
      }
      if (idx != i)
      {
         X[idx] = X[i];
         X[i] = imax;
      }
   }
}

/* procedure 19 */
static char *GetGoodGcc()
/*
 * Gets gcc path and name along with mandatory flags (-g/-m64/-pg,etc) by
 * querying Make.inc setting
 */
{
   static char gcc[2048];
   static int INIT=0;
   if (!INIT)
   {
      FILE *fpin;
      assert(system("make res/goodgcc.txt > /dev/null 2>&1") == 0);
      fpin = fopen("res/goodgcc.txt", "r");
      assert(fpin);
      assert(fscanf(fpin, "'%[^\']", gcc) == 1);
      fclose(fpin);
   }
   return(gcc);
}

/* procedure 20 */
static char *GetKCFlags(char pre)
/*
 * Gets flags being used for <pre>KCFLAGS
 */
{
   char ln[4096];
   FILE *fpin;
   int i;

   if (pre == 'z')
      pre = 'd';
   else if (pre == 'c')
      pre = 's';
   i = system("rm -f res/kcflags.txt");
   sprintf(ln, "grep \"%cKCFLAGS = \" Make.inc | sed s/%cKCFLAGS\\ =\\ // > res/kcflags.txt", toupper(pre), toupper(pre));
   assert(system(ln) == 0);
   fpin = fopen("res/kcflags.txt", "r");
   assert(fpin);
   assert(fgets(ln, 4096, fpin) != NULL);
   fclose(fpin);
/*
 * Get rid of trailing and leading whitespaces
 */
   for (i=0; ln[i]; i++);
   for (i--; isspace(ln[i]); i--);
   ln[i+1] = '\0';
   for (i=0; isspace(ln[i]); i++);
   return(DupString(ln+i));
}

/* procedure 21 */
static int *GF_GetIntList1(int ival)
/*
 * returns integer array with iarr[0] = 1, iarr[1] = ival
 */
{
   int *iarr;
   iarr = malloc(2*sizeof(int));
   assert(iarr);
   iarr[0] = 1;
   iarr[1] = ival;
   return(iarr);
}

/* procedure 22 */
static int *GF_GetIntList2(int ival1, int ival2)
/*
 * returns integer array with iarr[0] = 1, iarr[1] = ival1, ival[2] = ival2
 */
{
   int *iarr;
   iarr = malloc(3*sizeof(int));
   assert(iarr);
   iarr[0] = 2;
   iarr[1] = ival1;
   iarr[2] = ival2;
   return(iarr);
}

/* procedure 23 */
static int *GF_DupIntList(int *L)
/*
 * dups a list of integers L, whose data length is given by L[0];
 * list is this length+1, since 0'th location gets data length.
 */
{
   int *ip, n;
   if (!L)
      return(NULL);
   n = L[0] + 1;
   ip = malloc(n*sizeof(int));
   assert(ip);
   memcpy(ip, L, n*sizeof(int));
   return(ip);
}
#ifdef ATL_GETFLAGS
/* procedure 24 */
static double *GF_GetNDoubleArgs(int nargs, char **args, int i, int n)
/*
 * Reads in n doubles from commandline args to produce n-len double array.
 */
{
   int k;
   double *darr;
   void PrintUsage(char*, int, char*);

   if (n < 1)
      return(NULL);
   darr = malloc(sizeof(double)*n);
   assert(darr);

   for (k=0; k < n; k++)
   {
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      darr[k] = atof(args[i]);
   }
   return(darr);
}

/* procedure 25 */
static double *GF_GetDoubleList(int nargs, char **args, int i, int nmul)
/*
 * Gets a list of doubles, whose length is given by atoi(args[i])*nmul
 * list is this length+1, since 0'th location gets atoi(args[i])
 */
{
   int n, k;
   double *darr;
   void PrintUsage(char*, int, char*);

   if (++i >= nargs)
      PrintUsage(args[0], i, NULL);
   n = atoi(args[i]) * nmul;
   assert(n > 0);
   darr = malloc(sizeof(double)*(n+1));
   assert(darr);

   darr[0] = n / nmul;
   for (k=0; k < n; k++)
   {
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      darr[k+1] = atof(args[i]);
   }
   return(darr);
}

/* procedure 26 */
static int *GF_GetIntList(int nargs, char **args, int i, int nmul)
/*
 * Gets a list of integers, whose length is given by atoi(args[i])*nmul
 * list is this length+1, since 0'th location gets atoi(args[i])
 */
{
   int n, *iarr, k;
   void PrintUsage(char*, int, char*);

   if (++i >= nargs)
      PrintUsage(args[0], i, NULL);
   n = atoi(args[i]) * nmul;
   assert(n > 0);
   iarr = malloc(sizeof(int)*(n+1));
   assert(iarr);

   iarr[0] = n / nmul;
   for (k=0; k < n; k++)
   {
      if (++i >= nargs)
         PrintUsage(args[0], i, NULL);
      iarr[k+1] = atoi(args[i]);
   }
   return(iarr);
}
#endif

/* procedure 27 */
static int *GF_IntRange2IntList(int N0, int NN, int incN)
{
   int i, n;
   int *iarr;

   for (i=N0, n=0; i <= NN; i += incN) n++;
   iarr = malloc(sizeof(int)*(n+1));
   assert(iarr);
   iarr[0] = n;
   for (i=N0, n=1 ; i <= NN; i += incN, n++)
      iarr[n] = i;
   return(iarr);
}
/* procedure 28 */
static int GetNativeVLEN(char pre)
{
#ifdef ATL_AVXZ
   if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
      return(16);
   return(8);
#elif defined(ATL_AVX) || defined(ATL_AVXMAC)
   if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
      return(8);
   return(4);
#elif defined(ATL_GAS_ARM64)
      if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
         return(4);
      return(2);
#elif defined(ATL_GAS_ARM)
   #ifndef ATL_NONIEEE
      return(1);
   #else
      if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
         return(4);
      return(1);
   #endif
#elif defined(ATL_SSE2) || defined(ATL_VECARM1) || defined(ATL_VSX)
   if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
      return(4);
   return(2);
#elif defined(ATL_SSE1)
   if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
      return(4);
   return(1);
#else
   return(0);   /* 0 means: I don't know */
#endif
}

/* procedure 29 */
static int ATL_numBitsNeeded(int N)
/*
 * RETURNS: min number of bits to store value N
 * NOTE: assumes negative numbers are flags meaning not stored!
 */
{
   int i;
   if (N <= 0)
      return(0);
   for (i=0; (1<<i) <= N; i++);
   return(i);
}

/* procedure 30 */
static int GetNumVecRegs(char pre)
{
#ifdef ATL_AVXZ
   #ifdef ATL_USE64BITS
      return(32);
   #else
      return(8);
   #endif
#elif defined(ATL_SSE1)
   #ifdef ATL_USE64BITS
      return(16);
   #else
      return(8);
   #endif
#elif defined(ATL_VECARM1) || defined(ATL_FPV3D32MAC)
   return(32);
#elif defined(ATL_FPV3D16MAC)
   if (pre == 's' || pre == 'c' || pre == 'S' || pre == 'C')
      return(32);
    return(16);
#elif defined(ATL_VSX)
   return(64);
#else
   return(0);   /* 0 means: I don't know */
#endif
}

/*
 * The following are functions that allow us to operate on generic lists
 * of any structure.  The idea is that you pass the byte offset where
 * the next pointer is found, as well as the data item you want to access.
 * Since we are loading raw memory, note that your type must match precisely.
 * Eg., GetIntAtOff is only used with int (signed or unsigned), not any
 * intergral (it will fail for short or char).  Assuming you have a pointer
 * to the structure in question:
 *    GetOffset(&p->next, p);
 * can be used to calculate the offnxt these funcs require.
 */
/* procedure 31 */
static int GetOffset(void *v0, void *v1)
/*
 * Find (postive) gap between addresses.  Used inside structs, where dist
 * is short enough that size_t not required.
 */
{
   size_t i0=(size_t)v0, i1=(size_t)v1;
   if (v1 > v0)
      return(v1-v0);
   return(v0-v1);
}

/* procedure 32 */
static int GetIntAtOff(void *p, int off)
/*
 * RETURNS: derefence *(p+off) interpeted an integer address.
 */
{
   return(*((int*)(((char*)p)+off)));
}

/* procedure 33 */
static void *GetPtrAtOff(void *p, int off)
{
   return( *((void**)(((char*)p)+off)) );
}

/* procedure 34 */
static char *GetStrAtOff(void *p, int off)
{
   return( *((char**)(((char*)p)+off)) );
}

/* procedure 35 */
static double GetDoubleAtOff(void *mp, int off)
{
   return(*((double*)(((char*)mp)+off)));
}

/* procedure 36 */
static void SetPtrAtOff(void *vb, int off, void *p)
{
   *((void**)(((char*)vb)+off)) = p;
}

/* procedure 37 */
static int CountListEntries(void *vb, int nxtoff)
{
   int i=0;
   char *cp = vb;
   while (cp)
   {
      cp = GetStrAtOff(cp, nxtoff);
      i++;
   }
   return(i);
}

/* procedure 38 */
static int CountListMaskALL(void *vb, int msk, int nxtoff, int ioff)
/*
 * counts number of nodes where all bits set in msk are set in int at ioff
 */
{
   int i=0;
   char *cp = vb;
   while (cp)
   {
      int flag;
      flag = GetIntAtOff(cp, ioff);
      if ((flag&msk) == msk)
         i++;
      cp = GetStrAtOff(cp, nxtoff);
   }
   return(i);
}

/* procedure 39 */
static int GetIntMaxMinAtOff(void *vb, int nxtoff, int off, int *MAX, int *MIN)
/*
 * Find the max and min of any integer in vb+off in list, next ptr at nxtoff;
 * Works for any list with a next ptr, and any int stored as int native
 * (will not work for short/char/etc).
 * RETURNS: number of nodes searched.
 */
{
   char *cp=vb;
   register int max=0, min=0, n=0;

   if (vb)
   {
      min = max = GetIntAtOff(cp, off);
      n = 1;
      while( (cp = GetStrAtOff(cp, nxtoff)) )
      {
         register int i;
         i = GetIntAtOff(cp, off);
         max = (max >= i) ? max : i;
         min = (min <= i) ? min : i;
         n++;
      }
   }
   *MAX = max;
   *MIN = min;
   return(n);
}

/* procedure 40 */
void *iFindMinAtOff(void *vb, int nxtoff, int off)
/*
 * RETURNS: ptr to struct containing minimum integer at off
 */
{
   void *ret=NULL;
   if (vb)
   {
      char *cp=vb;
      int min;

      ret = vb;
      min = GetIntAtOff(cp, off);
      while( (cp = GetStrAtOff(cp, nxtoff)) )
      {
         int i;
         i = GetIntAtOff(cp, off);
         if (i < min)
         {
            ret = cp;
            min = i;
         }
      }
   }
   return(ret);
}

/* procedure 41 */
static int *GetIntsFromQ(void *vb, int nxtoff, int off)
/*
 * Given N nodes in vb, returns an N+1 int array, with 1st elt storing N,
 * and subsequent elts storing the integer found at off in the struct.
 * RETURNS: NULL for NULL vb, else above.
 */
{

   if (vb)
   {
      char *cp=vb;
      int *ip;
      int i, n;

      n = CountListEntries(vb, nxtoff);
      ip = malloc((n+1)*sizeof(int));
      assert(ip);
      ip[0] = n;
      for (i=0; i < n; i++)
      {
         ip[i+1] = GetIntAtOff(vb, off);
         vb = GetPtrAtOff(vb, nxtoff);
      }
      return(ip);
   }
   return(NULL);
}

/* procedure 42 */
void *ATL_FindLastNode(void *vp, int nxtoff)
{
   if (vp)
   {
      char *nxt;
      while ( (nxt = GetStrAtOff(vp, nxtoff)) )
         vp = nxt;
   }
   return(vp);
}

/* procedure 43 */
void *FindNodeWithMaskOR(void *p, int nxtoff, int bvoff, int mask)
/*
 * Assuming list with next ptr at nxtoff, and bitvec at bvoff, search
 * for any node with any bits present int mask set.
 * RETURNS: ptr to first node with bits matching any bit in mask, if none, NULL
 * NOTE: this function does a logical OR (any mask bit set) operation, not AND!
 */
{
   while (p)
   {
      int bv;
      bv = GetIntAtOff(p, bvoff);
      if (bv & mask)
         return(p);
      p = GetStrAtOff(p, nxtoff);
   }
   return(NULL);
}

/* procedure 44 */
void PrefixStrAllNodes(void *p, int nxtoff, int stoff, char *str)
/*
 * Prefixes all strings at stoff in queue starting at p if not already prefxd
 */
{
   unsigned int len;
   if (!str || !p)
      return;
   len = strlen(str);
   if (!len)
      return;
   while(p)
   {
      char *os, *ns;
      unsigned int olen;
      int DOPREF;

      os = GetStrAtOff(p, stoff);
      if (os)
      {
         olen = strlen(os);
         if (olen < len)
            DOPREF = 1;
         else
            DOPREF = strncmp(os, str, len);
         if (DOPREF)
         {
            ns = malloc(len+olen+1);
            assert(ns);
            strcpy(ns, str);
            strcpy(ns+len, os);
            free(os);
            SetPtrAtOff(p, stoff, ns);
         }
      }
      p = GetPtrAtOff(p, nxtoff);
   }
}

/* procedure 45 */
void *FindNodeWithIval(void *p, int nxtoff, int ivoff, int val)
{
   while (p)
   {
      int iv;
      iv = GetIntAtOff(p, ivoff);
      if (iv == val)
         return(p);
      p = GetStrAtOff(p, nxtoff);
   }
   return(NULL);
}

/* procedure 46 */
void *RemoveNodeFromList(void *vb, void *vp, int nxtoff)
/*
 * Finds struct ptr vp in vb, and removes it from singly-linked list
 * RETURNS: possibly changed vb (only in case vb == vp)
 */
{
   if (vb)
   {
      char *nxt, *prv=vb, *cur;
      nxt = GetStrAtOff(vb, nxtoff);
      if (vb == vp)
      {
         SetPtrAtOff(vb, nxtoff, NULL);
         return((void*)nxt);
      }
      cur = nxt;

      while (cur)
      {
         nxt = GetStrAtOff(cur, nxtoff);
         if (cur == vp)
         {
            SetPtrAtOff(prv, nxtoff, nxt);
            SetPtrAtOff(cur, nxtoff, NULL);
            return(vb);
         }
         prv = cur;
         cur = nxt;
      }
   }
   if (vb)
   {
      fprintf(stderr, "ERROR: cannot find %p in basefile %p!\n", vp, vb);
      exit(1);
   }
   return(NULL);
}

/* procedure 47 */
void *SortListByIval_G2L(void *vb, int nxtoff, int ivoff)
/*
 * RETURNS: base ptr of list bv sorted from greatest-to-least by int at ivoff
 * NOTE: original list is reordered, so vb must be replaced by return baseptr!
 */
{
   void *sb=NULL;
   while (vb)
   {
      void *minp;
      minp = iFindMinAtOff(vb, nxtoff, ivoff);
      vb = RemoveNodeFromList(vb, minp, nxtoff);
      SetPtrAtOff(minp, nxtoff, sb);
      sb = minp;
   }
   return(sb);
}
#endif /* end atlas_genparse.h guard */
