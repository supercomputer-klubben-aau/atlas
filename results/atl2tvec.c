#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

static char *SEPSTR = ",";
typedef struct wOrDs WORDS;
struct wOrDs
{
   char *word;
   WORDS *next;
};

WORDS *GetWord(char *wrd, int wlen)
/*
 * returns a WORD, copying wlen chars from wrd into malloced space
 */
{
   WORDS *wp;
   wp = malloc(sizeof(WORDS));
   assert(wp);
   wp->word = malloc((wlen+1)*sizeof(char));
   assert(wp->word);
   strncpy(wp->word, wrd, wlen);
   wp->word[wlen] = '\0';
   wp->next = NULL;
   return(wp);
}

WORDS *KillWord(WORDS *wp)
{
   WORDS *next=NULL;
   if (wp)
   {
      next = wp->next;
      if (wp->word) free(wp->word);
      free (wp);
   }
   return(next);
}

WORDS *KillAllWords(WORDS *wp0)
{
   while (wp0)
      wp0 = KillWord(wp0);
   return(NULL);
}

void PrintWordList(WORDS *wp)
{
   int i;
   if (!wp) fprintf(stderr, "NULL wordlist\n");
   for (i=0; wp; i++, wp = wp->next)
   {
      fprintf(stderr, "word %d: %s\n", i, wp->word);
   }
}

WORDS *GetWords(char *ln)
{
   int i;
   WORDS *wp0, *wp;
   wp = wp0 = GetWord(" ", 1);
   while (*ln)
   {
      while(isspace(*ln)) ln++;
      if (*ln)
      {
         for (i=0; ln[i] && !isspace(ln[i]); i++);
         wp->next = GetWord(ln, i);
         wp = wp->next;
         ln += i;
      }
   }
   wp->next = NULL;
   if (wp0 != wp) wp = wp0->next;
   else wp = NULL;
   KillWord(wp0);
   return(wp);
}

int WordIsFloat(char *wrd)
/*
 * Returns 0 if wrd is not a valid float, 1 if it is
 * Floats have form [+,-][#].[#][e/d<exp>]
 * where at least one [#] must appear
 */
{
   int GotDigit=0, GotDot=0, GotExp=0;
   int i;

   switch(*wrd)
   {
   case '0':
   case '1':
   case '2':
   case '3':
   case '4':
   case '5':
   case '6':
   case '7':
   case '8':
   case '9':
      GotDigit = 1;
      break;
   case '-':
   case '+':
      break;
   case '.':
      GotDot = 1;
      break;
   default:
      return(0);
   }
   for (i=1; wrd[i]; i++)
   {
      if (wrd[i] == '.') GotDot++;
      else if (isdigit(wrd[i])) GotDigit++;
      else if (wrd[i] != 'e' || wrd[i] == 'd') GotExp++;
      else if (wrd[i] == '-' || wrd[i] == '+')
      {
         if (wrd[i-1] != 'e' && wrd[i-1] != 'd') return(0);
      }
      else return(0);
   }
   return ( (GotDot == 1) && (GotDigit) && (GotExp==1 || GotExp==0) );
}

int WordIsInt(char *wrd)
/*
 * Returns 0 if wrd is not a valid integer.
 */
{
   int i;

   for (i=0; wrd[i]; i++)
      if (!isdigit(wrd[i]))
         return(0);
   return(1);
}
void UpCaseWords(WORDS *wb)
{
   int i;
   for (; wb; wb = wb->next)
   {
      if (wb->word)
         for (i=0; wb->word[i]; i++)
            wb->word[i] = toupper(wb->word[i]);
   }
}

typedef struct wOrDQ WORDQ;
struct wOrDQ
{
   WORDS *wlist;
   WORDQ *next;
};

static WORDQ *GetWordQ(WORDS *wrds, WORDQ *next)
{
   WORDQ *wp;
   wp = malloc(sizeof(WORDQ));
   assert(wp);
   wp->wlist = wrds;
   wp->next = next;
   return(wp);
}

static WORDQ *KillWordQ(WORDQ *killme)
{
   WORDQ *next;
   if (killme)
   {
      if (killme->wlist)
         KillAllWords(killme->wlist);
      next = killme->next;
      free(killme);
      return(next);
   }
   return(NULL);
}

static WORDQ *KillAllWordQ(WORDQ *killme)
{
   while(killme)
      killme = KillWordQ(killme);
   return(NULL);
}
WORDQ *GetWordQTrans(WORDQ *wqN)
{
  WORDQ *wqT, *qT, *qN;
  WORDS *wpN, *wpT;
  int M, N, i, j, J;

  if (!wqN)
     return(NULL);
/*
 * Find M rows of data with N outputs, stored one row per WORDQ
 */
  for (N=0, wpN=wqN->wlist; wpN; N++, wpN = wpN->next);
  for (M=0, qN = wqN; qN; M++, qN = qN->next);
/*
 * Allocate N WORDQ entries for transposed data
 */
   wqT = qT = GetWordQ(NULL, NULL);
   for (j=1; j < N; j++)
   {
      qT->next = GetWordQ(NULL, NULL);
      qT = qT->next;
   }
/*
 * Transpose data
 */
   for (qT=wqT, J=0; J < N; J++, qT = qT->next)
   {
      for (wpN=wqN->wlist, j=0; j < J; j++, wpN = wpN->next);
      qT->wlist = GetWord(wpN->word, strlen(wpN->word));
      wpT = qT->wlist;
      for (qN=wqN->next, i=1; i < M; qN = qN->next, i++)
      {
         for (wpN=qN->wlist, j=0; j < J; j++, wpN = wpN->next);
         wpT->next = GetWord(wpN->word, strlen(wpN->word));
         wpT = wpT->next;
      }
   }
   return(wqT);
}

WORDQ *ReadTimerFile(FILE* fpin)
/*
 * This routine scans an ATLAS timer output file and returns a queue of
 * word lists, where each entry in the WORDQ is a row, and the entries
 * of WORDQ->wlist are column entries
 * Upon return, the first entry in WORDQ should be the column headings
 * (labels), and subsequent entries are values
 */
{
   char ln0[1024], ln1[1024], *nxtln, *prvln, *ln;
   WORDQ *wq, *wbase;
   WORDS *wrds;
/*
 * Search for line of === indicating previous line are headers
 */
   if (!fgets(ln0, 1024, fpin))
      return(NULL);
   prvln = ln0;
   nxtln = ln1;
   while (fgets(nxtln, 1024, fpin))
   {
      if (strstr(nxtln, "====")) break;
      ln = nxtln;
      nxtln = prvln;
      prvln = ln;
   }
   assert(strstr(nxtln, "===="));
   wbase = wq = GetWordQ(GetWords(prvln), NULL);
   assert(wbase->wlist);
   ln = ln0;
   while (fgets(ln, 1024, fpin))
   {
      wrds = GetWords(ln);
      if (wrds)
      {
         wq->next = GetWordQ(wrds, NULL);
         wq = wq->next;
      }
   }
   return(wbase);
}


void PrintWordsWithSep(FILE *fpout, WORDS *wrds, char *suff, char *sep)
{
   if (!wrds)
      return;
   fprintf(fpout, "%s", wrds->word ? wrds->word : "ERROR");
   if (suff)
      fprintf(fpout, "%s", suff);
   for (wrds = wrds->next; wrds; wrds = wrds->next)
      fprintf(fpout, "%s%s", sep, wrds->word ? wrds->word : "ERROR");
   fprintf(fpout, "\n");
}

void PrintWordLine(FILE *fpout, WORDS *wrds)
{
   PrintWordsWithSep(fpout, wrds, NULL, SEPSTR);
}

void PrintWordQ(FILE *fpout, char *name, WORDQ *wq)
{
   fprintf(fpout, "*****%s*****\n", name);
   while (wq)
   {
      PrintWordsWithSep(fpout, wq->wlist, NULL, SEPSTR);
      wq = wq->next;
   }
}

WORDS *GetColByNum(WORDQ *wqb, int I)
/*
 * Returns WORDS list of Ith column in row-major struct wqb
 */
{
   WORDQ *wq;
   WORDS *wp, *cb=NULL, *cp;
   int i;

   if (!wqb)
      return(NULL);
   for (wq=wqb; wq; wq = wq->next)
   {
      for (wp=wq->wlist, i=0; i < I; i++, wp = wp->next);
      if (cb)
      {
         cp->next = GetWord(wp->word, strlen(wp->word));
         cp = cp->next;
      }
      else
         cp = cb = GetWord(wp->word, strlen(wp->word));
   }
   return(cb);
}

WORDS *GetColByName(WORDQ *wqb, char *name)
/*
 * Finds column with matching name, assuming row-major wt labels in row 0
 */
{
   int i;
   WORDS *wp;

   if (!wqb)
      return(NULL);
   for (i=0, wp=wqb->wlist; wp && strcmp(name, wp->word); wp = wp->next, i++);
   if (!wp)
      return(NULL);
   return(GetColByNum(wqb, i));
}

double *Words2Doubles(WORDS *wpb, int *N)
/*
 * Translates words to double array, N is length of array
 */
{
   double *d;
   WORDS *wp;
   int n;

   for (n=0, wp=wpb; wp; n++, wp = wp->next);   /* count words */
   d = malloc(n*sizeof(double));
   for (n=0, wp=wpb; wp; n++, wp = wp->next)    /* convert to doubles */
      d[n] = atof(wp->word);
   *N = n;
   return(d);
}

WORDS *DupWordList(WORDS *wpb0)
/*
 * Duplicates wordlist wpb0
 */
{
   WORDS *wpb1, *wp1, *wp0;

   if (!wpb0)
      return(NULL);
   wpb1 = wp1 = GetWord(wpb0->word, strlen(wpb0->word));
   for (wp0=wpb0->next; wp0; wp0 = wp0->next)
   {
      wp1->next = GetWord(wp0->word, strlen(wp0->word));
      wp1 = wp1->next;
   }
   return(wpb1);
}

WORDQ *GetStridedQ(WORDQ *wqb, int stride)
/*
 * Allocates a new queue where inc elts of wqb are skipped between
 * each new queue
 */
{
   WORDQ *wqbS, *wqS, *wq;    /* S for strided */
   const int gap = stride-1;
   int i;

   if (!wqb)
      return(NULL);
   wqbS = wqS = GetWordQ(DupWordList(wqb->wlist), NULL);
   for (wq=wqb->next; wq; wq = wq->next)
   {
      for (i=0; i < gap && wq; i++)
         wq = wq->next;
      if (!wq)
         return(wqbS);
      wqS->next = GetWordQ(DupWordList(wq->wlist), NULL);
      wqS = wqS->next;
   }
   return(wqbS);
}

WORDQ *KillBadLines(WORDQ *wqb, int verb)
/*
 * Kills any WORDQ entry that doesn't have the same number of columns as wqb
 */
{
   int i, n;
   WORDQ *wq, *prev;
   WORDS *wrds;

   for (n=0,wrds=wqb->wlist; wrds; n++,wrds=wrds->next);  /* # of col labels */

   prev = wqb;
   wq = wqb->next;
   do
   {
      for (i=0,wrds=wq->wlist; wrds; i++,wrds=wrds->next);
      if (i != n)
      {
         if (verb)
         {
            fprintf(stderr, "killing anomolous line :");
            PrintWordLine(stderr, wq->wlist);
         }
         prev->next = wq = KillWordQ(wq);
      }
      else
      {
         prev = wq;
         wq = wq->next;
      }
   }
   while(wq);
   return(wqb);
}

void PrintDoublesWithSep(FILE *fpout, char *name, int N, double *d, char *sep)
{
   int i;

   fprintf(fpout, "%s", name);
   for (i=0; i < N; i++)
      fprintf(fpout, "%s%lf", sep, d[i]);
   fprintf(fpout, "\n");
}

void PrintDoublesWithSepRep(FILE *fpout, char *name, int n, int nrep,
                            double *d, char *sep)
/*
 * Print all of d, which is in nrep-major storage, length >= n*nrep
 */
{
   int i, j;

   for (i=0; i < nrep; i++)
   {
      fprintf(fpout, "%s_%d", name, i);
      for (j=0; j < n; j++)
         fprintf(fpout, "%s%lf", sep, d[i+j*nrep]);
      fprintf(fpout, "\n");
   }
}

/*
 * This meanDiff function provided by Tony Castaldo
 */
/******************************************************************************
 * meanDiff.c: Tests whether the means on two timing lists are different, and
 * returns the largest confidence interval of difference it can achieve. This
 * is accomplished with a standard T-Test.
 ******************************************************************************/
#include  <math.h>

#define MeanDiff_Rows 14                      /* Values per column */
#define MeanDiff_Cols 34                      /* Columns (DOF values) */
/*
 *-----------------------------------------------------------------------------
 * These are the critical values of the t-Distribution.
 * Consider this in column major order; the first 14 values are the first
 * column; the next 14 values are the next column, etc.
 * Each column corresponds to one of the 34 DOF values {1-30,40,60,120,infty}.
 * We make infty equivalent to 1024 or more.
 *-----------------------------------------------------------------------------
 */
double meanDiffT_Table[]={
  0.40,    0.30,    0.20,    0.15,    0.10,    0.05,    0.025,   0.02,
  0.015,   0.01,    0.0075,  0.005,   0.0025,  0.0005,

  0.325,   0.727,   1.376,   1.963,   3.078,   6.314,  12.706,  15.895, /* v=1*/
 21.205,  31.821,  42.434,  63.657, 127.322, 636.59,

  0.289,   0.617,   1.061,   1.386,   1.886,   2.920,   4.303,   4.849, /* v=2*/
  5.643,   6.965,   8.073,   9.925,  14.089,  31.598,

  0.277,   0.584,   0.978,   1.250,   1.638,   2.353,   3.182,   3.482, /* v=3*/
  3.896,   4.541,   5.047,   5.841,   7.453,  12.924,

  0.271,   0.569,   0.941,   1.190,   1.533,   2.132,   2.776,   2.999, /* v=4*/
  3.298,   3.747,   4.088,   4.604,   5.598,   8.610,

  0.267,   0.559,   0.920,   1.156,   1.476,   2.015,   2.571,   2.757, /* v=5*/
  3.003,   3.365,   3.634,   4.032,   4.773,   6.869,

  0.265,   0.553,   0.906,   1.134,   1.440,   1.943,   2.447,   2.612, /* v=6*/
  2.829,   3.143,   3.372,   3.707,   4.317,   5.959,

  0.263,   0.549,   0.896,   1.119,   1.415,   1.895,   2.365,   2.517, /* v=7*/
  2.715,   2.998,   3.203,   3.499,   4.029,   5.408,

  0.262,   0.546,   0.889,   1.108,   1.397,   1.860,   2.306,   2.449, /* v=8*/
  2.634,   2.896,   3.085,   3.355,   3.833,   5.041,

  0.261,   0.543,   0.883,   1.100,   1.383,   1.833,   2.262,   2.398, /* v=9*/
  2.574,   2.821,   2.998,   3.250,   3.690,   4.781,

  0.260,   0.542,   0.879,   1.093,   1.372,   1.812,   2.228,   2.359, /*v=10*/
  2.527,   2.764,   2.932,   3.169,   3.581,   4.587,

  0.260,   0.540,   0.876,   1.088,   1.363,   1.796,   2.201,   2.328, /*v=11*/
  2.491,   2.718,   2.879,   3.106,   3.497,   4.437,

  0.259,   0.539,   0.873,   1.083,   1.356,   1.782,   2.179,   2.303, /*v=12*/
  2.461,   2.681,   2.836,   3.055,   3.428,   4.318,

  0.259,   0.537,   0.870,   1.079,   1.350,   1.771,   2.160,   2.282, /*v=13*/
  2.436,   2.650,   2.801,   3.012,   3.372,   4.221,

  0.258,   0.537,   0.868,   1.076,   1.345,   1.761,   2.145,   2.264, /*v=14*/
  2.415,   2.624,   2.771,   2.977,   3.326,   4.140,

  0.258,   0.536,   0.866,   1.074,   1.341,   1.753,   2.131,   2.249, /*v=15*/
  2.397,   2.602,   2.746,   2.947,   3.286,   4.073,

  0.258,   0.535,   0.865,   1.071,   1.337,   1.746,   2.120,   2.235, /*v=16*/
  2.382,   2.583,   2.724,   2.921,   3.252,   4.015,

  0.257,   0.534,   0.863,   1.069,   1.333,   1.740,   2.110,   2.224, /*v=17*/
  2.368,   2.567,   2.706,   2.898,   3.222,   3.965,

  0.257,   0.534,   0.862,   1.067,   1.330,   1.734,   2.101,   2.214, /*v=18*/
  2.356,   2.552,   2.689,   2.878,   3.197,   3.922,

  0.257,   0.533,   0.861,   1.066,   1.328,   1.729,   2.093,   2.205, /*v=19*/
  2.346,   2.539,   2.674,   2.861,   3.174,   3.883,

  0.257,   0.533,   0.860,   1.064,   1.325,   1.725,   2.086,   2.197, /*v=20*/
  2.336,   2.528,   2.661,   2.845,   3.153,   3.849,

  0.257,   0.532,   0.859,   1.063,   1.323,   1.721,   2.080,   2.189, /*v=21*/
  2.328,   2.518,   2.649,   2.831,   3.135,   3.819,

  0.256,   0.532,   0.858,   1.061,   1.321,   1.717,   2.074,   2.183, /*v=22*/
  2.320,   2.508,   2.639,   2.819,   3.119,   3.792,

  0.256,   0.532,   0.858,   1.060,   1.319,   1.714,   2.069,   2.177, /*v=23*/
  2.313,   2.500,   2.629,   2.807,   3.104,   3.768,

  0.256,   0.531,   0.857,   1.059,   1.318,   1.711,   2.064,   2.172, /*v=24*/
  2.307,   2.492,   2.620,   2.797,   3.091,   3.745,

  0.256,   0.531,   0.856,   1.058,   1.316,   1.708,   2.060,   2.167, /*v=25*/
  2.301,   2.485,   2.612,   2.787,   3.078,   3.725,

  0.256,   0.531,   0.856,   1.058,   1.315,   1.706,   2.056,   2.162, /*v=26*/
  2.296,   2.479,   2.605,   2.779,   3.067,   3.707,

  0.256,   0.531,   0.855,   1.057,   1.314,   1.703,   2.052,   2.158, /*v=27*/
  2.291,   2.473,   2.598,   2.771,   3.057,   3.690,

  0.256,   0.530,   0.855,   1.056,   1.313,   1.701,   2.048,   2.154, /*v=28*/
  2.286,   2.467,   2.592,   2.763,   3.047,   3.674,

  0.256,   0.530,   0.854,   1.055,   1.311,   1.699,   2.045,   2.150, /*v=29*/
  2.282,   2.462,   2.586,   2.756,   3.038,   3.659,

  0.256,   0.530,   0.854,   1.055,   1.310,   1.697,   2.042,   2.147, /*v=30*/
  2.278,   2.457,   2.581,   2.750,   3.030,   3.646,

  0.255,   0.529,   0.851,   1.050,   1.303,   1.684,   2.021,   2.125, /*v=40*/
  2.250,   2.423,   2.542,   2.704,   2.971,   3.551,

  0.254,   0.527,   0.848,   1.045,   1.296,   1.671,   2.000,   2.099, /*v=60*/
  2.223,   2.390,   2.504,   2.660,   2.915,   3.460,

  0.254,   0.526,   0.845,   1.041,   1.289,   1.658,   1.980,   2.076,/*v=120*/
  2.196,   2.358,   2.468,   2.617,   2.860,   3.373,

  0.253,   0.524,   0.842,   1.036,   1.282,   1.645,   1.960,   2.054,/*v=inf*/
  2.170,   2.326,   2.432,   2.576,   2.807,   3.291};

int meanDiffT_DOF[MeanDiff_Cols]=
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     19,20,21,22,23,24,25,26,27,28,29,30,40,60,120,
     1024};

/*----------------------------------------------------------------------------
 * meanDiff: Finds the appropriate confidence interval.
 * Negative means A1 is better, Positive means A2 is better.
 *
 * -1, +1: Undecided; sign based on least minimum.
 * +/-.40, .30, .20, .15, .10, .05, .025, .02, .015, .01, .0075, .0025, .0005
 * All mean we have that level of confidence (so far) that A1 or A2 has the
 * lesser mean.
 *----------------------------------------------------------------------------
 * We assume variances are NOT equal, but we demand the counts ARE equal.
 * We compute a confidence interval for (meanA1 - meanA2), which requires
 * T=t_{\alpha/2} and a standard error, SE=sqrt((s_1^2+s_2^2)/count); where
 * s_1^2 and s_2^2 are the variances of each population.
 *
 * The point estimate is PE=(meanA1-meanA2), and the confidence interval
 * is PE +/- T*SE. To compute DOF on T, we must truncate the value
 * DOF = ((s_1^2+s_2^2)^2) / ((count-1)*(s_1^4 + s_2^4))  to the lowest
 * available DOF in the table. Since we have v_1 and v_2,
 * DOF = ((v_1+v_2)^2) / ((count-1.)*(v_1^2 + v_2^2))
 *----------------------------------------------------------------------------
 */
double meanDiff(double *A1, double *A2, int count,
                double *stddev1, double *stddev2) /* standard deviations */
{
  int i;
  double x, PE, SE, DOF, meanA1=0.0, meanA2=0.0;
  double varA1=0.0, varA2=0.0, lowerB, upperB;
  double* T;

  if (count < 1) return(1.);                          /* dumb call. */
  for (i=0; i<count; i++)                             /* Compute means. */
  {
    meanA1 += A1[i];                                  /* Add to mean. */
    meanA2 += A2[i];                                  /* .. */
  }

  meanA1 /= (double) count;                            /* Divide by count. */
  meanA2 /= (double) count;                            /* .. */

  if (count < 2)                                      /* If no DOF, */
  {
    if (meanA1 < meanA2) return(-1.);                 /* Go by only value. */
    return( 1.);                                      /* .. */
  }

  /* Must compute variances. */
  for (i=0; i<count; i++)
  {
    x = (A1[i]-meanA1);                               /* Find difference. */
    varA1 += x*x;                                     /* Add to variance. */
    x = (A2[i]-meanA2);                               /* Find difference. */
    varA2 += x*x;                                     /* Add to variance. */
  }

  varA1 /= (double) (count-1);                        /* Compute variance. */
  varA2 /= (double) (count-1);                        /* .. */
  *stddev1 = sqrt(varA1);                             /* standard deviation */
  *stddev2 = sqrt(varA2);                             /* standard deviation */

  SE = sqrt((varA1+varA2)/count);                 /* Compute standard error. */
  PE = meanA1 - meanA2;                           /* Point estimate of diff. */

  /*--------------------------------------------------
   * Compute degrees of freedom.
   *-------------------------------------------------- */
  x = (varA1+varA2);                                  /* Add variances. */
  DOF = x*x*(count-1.);                               /* Set the numerator. */
  x = (varA1*varA1 + varA2*varA2);                    /* Compute denominator. */
  DOF /= x;                                           /* Find the DOF. */
  for (i=1; i<MeanDiff_Cols; i++)                     /* Search widths. */
    if (DOF < meanDiffT_DOF[i]) break;           /* Exit when greater found. */

  /*--------------------------------------------------
   * Recall that the first column is the actual conf,
   * so we want 'i' to be one-relative, not zero rel.
   *-------------------------------------------------- */
  T = meanDiffT_Table + i*MeanDiff_Rows;              /* Point at our column. */

  /*--------------------------------------------------
   * Now, in reverse, we form the confidence interval
   * at every level of certainty; and when we find
   * one that does NOT include zero, we know there is
   * a difference in the means. If we don't find one,
   * we can't tell.
   */
  /*-------------------------------------------------- */

  for (i=(MeanDiff_Rows-1); i>0; i--)                 /* Scan backwards. */
  {
    lowerB = PE - T[i]*SE;                  /* Find lower bound of conf int. */
    upperB = PE + T[i]*SE;                  /* Find upper bound of conf int. */
    if (lowerB > 0.0 || upperB < 0.0) break;     /* If zero impossible, exit. */
  }

  if (meanA1 < meanA2)                    /* If we are returning a negative, */
  {
    if (i == 0) return(-1.0);             /* Could not decide. */
    x = meanDiffT_Table[i];               /* Get the confidence interval. */
    x = 1.-2.*x;                          /* Find our surety level. */
    return(0.-x);                         /* Exit with result. */
  }                                       /* End if returning a negative. */

  /* We are returning a positive. */
  if (i==0) return(1.0);                  /* Could not decide.*/
  x = meanDiffT_Table[i];                 /* Get the confidence interval. */
  x = 1.-2.*x;                            /* Find our surety level. */
  return(x);                              /* Exit with A2 the greater. */

} /* END *** meanDiff *** */

FILE *OpenInputFile(char *file)
{
   FILE *fp;
   if (!strcmp(file, "stdin"))
      fp = stdin;
   else if (!strcmp(file, "stderr"))
      fp = stderr;
   else
   {
      fp = fopen(file, "r");
      if (!fp)
      {
         fprintf(stderr, "Unable to open input file %s!\n", file);
         exit(-1);
      }
   }
   return(fp);
}

FILE *OpenOutputFile(char *file)
{
   FILE *fp;
   if (!strcmp(file, "stdout"))
      fp = stdout;
   else if (!strcmp(file, "stderr"))
      fp = stderr;
   else
   {
      fp = fopen(file, "w");
      if (!fp)
      {
         fprintf(stderr, "Unable to open output file %s!\n", file);
         exit(-1);
      }
   }
   return(fp);
}

void PrintUsage(char *name, char *arg, int i)
{
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' ON %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s <flags> ; flags include:\n", name);
   fprintf(stderr,
           "   -# : # of reps of each timing consecutive output\n");
   fprintf(stderr,
           "   -i <file> : input file for processed results [stdin]\n");
   fprintf(stderr,
           "   -o <file> : output file for processed results [stdout]\n");
   fprintf(stderr, "   -c \"comment line\"\n");
   fprintf(stderr,"   -s <sepstr> : cols in output seperated by this string\n");
   fprintf(stderr,"   -H # <h1> ... <h#> : grab columns with header strings\n");
   fprintf(stderr,
           "   -C # <c1> ... <c#> : grab provided columns numbers (1-based)\n");
   fprintf(stderr,
           "   -S <start> <stride> : grab strided rows starting at start\n");
   exit (i ? i : -1);
}

#define DWAVG 0  /* do what */
#define DWMIN 1
#define DWMAX 2

int GetFlags  /* returns: DOWHAT (AVG,MIN,MAX) */
(
   int nargs, char **args,
   FILE **fpin,         /* input file (stdin) */
   FILE **fpout,        /* output file (stdout) */
   char **cmnt,         /* header comment ("No header specified") */
   int *nrep,           /* # of repititions (1) */
   int *start,          /* outline line to start at (1) */
   int *stride,         /* stride between lines (1) */
   int *ncols,          /* # of cols to make into vectors */
   int **cols,          /* NULL: specced by name, else list of col #s to use */
   char ***colnams      /* NULL: specced by #, else list of headers to get */
)
{
   int i, dow = DWAVG;
   *fpin = stdin;
   *fpout = stdout;
   *cmnt = "No  header specified";
   *nrep = *start = *stride = 1;
   *ncols = 0;
   *cols = NULL;
   *colnams = NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "no '-' preceding flag!", i);
      switch(args[i][1])
      {
         case 'c':
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -c ", i-1);
            *cmnt = args[i];
            break;
         case 's':
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -s ", i-1);
            SEPSTR = args[i];
            break;
         case '#':
            if (args[i][2] == '>')
               dow = DWMAX;
            else if (args[i][2] == '<')
               dow = DWMIN;
            else if (args[i][2] == '+')
               dow = DWAVG;
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -# ", i-1);
            *nrep = atoi(args[i]);
            break;
         case 'o':
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -# ", i-1);
            *fpout = OpenOutputFile(args[i]);
            break;
         case 'S':
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -S", i-1);
            *start = atoi(args[i]);
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -S", i-1);
            *stride = atoi(args[i]);
            break;
         case 'i':
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -i ", i-1);
            *fpin = OpenInputFile(args[i]);
            break;
         case 'C':  /* col numbers: -C # <c1> ... <c#> */
         {
            int j, n, *cn;
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -C ", i-1);
            n = atoi(args[i]);
            cn = malloc(sizeof(int)*n);
            assert(cn);
            for (j=0; j < n; j++)
            {
               if (++i >= nargs)
                  PrintUsage(args[0], "out of flags in -C ", i-1);
               cn[j] = atoi(args[i]) - 1;
            }
            if (*colnams)
                free(*colnams);  /* they are exclusive */
            *colnams = NULL;
            *ncols = n;
            *cols = cn;
         }
            break;
         case 'H':  /* column header strings: -H # <h1> ... <h#> */
         {
            int j, n;
            char **cn;

            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -H ", i-1);
            n = atoi(args[i]);
            cn = malloc(sizeof(char*)*n);
            assert(cn);
            for (j=0; j < n; j++)
            {
               if (++i >= nargs)
                  PrintUsage(args[0], "out of flags in -H ", i-1);
               cn[j] = args[i];
            }
            if (*cols)
               free(*cols);
            *cols = NULL;   /* exlusive */
            *ncols = n;
            *colnams = cn;
         }
            break;
         default :
            PrintUsage(args[0], args[i], i);
      }                 /* end of case over flags */
   }                    /* end of for over flags */
   return(dow);
}
double *GetMaxMinAvgSample(int n, int nrep, double *mf0)
/*
 * Takes an n*nrep length nrep-major array, and returns an 3*n-length array
 * containing the minimum in the first $n$ elts, the avg in the next $n$ elts,
 * and the max in the last $n$ elts.
 */
{
   double *mfS, *mfL, *mfA, max, min, avg;
   const double dreps = 1.0 / (double)nrep;
   int i, k;

   mfS = malloc(3*n*sizeof(double));
   assert(mfS);
   mfA = mfS + n;
   mfL = mfA + n;

   for (i=0; i < n; i++)
   {
      min = avg = max = *mf0++;
      for (k=1; k < nrep; k++, mf0++)
      {
         if (max < *mf0)
            max = *mf0;
         if (min > *mf0)
            min = *mf0;
         avg += *mf0;
      }
      mfL[i] = max;
      mfS[i] = min;
      mfA[i] = avg * dreps;
   }
   return(mfS);
}

int WordsAreDiff(WORDS *wpA, WORDS *wpB)
/*
 * RETURNS word with difference (starting from 1), 0 if same
 */
{
   int i;
   for(i=1; wpA && wpB; i++, wpA = wpA->next, wpB = wpB->next)
      if (strcmp(wpA->word, wpB->word))
         return(i);
   return(wpA == wpB ? 0 : i);  /* not same if one longer */
}

WORDQ *ConsolidateDupWordQ(WORDQ **wqA0, WORDQ **wqB0)
/*
 * Finds duplicate WORDQs in A & B, and removes them from their data
 * structures and puts them on a new consolicated WORDQ, which is returned
 */
{
   WORDQ *wqA, *wqB;            /* gen ptr for A & B Qs */
   WORDQ *wqAp, *wqBp;          /* prev ptrs for A & B Qs */
   WORDQ *wqD=NULL, *wqd;      /* wqD (wqd) is base (gen) ptr for dup Q */

/*
 * Remove all common Qs at front; this requires changing wqA/B
 */
   wqA = *wqA0;
   wqB = *wqB0;
   while (wqA && wqB && !WordsAreDiff(wqA->wlist, wqB->wlist))
   {
      if (!wqD)
         wqD = wqA;
      else
         wqd->next = wqA;
      wqd = wqA;
      wqA = wqA->next;
      wqd->next = NULL;
      wqB = KillWordQ(wqB);
   }
   *wqA0 = wqA;
   *wqB0 = wqB;
   if (!wqA || !wqB)
      return(wqD);
/*
 * Search remainder of Q for dupls; this does not affect 1st ptrs of input Qs
 */
   wqAp = wqA;
   wqA = wqA->next;
   wqBp = wqB;
   wqB = wqB->next;
   while (wqA && wqB)
   {
      if (!WordsAreDiff(wqA->wlist, wqB->wlist))
      {
         if (!wqD)
            wqD = wqA;
         else
            wqd->next = wqA;
         wqd = wqA;
         wqAp->next = wqA = wqA->next;
         wqd->next = NULL;
         wqBp->next = wqB = KillWordQ(wqB);
      }
      else
      {
         wqAp = wqA;
         wqA = wqA->next;
         wqBp = wqB;
         wqB = wqB->next;
      }
   }
   return(wqD);
}

WORDQ *KillWordQByName(WORDQ *wq, char *name)
{
   WORDQ *wqb = wq, *wqP=NULL;
   while (wq)
   {
      if (strcmp(wq->wlist->word, name))
      {
         wqP = wq;
         wq = wq->next;
      }
      else /* names match, so kill this WORDQ */
      {
         if (wqP)
         {
            wq = wqP->next = KillWordQ(wq);
            wqP = wq;
         }
         else
            wq = wqb = KillWordQ(wq);
      }
   }
   return(wqb);
}

void PrintWordQsWithSep(FILE *fpout, WORDQ *wqb, char *suff, char *sep)
{
   while(wqb)
   {
      PrintWordsWithSep(fpout, wqb->wlist, suff, sep);
      wqb = wqb->next;
   }
}

WORDQ *GetMassagedQ(FILE *fpin, int start, int stride, int verb)
/*
 * Reads in file, then massages it into standard format:
 * 1. Labels/names (row 1) are upcased
 * 2. Lines w/o the proper number of cols are deleted
 * 3. Stride and start are taken care of
 * RETURNS: massaged queue.
 */
{
   int i;
   WORDQ *wqb, *wq;    /* base and general ptr, resp. */

   wqb = ReadTimerFile(fpin);
   if (!wqb)
      return(NULL);
   UpCaseWords(wqb->wlist);        /* upcase labels for ease of strcmp */
   wqb = KillBadLines(wqb, verb);
   if (start == 1 && stride == 1)
      return(wqb);
   for (i=0,wq=wqb; i < start; i++, wq = wq->next);
   wq = GetStridedQ(wq, stride);
   wq = GetWordQ(DupWordList(wqb->wlist), wq);  /* put headers back on */
   KillAllWordQ(wqb);
   return(wq);
}

int WordNumberInList
(
   WORDS *wl,           /* list of words to search through */
   char *mtch           /* word to search for in wl */
)
{
   int i;

   for (i=0; wl; i++, wl = wl->next)
   {
      if (!strcmp(mtch, wl->word))
         return(i);
   }
   return(-1);
}

int *ColNames2ColNums   /* returns N-len array of column #s to fetch */
(
   WORDS *headb,        /* headers of columns */
   int N,               /* # of names */
   char **names
)
{
   int i, *nums, uplen=0;
   char *upnam=NULL;

   nums = malloc(sizeof(int)*N);
   assert(nums);
   for (i=0; i < N; i++)
   {
      int k, len;

      len = strlen(names[i]) + 1;
      if (len > uplen)
      {
         if (upnam)
            free(upnam);
         upnam = malloc(sizeof(char)*len);
      }
      for (k=0; upnam[k] = toupper(names[i][k]); k++);
      nums[i] = WordNumberInList(headb, upnam);
      if (nums[i] < 0)
      {
         fprintf(stderr, "NO SUCH COLUMN HEADER = '%s'\n", names[i]);
         free(nums);
         exit(-1);
      }
   }
   if (upnam)
   free(upnam);
   return(nums);
}

WORDQ *GetColumnsByNums
(
   WORDQ *rowq,    /* original row-wise queue from ATLAS output file */
   int N,          /* # of columns to grab */
   int *cnums      /* N-len array of column indexes to grab */
)
{
   int j;
   WORDQ *wqb, *wqp;

   if (N < 1)
      return(NULL);
   wqb = wqp = GetWordQ(GetColByNum(rowq, cnums[0]), NULL);
   for (j=1; j < N; j++)
   {
      wqp->next = GetWordQ(GetColByNum(rowq, cnums[j]), NULL);
      wqp = wqp->next;
   }
   return(wqb);
}

int CountWords(WORDS *words)
{
   int i;
   for (i=0; words; words = words->next, i++);
   return(i);
}

char GetStringType  /* returns c for char, s for string */
(
   WORDS *wp
)
{
   for (; wp; wp = wp->next)
   {
      if (wp->word[1] != '\0')
         return('s');
   }
   return('c');
}

void PrintOutVectors
(
   int N,               /* # of vectors */
   int nrep,            /* # of repititions within the vectors */
   WORDQ *wqb,          /* N-len queue of vectors stored as WORDS */
   FILE *fpout,         /* stream to print to */
   char *cmnt           /* comment to stick as 1st line of file */
)
/*
 * Prints all vectors to fpout.  The output looks like:
 * #<cmnt>
 * <nvec>      # number of vectors
 * <vec1>
 * <vec2>
 *
 * Each vector consists of:
 * <name>
 * <len> <nrep> [i,d,s,c]    # char gives data type (int, double, string, char)
 * <elt1>
 * <elt2>
 */
{
   int j;
   char pre;
   WORDQ *wq=wqb;

   fprintf(fpout, "#%s\n", cmnt);
   fprintf(fpout, "%d\n", N);
   for (j=0; j < N; j++)
   {
      WORDS *wp = wq->wlist;
      if (WordIsFloat(wp->next->word))
         pre = 'd';
      else if (WordIsInt(wp->next->word))
         pre = 'i';
      else
         pre = GetStringType(wp->next);
      fprintf(fpout, "%s\n%d %d %c\n", wp->word, CountWords(wp->next),
              nrep, pre);
      for (wp=wp->next; wp; wp = wp->next)
          fprintf(fpout, "%s\n", wp->word);
      wq = wq->next;
   }
}

int main (int nargs, char **args)
{
   int dw, nrep, start, stride, ncols;
   int *cols;
   FILE *fpin, *fpout;
   char *cmnt, **colnams;
   WORDQ *wqb, *wq;

   dw = GetFlags(nargs, args, &fpin, &fpout, &cmnt, &nrep, &start, &stride,
                 &ncols, &cols, &colnams);
   wqb = GetMassagedQ(fpin, start, stride, 0);
   if (colnams)
   {
      assert(!cols);
      cols = ColNames2ColNums(wqb->wlist, ncols, colnams);
      free(colnams);
      colnams = NULL;
   }
   assert(cols);
   wq = GetColumnsByNums(wqb, ncols, cols);
   KillAllWordQ(wqb);             /* into only selected column-major queues */
   PrintOutVectors(ncols, nrep, wq, fpout, cmnt);
   KillAllWordQ(wq);
   if (fpin != stdin)
      fclose(fpin);
   free(cols);
   return(0);
}
