/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2011 R. Clint Whaley
 */
#include <stdio.h>
#include <assert.h>
#include "atlconf_misc.h"
void PrintUsage(char *name, int iarg, char *flag)
{
   fprintf(stderr, "Unknown flag '%s' in position %d!\n", flag, iarg);
   fprintf(stderr, "USAGE: [-l <lvl>] <gcc candidate>\n");
   exit(iarg);
}
int main(int nargs, char **args)
{
   int lvl=0;  /* 0: is gcc, 1: is gcc 4 but not apple gcc 2: gcc 4.x, with x >= 4 */
   int i;
   char *comp=NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] == '-')
      {
         if (args[i][1] == 'l')
	 {
	    if (++i >= nargs)
	       PrintUsage(args[0], i, "out of arguments");
	    lvl = atoi(args[i]);
	 }
	 else
	    PrintUsage(args[0], i, args[i]);
      }
      else
         comp = args[i];
   }
   assert(comp);
   if (!CompIsGcc(comp))
      return(1);
   if (lvl)
   {
      int icmp, major, minor, patch;

      GetGccVers(comp, &icmp, &major, &minor, &patch);
      #if 0
         fprintf(stderr, "comp='%s': cmp=%d, major=%d, minor=%d, patch=%d\n",
	         comp, icmp, major, minor, patch);
      #endif
      if (icmp || major < 4)
         return(2);
      if (lvl > 1)
         if (minor < 4)
	    return(3);
   }
   printf("%s\n", comp);
   return(0);
}
