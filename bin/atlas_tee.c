/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */

#include <stdio.h>
#include <stdlib.h>
int main(int nargs, char *args[])
{
   char ln[512];
   FILE *fp;
   int i;

   if (nargs != 2)
   {
      fprintf(stderr, "USAGE : %s <fnam>\n", args[0]);
      exit(-1);
   }
   fp = fopen(args[1], "a");
   if (fp == NULL)
   {
      fprintf(stderr, "%s: unable to open file %s\n", args[0], args[1]);
      exit(-1);
   }
   while(fgets(ln, 512, stdin))
   {
      fprintf(stdout, "%s", ln);
      fprintf(fp, "%s", ln);
   }
   i = ferror(fp);
   fclose(fp);
   return(i);
}
