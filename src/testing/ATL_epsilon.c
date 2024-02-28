/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1998 Jeff Horner
 * Code contributers : Jeff Horner, R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_tst.h"

TYPE Mjoin(PATL,epsilon)(void)
{
   static TYPE eps;
   const TYPE half=0.5;
   volatile TYPE maxval, f1=0.5;

   do
   {
      eps = f1;
      f1 *= half;
      maxval = 1.0 + f1;
   }
   while (maxval != 1.0);
   return(eps);
}

