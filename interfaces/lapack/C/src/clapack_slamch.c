/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2014, 2009 R. Clint Whaley
 */
/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 R. Clint Whaley
 */
#define SREAL
#include "clapack.h"
#include "atlas_slamch.h"

float clapack_slamch(char what)
{
   switch(what)
   {
   case 'r':
   case 'R':
      return(ATL_slaROUND);
   case 's':
   case 'S':
      return(ATL_slaSAFMIN);
   case 'o':
   case 'O':
      return(ATL_slaOVERTHRESH);
   case 'u':
   case 'U':
      return(ATL_slaUNDERTHRESH);
   case 'l':
   case 'L':
      return(ATL_slaMAXEXP);
   case 'm':
   case 'M':
      return(ATL_slaMINEXP);
   case 'n':
   case 'N':
      return(ATL_slaMANTDIG);
   case 'p':
   case 'P':
      return(ATL_slaPRECISION);
   case 'b':
   case 'B':
      return(ATL_slaBASE);
   case 'e':
   case 'E':
      return(ATL_slaEPSILON);
   default:
      return(0.0);
   }
   return(0.0);
}
