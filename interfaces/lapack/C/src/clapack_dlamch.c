/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2014, 2009 R. Clint Whaley
 */
/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2009 R. Clint Whaley
 */
#define DREAL
#include "clapack.h"
#include "atlas_dlamch.h"

double clapack_dlamch(char what)
{
   switch(what)
   {
   case 'r':
   case 'R':
      return(ATL_dlaROUND);
   case 's':
   case 'S':
      return(ATL_dlaSAFMIN);
   case 'o':
   case 'O':
      return(ATL_dlaOVERTHRESH);
   case 'u':
   case 'U':
      return(ATL_dlaUNDERTHRESH);
   case 'l':
   case 'L':
      return(ATL_dlaMAXEXP);
   case 'm':
   case 'M':
      return(ATL_dlaMINEXP);
   case 'n':
   case 'N':
      return(ATL_dlaMANTDIG);
   case 'p':
   case 'P':
      return(ATL_dlaPRECISION);
   case 'b':
   case 'B':
      return(ATL_dlaBASE);
   case 'e':
   case 'E':
      return(ATL_dlaEPSILON);
   default:
      return(0.0);
   }
   return(0.0);
}
