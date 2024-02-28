/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#define ATL_WANT_ILCM
#include "atlas_iopt.h"
int ATL_lcm(const int M, const int N)
{
   return(ATL_iLCM(M,N));
}
