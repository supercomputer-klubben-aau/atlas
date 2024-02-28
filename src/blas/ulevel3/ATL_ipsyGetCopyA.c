/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_amm.h"
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_aX.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aN.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_aX.h))
#ifdef TCPLX
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_a1.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aN.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromACg_aX.h))
#endif
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
cm2am_t Mjoin(PATL,ipsyGetCopyA)(ipinfo_t *ip, ATL_CUINT flg, cm2am_t *CPT)
{
   cm2am_t cpN, cpT;
   ATL_UINT icpA;
   const SCALAR alpha = ip->alpA;

   icpA = ATL_GetViewA2BLK(ip->idx);
   icpA += icpA;
   #ifdef TCPLX
      if (alpha[1] == ATL_rzero)
      {
         const TYPE ralp = *alpha;
         if (ralp == ATL_rone)
         {
            cpN = (cm2am_t) ATL_CpyFromATg_a1[icpA];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromACg_a1[icpA]:ATL_CpyFromANg_a1[icpA]);
         }
         else if (ralp == ATL_rnone)
         {
            cpN = (cm2am_t) ATL_CpyFromATg_aN[icpA];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromACg_aN[icpA]:ATL_CpyFromANg_aN[icpA]);
         }
         else
         {
            cpN = (cm2am_t) ATL_CpyFromATg_aX[icpA];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromACg_aX[icpA]:ATL_CpyFromANg_aX[icpA]);
         }
      }
      else
      {
         cpN = (cm2am_t) ATL_CpyFromATg_aX[icpA];
         cpT = (cm2am_t) ((flg&2) ?
                          ATL_CpyFromACg_aX[icpA]:ATL_CpyFromANg_aX[icpA]);
      }
   #else /* real */
      if (alpha == ATL_rone)
      {
         cpN = (cm2am_t) ATL_CpyFromATg_a1[icpA];
         cpT = (cm2am_t) ATL_CpyFromANg_a1[icpA];
      }
      else if (alpha == ATL_rnone)
      {
         cpN = (cm2am_t) ATL_CpyFromATg_aN[icpA];
         cpT = (cm2am_t) ATL_CpyFromANg_aN[icpA];
      }
      else
      {
         cpN = (cm2am_t) ATL_CpyFromATg_aX[icpA];
         cpT = (cm2am_t) ATL_CpyFromANg_aX[icpA];
      }
   #endif
   *CPT = cpT;
   return(cpN);
}
