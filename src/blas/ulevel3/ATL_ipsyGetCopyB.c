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
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_a1.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aN.h))
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_aX.h))
#endif
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
cm2am_t Mjoin(PATL,ipsyGetCopyB)
(
   ipinfo_t *ip,
   ATL_CUINT flg, /* bitvec: 0:Upper?; 1:HEMM? */
   cm2am_t *CPT
)
{
   cm2am_t cpN, cpT;
   ATL_UINT icpB;
   const SCALAR alpha = ip->alpB;

   icpB = ATL_GetViewB2BLK(ip->idx);
   icpB += icpB;
   #ifdef TCPLX
      if (alpha[1] == ATL_rzero)
      {
         const TYPE ralp = *alpha;
         if (ralp == ATL_rone)
         {
            cpN = (cm2am_t) ATL_CpyFromANg_a1[icpB];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromAHg_a1[icpB]:ATL_CpyFromATg_a1[icpB]);
         }
         else if (ralp == ATL_rnone)
         {
            cpN = (cm2am_t) ATL_CpyFromANg_aN[icpB];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromAHg_aN[icpB]:ATL_CpyFromATg_aN[icpB]);
         }
         else
         {
            cpN = (cm2am_t) ATL_CpyFromAHg_aX[icpB];
            cpT = (cm2am_t) ((flg&2) ?
                             ATL_CpyFromAHg_aX[icpB]:ATL_CpyFromATg_aX[icpB]);
         }
      }
      else
      {
         cpN = (cm2am_t) ATL_CpyFromANg_aX[icpB];
         cpT = (cm2am_t) ((flg&2) ?
                          ATL_CpyFromAHg_aX[icpB]:ATL_CpyFromATg_aX[icpB]);
      }
   #else
      if (alpha == ATL_rone)
      {
         cpN = (cm2am_t) ATL_CpyFromANg_a1[icpB];
         cpT = (cm2am_t) ATL_CpyFromATg_a1[icpB];
      }
      else if (alpha == ATL_rnone)
      {
         cpN = (cm2am_t) ATL_CpyFromANg_aN[icpB];
         cpT = (cm2am_t) ATL_CpyFromATg_aN[icpB];
      }
      else
      {
         cpN = (cm2am_t) ATL_CpyFromANg_aX[icpB];
         cpT = (cm2am_t) ATL_CpyFromATg_aX[icpB];
      }
   #endif
   *CPT = cpT;
   return(cpN);
}

