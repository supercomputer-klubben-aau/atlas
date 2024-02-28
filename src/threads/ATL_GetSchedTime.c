#ifdef ATL_MUT
   #define ATL_GetSchedTime ATL_GetSchedTime_mut
   #define ATL_TIMEs ATL_AMM_SCHEDmut
   #define NOSCHEDlac 1
   #define NOSCHEDpub 1
   #define NOSCHEDprv 1
   #define NOSCHEDmix 1
#elif defined(ATL_LAC)
   #define ATL_GetSchedTime ATL_GetSchedTime_lac
   #define ATL_TIMEs ATL_AMM_SCHEDlac
   #define NOSCHEDmut 1
   #define NOSCHEDpub 1
   #define NOSCHEDprv 1
   #define NOSCHEDmix 1
#elif defined(ATL_PUB)
   #define ATL_GetSchedTime ATL_GetSchedTime_pub
   #define ATL_TIMEs ATL_AMM_SCHEDpub
   #define NOSCHEDmut 1
   #define NOSCHEDlac 1
   #define NOSCHEDprv 1
   #define NOSCHEDmix 1
#elif defined(ATL_MIX)
   #define ATL_GetSchedTime ATL_GetSchedTime_mix
   #define ATL_TIMEs ATL_AMM_SCHEDmix
   #define NOSCHEDmut 1
   #define NOSCHEDpub 1
   #define NOSCHEDprv 1
   #define NOSCHEDlac 1
#elif defined(ATL_PRV)
   #define ATL_GetSchedTime ATL_GetSchedTime_prv
   #define ATL_TIMEs ATL_AMM_SCHEDprv
   #define NOSCHEDmut 1
   #define NOSCHEDpub 1
   #define NOSCHEDlac 1
   #define NOSCHEDmix 1
#endif
#include "atlas_tsched_time.h"
#ifdef ATL_AMM_NSCHED
float ATL_GetSchedTime(unsigned int idx)
{
   return(ATL_TIMEs[idx]);
}
#endif
