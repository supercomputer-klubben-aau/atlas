#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_UPR,amm_kern.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_sqsyrk.h))
#include Mstr(Mjoin(ATLAS_PRE,amm_umsyrk.h))
#include Mstr(Mjoin(ATLAS_PRE,ipmen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,syrk_view.h))
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromANg_a1.h))
#ifdef TCPLX
   #include Mstr(COPY/Mjoin(ATLAS_PRE,FromAHg_a1.h))
#endif
#include Mstr(COPY/Mjoin(ATLAS_PRE,FromATg_a1.h))
/* #define DEBUG 1 */
double Mjoin(PATL,ipsyrkInfo)
(
   ipinfo_t *ip,        /* output */
   int flag,            /* 0:HERK */
   enum ATLAS_TRANS TA,
   size_t N,
   size_t K,
   size_t lda,
   size_t ldc,
   const SCALAR alpha,
   const SCALAR beta
)
{
   size_t nfnblks;
   unsigned int npnblks, mb, nb, pnb, pmb, mu, nu;
   unsigned int idB=ATL_VIEW_BEST_IDX, nbB=0;
   int i;
   double timB=N*N*K;
   #ifdef TCPLX
      const int CONJ = flag&1;
      enum ATLAS_TRANS TB;
   #else
      const enum ATLAS_TRANS TB = (TA==AtlasNoTrans) ? AtlasTrans:AtlasNoTrans;
   #endif
   for (i=ATL_VIEW_BEST_IDX; i >= 0; i--)
   {
      double tim;
      float *tpf, *stpf;
      ATL_iptr_t nfdiag, ngblk;
      unsigned int nb, pnb;

      nb = ATL_GetVWipmenNB(i);
      if (i && nb+nb > N)
         continue;
      nfdiag = N / nb;
      pnb = N - nb*nfdiag;
      ngblk = ((nfdiag-1)*nfdiag)>>1;
      tim = nb*nb;
      tim *= K;
      tpf  = (float*)(ATL_VIEW_ipmen+ATL_VWipmen_IDXMUL(i));
      stpf = (float*)(ATL_VIEW_syrk+ATL_VWsyrk_IDXMUL(i));
      tim *= (*tpf * 2.0*ngblk) + (*stpf * nfdiag);
/*
 *    For partial block, guess performance poor (half normal) if pnb <= 3*U,
 *    and good (.85 normal) if bigger than that.
 */
      if (pnb)
      {
         unsigned int imm;
         double tg, ts, nf;
         imm = ATL_GetVWipmenAMMI(i);
         ATL_AMM_GetMNU(imm, mu, nu);
         nf = (pnb*1.0)*pnb * K;
         tg = (pnb >= 3*Mmax(mu,nu)) ? .85 : .5;
         ts = (pnb >= 3*ATL_SYRKK_NU) ? .85 : .5;
         tg *= *tpf * (nf+nf);
         ts *= *stpf * nf;
         tim += tg + ts;
      }
      #ifdef DEBUG
         printf("%u: D=(%u,%u), NB=%u(%u) tim=%e\n", i, N, K, nb, pnb, tim);
      #endif
      if (tim > timB)  /* don't keep looking if we are getting slower */
         break;        /* assume get smaller to reduce SYRK cost */
      {
         idB = i;
         nbB = nb;
         nfnblks = nfdiag;
         timB = tim;
      }
   }
   i = ATL_GetVWipmenAMMI(idB);
   ATL_AMM_GetMNU(i, mu, nu);
   if (nbB >= N)
   {
      npnblks = 1;
      nfnblks = 0;
      pmb = mb = ((N+mu-1)/mu)*mu;
      pnb = nb = ((N+nu-1)/nu)*nu;
   }
   else
   {
      mb = nb = nbB;
      pnb = N - nfnblks*nbB;
      if (pnb)
      {
         npnblks = 1;
         pnb = N - nfnblks*nbB;
         pmb = ((pnb+mu-1)/mu)*mu;
         pnb = ((pnb+nu-1)/nu)*nu;
      }
      else
         npnblks = pnb = pmb = 0;
   }
   #ifdef DEBUG
      printf("F%u: D=(%u,%u), NB=%u(%u)\n", idB, N, K, nbB, pnb);
   #endif
   #ifdef TCPLX
      if (CONJ)
         TB = (TA == AtlasNoTrans) ? AtlasConjTrans : AtlasNoTrans;
      else
         TB = (TA == AtlasNoTrans) ? AtlasTrans : AtlasNoTrans;
   #endif
/*
 * Consider if it would be better to use inner-product version that uses
 * only syrk instead of both syrk+gemm.  Need to include copy cost in time
 * estimate, since syrk+gemm requires up to 3 A/B copies, while syrk 1.
 */
   if (N <= ATL_VWsyrk_MAX_NB)
   {
      double timS;
      unsigned int nbS;
      float *tpf;
      for (i=ATL_VWsyrk_NCASES-1; i; i--)
      {
         nbS = ATL_GetVWipmenNB(i);
         if (nbS <= N)
            break;
      }
      if (nbS != N) /* see if this or bigger NB closer to N */
      {
         unsigned int nbL;
         nbL = ATL_GetVWipmenNB(i+1);
         i = (nbL-nb > nb-nbS) ? i+1 : i;
         nbS = ((N+ATL_SYRKK_NU-1)/ATL_SYRKK_NU)*ATL_SYRKK_NU;
      }
      if (nbS >= N)
         return(-1.0);
      tpf = (float*)(ATL_VIEW_syrk+ATL_VWsyrk_IDXMUL(i));
      timS = (*tpf*nbS)*nbS*K;
      #ifdef DEBUG
         printf("tim[G,S]=%e, %e\n", timB, timS);
      #endif
      if (timS <= timB)  /* if syrk wins w/o including copy cost */
         return(-timS);  /* save expense of computing copy cost */
      else               /* see if syrk wins incl copy A/B estimate */
      {
         const ATL_cparr_t *cpA;
         ATL_cpflt_t *cpT;
         double ne, timCP;

         cpA = (TA == AtlasNoTrans) ? ATL_CpyFromANg_a1 : ATL_CpyFromATg_a1;
         cpT = (ATL_cpflt_t*)(cpA+ATL_VWsyrk_MIN_A2BLK+ATL_VWsyrk_MIN_A2BLK+1);
         timCP = *cpT * nbS*K;
         #ifdef DEBUG
            printf("timS=%e+%e\n", timS, timCP);
         #endif
         timS += timCP;
         nbS = ((N+nu-1)/nu)*nu;
         ne = nbS;
         ne *= K;
         ATL_GetVWipmenA2BLK(i);
         cpT = (ATL_cpflt_t*)(cpA+i+i+1);
         timCP = *cpT * ne * 1.1;
         #ifdef DEBUG
            printf("timB=%e+%e\n", timB, timCP);
         #endif
         timB += timCP;
         if (timS <= timB)
            return(-timS);
      }
   }
   Mjoin(PATL,ipmenInfoPop)(ip, idB, TA, TB, N, N, K, lda, lda, ldc, alpha,
                             beta, nfnblks, npnblks, mb, pmb, nfnblks,
                             npnblks, nb, pnb);
   return(timB);
}
