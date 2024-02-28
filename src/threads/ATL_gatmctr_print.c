#include "atlas_gatmctr.h"  /* see to understand usage & data struct */
void ATL_gatmctr_print(FILE *fpout, void *ac)
{
   long *gp=ac, *tp;
   const long P=gp[ATL_GAC_P], flg=gp[ATL_GAC_FLG];
   long i;

   if (GAC_FLAG_SET(flg, ATL_GAC_NOPRV))
   {
      long *lp = ATL_IncBySafeLS(gp), *lac;
      const long incLow = gp[2], incCnt=gp[3];
      fprintf(fpout,
              "gatmctr-PUB, P=%ld, iLow=%ld, iCnt=%ld, off=%ld, flg=%lx\n",
              P, incLow, incCnt, gp[4], flg);
      lac = ATL_AddBytesPtr(lp, incLow);
      for (i=0; i < P; i++, lac=ATL_AddBytesPtr(lac, incCnt))
         fprintf(fpout, "   %ld: low=%ld, cnt=%ld\n", i, lp[i], *lac);
      fprintf(fpout, "   L1: low=0, cnt=%ld\n", *lac);
   }
   else if (GAC_FLAG_SET(flg, ATL_GAC_NOPUB))
   {
      long *lp = ATL_IncBySafeLS(gp);
      fprintf(fpout, "gatmctr-PRV, P=%ld, off=%ld, flg=%lx\n", P, gp[2], flg);
      for (i=0; i < P; i++, lp=ATL_IncBySafeLS(lp))
         fprintf(fpout, "   %ld: low=%ld, cnt=%ld\n", i, lp[1], *lp);
   }
   else
   {
      const long inc=gp[ATL_GAC_INC];
      fprintf(fpout, "gAtmCtr, P=%ld, inc=%ld, N=%ld, flg=%lx:\n",
              P, inc, gp[ATL_GAC_N], flg);
      tp = ATL_AddBytesPtr(gp, ATL_SAFELS);
      for (i=0; i < P; i++)
      {
         long *prv=tp, *lowPubP=ATL_IncBySafeLS(tp),
                      *cntPubP=ATL_IncBySafeLS(lowPubP);
         fprintf(fpout,
            "   %ld: prv=[%ld, %ld], pub=[%ld, %ld], lastPubCnt=%ld, lck=%p\n",
                 i, *prv, *lowPubP-1, *lowPubP, *cntPubP, prv[1],
                 ATL_IncBySafeLS(cntPubP));
         tp = ATL_AddBytesPtr(tp, inc);
      }
   }
}
