#define ATL_GLOBIDX 1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_tlvl3.h"
#include "atlas_bitvec.h"
#include "atlas_cbc.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,ipgen_view.h))
#include Mstr(Mjoin(ATLAS_PRE,ipmen_view.h))

/*
 * recurs on any standard GEMM interface routine, passed as a function
 * pointer in amm.  flg is a bitflag, meaning if set:
 * 0    1 : Do divide M
 * 1    2 : Do divide N
 * 2    4 : Do divide K
 * 3    8 : stop dividing M when it is < 4*MAXMB
 * 4   16 : stop dividing N when it is < 4*MAXNB
 * 5   32 : stop dividing K when it is < 3*MAXKB
 */
int Mjoin(PATL,ammm_REC)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   size_t M,
   size_t N,
   size_t K,
   const SCALAR alpha,
   const TYPE *A,
   size_t lda,
   const TYPE *B,
   size_t ldb,
   const SCALAR beta,
   TYPE *C,
   size_t ldc,
   ATL_UINT flg,
   int (*amm)(enum ATLAS_TRANS,enum ATLAS_TRANS, size_t, size_t, size_t,
              const SCALAR, const TYPE*, size_t,  const TYPE*, size_t,
              const SCALAR, TYPE*, size_t)
)
{
   size_t Dk, Dm;
   if (amm(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc))
   {
      if ((flg&8) && M <= (ATL_VWipgen_LAST_MB<<2))
         flg &= ~1;
      if ((flg&16) && N <= (ATL_VWipgen_LAST_NB<<2))
         flg &= ~2;
      if ((flg&32) && K <= 3*ATL_VWipgen_LAST_KB)
         flg &= ~4;
/*
 *    Stopping criteria in case something is horribly wrong
 */
      if (K <= ATL_VWipgen_LAST_KB && M <= ATL_VWipgen_LAST_KB &&
          N <= ATL_VWipgen_LAST_MB || !(flg&7))
      {
         printf("ATLAS warning: out-of-workspace causes serial execution!\n");
         Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
         return(0);
      }
/*
 *    Divide K first: it cuts space from both A & B
 */
      if ((flg | 3) == flg)
         Dk = Mmax(M,N);
      else if ((flg&3) == 0)
         Dk = 0;
      else
         Dk = (flg&1) ? M : N;
      Dm = (flg&2) ? N : 0;
      if (K+K >= Dk && (flg&4))
      {
         const int KL = K>>1, KR = K-KL;
         #ifdef TCPLX
            TYPE ONE[2] = {ATL_rone, ATL_rzero};
         #else
            #define ONE ATL_rone
         #endif
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, M, N, KL, alpha, A, lda,
                                          B, ldb, beta, C, ldc, flg, amm));
         A += ((TA == AtlasNoTrans) ? KL*lda : KL)SHIFT;
         B += ((TB == AtlasNoTrans) ? KL : KL*ldb)SHIFT;
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, KR, alpha, A, lda,
                                     B, ldb, ONE, C, ldc, flg, amm));
         #ifndef TCPLX
            #undef ONE
         #endif
      }
/*
 *    If M largest dim (twice K), cut it instead
 */
      else if (M >= Dm && (flg&1))
      {
         const int ML = M>>1, MR = M-ML;
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, ML, N, K, alpha, A, lda,
                                          B, ldb, beta, C, ldc, flg, amm));
         A += ((TA == AtlasNoTrans) ? ML : ML*lda)SHIFT;
         return(Mjoin(PATL,ammm_REC)(TA, TB, MR, N, K, alpha, A, lda,
                                     B, ldb, beta, C+(ML SHIFT), ldc, flg,amm));
      }
/*
 *    Otherwise, cut N
 */
      else if (flg&2)
      {
         const int NL = N>>1, NR = N-NL;
         ATL_assert(!Mjoin(PATL,ammm_REC)(TA, TB, M, NL, K, alpha, A, lda,
                                          B, ldb, beta, C, ldc, flg, amm));
         B += ((TB == AtlasNoTrans) ? NL*ldb : NL)SHIFT;
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, NR, K, alpha, A, lda, B, ldb,
                                     beta, C+NL*(ldc SHIFT), ldc, flg, amm));
      }
      else  /* can't cut anymore, give up and use serial */
      {
         printf("ATLAS warning: out-of-workspace causes serial execution!\n");
         Mjoin(PATL,ammm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
         return(0);
      }
   }
   return(0);
}

int Mjoin(PATL,tammm)
(
   enum ATLAS_TRANS TA,
   enum ATLAS_TRANS TB,
   ATL_CSZT M,
   ATL_CSZT N,
   ATL_CSZT K,
   const SCALAR alpha,
   const TYPE *A,
   ATL_CSZT lda,
   const TYPE *B,
   ATL_CSZT ldb,
   const SCALAR beta,
   TYPE *C,
   ATL_CSZT ldc
)
{
   size_t nblks;
/*
 * No threading for Level-2 operations, or support for K=1,2
 */
   if (K <= 2 || M < 2 || N < 2)
      return(1);  /* signal failure */
   #if 0  /* need diff views before this makes a lot of sense */
   if (M <= ATL_VWopgen_BEST_MB && N <= ATL_VWopgen_BEST_NB &&
       K <= ATL_VWopgen_BEST_KB)
      return(1);
   #endif
/*
 * _tNK & _tMN require minimal wrkspc, so just return if they fail
 */
   #if 0  /* Need opdNK view for this! */
   if (N <= ATL_rkAMM_LASTNB && K <= ATL_VWopgen_MAX_KB)
      return(Mjoin(PATL,tammm_tNK)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                   beta, C, ldc));
/*
 * Eventually, need to rule out general case isn't better for large M/N,
 * but for now, use this case anytime it would work
 */
   if (N <= ATL_VWipmen_LAST_NB && M <= ATL_VWipmen_LAST_MB &&
       K > (ATL_VWipmen_LAST_KB<<1))
      return(Mjoin(PATL,tammm_tMN)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                   beta, C, ldc));
   #endif
/*
 * Routines that copy one entire matrix may need to cut varying dims in order
 * to fit within workspace requirements.  Use ammm_REC for this.
 */
   if (K <= ATL_VWopgen_MAX_KB)  /* outer product / rank-K with N > LASTNB */
   {
      return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc, 1, Mjoin(PATL,tammm_tK)));
   }

   if (N < (ATL_VWipgen_BEST_NB<<2) && K < (ATL_VWipgen_BEST_KB<<2) &&
       (M>>ATL_NTHRPOW2) > ATL_VWipgen_BEST_MB)
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                     beta, C, ldc, 6, Mjoin(PATL,tammm_sNK)));
   if (M < (ATL_VWipgen_BEST_MB<<2) && K < (ATL_VWipgen_BEST_KB<<2) &&
       (N>>ATL_NTHRPOW2) > ATL_VWipgen_BEST_NB)
         return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                     beta, C, ldc, 5, Mjoin(PATL,tammm_sMK)));
   if (K >= ATL_VWipgen_BEST_KB && M >= 12 && N >= 12)
      return(Mjoin(PATL,ammm_REC)(TA, TB, M, N, K, alpha, A, lda, B, ldb,
                                  beta, C, ldc, 63, Mjoin(PATL,tammm_G)));
   return(1);
}
