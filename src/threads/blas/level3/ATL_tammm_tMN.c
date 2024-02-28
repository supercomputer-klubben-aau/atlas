#define ATL_GLOBIDX  1
#include "atlas_misc.h"
#define ATL_ESTNCTR 1
#include "atlas_level1.h"
#include "atlas_tlvl3.h"
#include Mstr(Mjoin(ATLAS_PRE,ipmen_view.h))

#ifdef TCPLX
static void INLINE Mjoin(PATL,amm_tMN_k)
   (ipinfo_t *ip, size_t k, size_t bet, const TYPE *A, const TYPE *B,
    TYPE *rA, TYPE *iA, TYPE *rB, TYPE *iB, TYPE *rC, TYPE *iC)
{
   const size_t nfkblks = ip->nfkblks;
   ATL_CUINT nmu=ip->nmuF, nnu=ip->nnuF;
   ATL_UINT kb = ip->kb;
   if (k != nfkblks || kb == ip->kb0)  /* normal kb size */
   {
      ammkern_t amm_b1=ip->amm_b1, amm_bn=ip->amm_bn;
      ip->b2blk(kb, ip->nF, ip->ONE, B+ip->incBk*k, ip->ldb, rB, iB);
      ip->a2blk(kb, ip->mF, ip->ONE, A+ip->incAk*k, ip->lda, rA, iA);
      if (bet)
      {
         amm_bn(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm_b1(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
      }
      else
      {
         ammkern_t amm=ip->amm_b0;
         amm(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
      }
      amm_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
      amm_b1(nmu, nnu, kb, iA, rB, iC, iA, iB, rC);
   }
   else /* KB0/kb0 size */
   {
      ammkern_t amm_b1=ip->ammK1_b1, amm_bn=ip->ammK1_bn;
      kb = ip->kb0;
      ip->b2blk(kb, ip->nF, ip->ONE, B+ip->incBk*nfkblks, ip->ldb, rB, iB);
      ip->a2blk(kb, ip->mF, ip->ONE, A+ip->incAk*nfkblks, ip->lda, rA, iA);
      kb = ip->KB0;
      if (bet)
      {
         amm_bn(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm_b1(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
      }
      else
      {
         ammkern_t amm=ip->ammK1_b0;
         amm(nmu, nnu, kb, iA, iB, rC, rA, iB, iC);
         amm(nmu, nnu, kb, rA, iB, iC, rA, rB, rC);
      }
      amm_bn(nmu, nnu, kb, rA, rB, rC, iA, rB, iC);
      amm_b1(nmu, nnu, kb, iA, rB, iC, iA, iB, rC);
   }
}
#else
static void INLINE Mjoin(PATL,amm_tMN_k)
   (ipinfo_t *ip, size_t k, size_t bet, const TYPE *A, const TYPE *B,
    TYPE *pA, TYPE *pB, TYPE *pC)
{
   const size_t nfkblks = ip->nfkblks;
   ATL_UINT kb = ip->kb;
   if (k != nfkblks || kb == ip->kb0)  /* normal kb size */
   {
      ip->b2blk(kb, ip->nF, ATL_rone, B+ip->incBk*k, ip->ldb, pB);
      ip->a2blk(kb, ip->mF, ATL_rone, A+ip->incAk*k, ip->lda, pA);
      if (bet)
         ip->amm_b1(ip->nmuF, ip->nnuF, kb, pA, pB, pC, pA, pB, pC);
      else
         ip->amm_b0(ip->nmuF, ip->nnuF, kb, pA, pB, pC, pA, pB, pC);
   }
   else /* KB0/kb0 size */
   {
      ip->b2blk(ip->kb0, ip->nF, ATL_rone, B+ip->incBk*nfkblks, ip->ldb, pB);
      ip->a2blk(ip->kb0, ip->mF, ATL_rone, A+ip->incAk*nfkblks, ip->lda, pA);
      if (bet)
         ip->ammK1_b1(ip->nmuF, ip->nnuF, ip->KB0, pA, pB, pC, pA, pB, pC);
      else
         ip->ammK1_b0(ip->nmuF, ip->nnuF, ip->KB0, pA, pB, pC, pA, pB, pC);
   }
}
#endif

void Mjoin(PATL,DoWork_amm_tMN)(void *vpp, int rank, int vrank)
{
   ATL_tpool_t *pp=vpp;
   ATL_tamm_tMN_t *pd = pp->PD;
   ipinfo_t *ip = pd->ip;
   const size_t nUblks=pd->nUblks, nblksKp=pd->nblksKp, nKp=pd->nKp;
   const size_t nPblks=nKp*nblksKp;
   size_t namm=0;
   const TYPE *A=pd->A, *B=pd->B;
   #ifdef TCPLX
      TYPE *iA = pd->w + vrank*pd->wsz, *rA, *iB, *rB, *iC, *rC;
   #else
      TYPE *pA = pd->w + vrank*pd->wsz, *pB, *pC;
   #endif

   #ifdef TCPLX
      rA = iA + ip->szA;
      iB = rA + ip->szA;
      iB = ATL_AlignPtr(iB);
      rB = iB + ip->szB;
      iC = rB + ip->szB;
      iC = ATL_AlignPtr(iC);
      rC = iC + ip->szC;
   #else
      pB = pA + ip->szA;
      pB = ATL_AlignPtr(pB);
      pC = pB + ip->szB;
      pC = ATL_AlignPtr(pC);
   #endif
   if (nPblks)
   {
      size_t kctr;
      while ( (kctr = ATL_DecGlobalAtomicCount(pd->PartKctr, vrank)) )
      {
         size_t k = (nKp - kctr)*nblksKp, kmax;
         for (kmax=k+nblksKp; k < kmax; k++)
         {
            #ifdef TCPLX
               Mjoin(PATL,amm_tMN_k)(ip, k, namm, A, B, rA, iA, rB, iB, rC, iC);
            #else
               Mjoin(PATL,amm_tMN_k)(ip, k, namm, A, B, pA, pB, pC);
            #endif
            namm++;
         }
      }
   }
   if (nUblks)
   {
      size_t kctr, kk = nPblks + nUblks;
      while ( (kctr = ATL_DecGlobalAtomicCount(pd->blkUctr, vrank)) )
      {
         size_t k = kk - kctr;
         #ifdef TCPLX
            Mjoin(PATL,amm_tMN_k)(ip, k, namm, A, B, rA, iA, rB, iB, rC, iC);
         #else
            Mjoin(PATL,amm_tMN_k)(ip, k, namm, A, B, pA, pB, pC);
         #endif
         namm++;
      }
   }
   if (namm)
   {
      ATL_mutex_lock(pd->Cmut);
      if (pd->ndone)             /* beta already applied */
      #ifdef TCPLX
         ip->blk2c_b1(ip->mF, ip->nF, ip->alpC, rC, iC, ip->ONE, pd->C,ip->ldc);
      #else
         ip->blk2c_b1(ip->mF, ip->nF, ip->alpC, pC, ATL_rone, pd->C, ip->ldc);
      #endif
      else                        /* I'm first, apply real beta */
      #ifdef TCPLX
         ip->blk2c(ip->mF, ip->nF, ip->alpC, rC, iC, pd->beta, pd->C, ip->ldc);
      #else
         ip->blk2c(ip->mF, ip->nF, ip->alpC, pC, pd->beta, pd->C, ip->ldc);
      #endif
      (pd->ndone)++;
      ATL_mutex_unlock(pd->Cmut);
   }
   #ifdef DEBUG
   printf("%d: namm=%lu, wrk=[%p,%p]\n", vrank, (unsigned long) namm,
          pA, (pA+pd->wsz));
   #endif
}
/*
 * This routine should be used only when M <= LASTMB && N <= LASTNB, and it
 * therefore handles a simple inner-product to update a single block of C.
 * It requires minimal workspace of (szA+szB+szC)*P.
 */
int Mjoin(PATL,tammm_tMN)
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
   size_t ldc
)
{
   void *vp;
   size_t sz;
   ATL_UINT P, p;
   ATL_tamm_tMN_t pd;
   ipinfo_t ip;

   if (M > ATL_VWipgen_LAST_MB || N > ATL_VWipgen_LAST_NB ||
       K < Mmin(ATL_NTHREADS,8)*ATL_VWipgen_LAST_KB)
      return(2);  /* not correct shape! */

   P = Mjoin(PATL,tGetIPInfo_tMN)(&ip, ATL_NTHREADS, TA, TB, M, N, K, alpha,
                                  lda, ldb, beta, ldc);
   #if 0
      pd.nKp = 1;
      pd.nblksKp = ip.nfkblks;
   #elif 0
      pd.nKp = pd.nblksKp = 0;
   #else
      sz = (ip.nfkblks + 1)>>1;
      if (sz < P)
         pd.nKp = pd.nblksKp = 0;
      else if (sz >= (P<<4))
      {
         pd.nKp = 8;
         pd.nblksKp = sz>>3;
      }
      else if (sz >= (P<<3))
      {
         pd.nKp = 4;
         pd.nblksKp = sz>>2;
      }
      else
      {
         pd.nblksKp = 2;
         pd.nKp = sz>>1;
      }
   #endif
   pd.nUblks = (ip.nfkblks+1) - pd.nKp*pd.nblksKp;
   p = pd.nKp + pd.nUblks;
   P = Mmin(P,p);
   #ifdef DEBUG
printf("     P=%d, nPblks=%d, nKp=%d, nblksP=%d, nUblks=%d\n", P,
       (int)(pd.nKp*pd.nblksKp), (int)pd.nKp, (int)pd.nblksKp,(int)(pd.nUblks));
   #endif
   pd.wsz = ATL_MulBySize(ip.szA + ip.szB + ip.szC) + 3*ATL_Cachelen-1;
   pd.wsz = ATL_MulByCachelen(ATL_DivByCachelen(pd.wsz));
   vp = malloc(ATL_Cachelen + P*pd.wsz + ATL_MulBySize((ip.mu+ip.mu)*ip.nu));
   if (!vp)
      return(1);
   pd.wsz = ATL_DivBySize((pd.wsz SHIFT));
   pd.w = ATL_AlignPtr(vp);
   pd.ip = &ip;
   pd.beta = beta;
   pd.A = A;
   pd.B = B;
   pd.C = C;
   pd.Cmut = ATL_mutex_init();
   pd.ndone = 0;
   if (pd.nKp)
      pd.PartKctr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nKp, P), pd.nKp, 0);
   else
      pd.PartKctr = NULL;
   if (pd.nUblks)
      pd.blkUctr = ATL_SetGlobalAtomicCount(ATL_EstNctr(pd.nUblks, P),
                                            pd.nUblks, 0);
   else
      pd.blkUctr = NULL;

   ATL_goParallel(P, Mjoin(PATL,DoWork_amm_tMN), NULL, &pd, NULL);
   if (pd.PartKctr)
      ATL_FreeGlobalAtomicCount(pd.PartKctr);
   if (pd.blkUctr)
      ATL_FreeGlobalAtomicCount(pd.blkUctr);
   ATL_mutex_free(pd.Cmut);
   free(vp);
   return(0);
}
