/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_bitvec.h"
#define ATL_GLOBELT 1
#include "atlas_amm.h"
#undef ATL_GLOBELT
/*
 * Returns a bitvec where set bits are errors above tolerance tol.
 * Entries are ordered row-major for ease of printing
 * Combine with ATL_print2dBV for pictoral error report.
 * flg is bitvec:
 *  0 : if set, a is matrix A, else B or C
 *  1 : if set, a is matrix B, else A or C
 *  2 : if set, A is Upper-, else Lower-symmetric/hermitian
 *  3 : if set, A Hermitian, else symmetric
 */
void *Mjoin(PATL,ipsycmpBV)(int verb, double tol, ipinfo_t *ip, int flg,
                            const TYPE *a, ATL_iptr_t N, const TYPE *A)
{
   ATL_BV_t *bv=NULL;

   if (tol < 0.0)
      tol = -tol;
   if (flg & 1) /* a points to gemm's A matrix in ip */
   {
      ATL_iptr_t i, k;  /* block ctr along K dim */
      const ATL_iptr_t lda2 = ip->lda SHIFT, kb=ip->kb, nfmblks=ip->nfmblks;
      const ATL_iptr_t nkblks = ip->nfkblks, nmblks=nfmblks+ip->npmblks;
      for (k=0; k < N; k++)
      {
         for (i=0; i < N; i++)
         {
         #ifdef TCPLX
            const TYPE *aa;
            double diff, idiff, rc, ic;
            TYPE da[2];

            IdxAwElt_ip(ip, a, i, k, da);
            if (flg&4) /* Upper matrix */
            {
               if (i <= k) /* In upper portion */
               {
                  rc = A[i+i+k*lda2];
                  ic = A[1+i+i+k*lda2];
               }
               else
               {
                  rc = A[k+k+i*lda2];
                  if (flg&8) /* Hermitian, not symmetric! */
                     ic = -A[1+k+k+i*lda2];
                  else
                     ic = A[1+k+k+i*lda2];
               }
            }
            else /* Lower matrix */
            {
               if (i >= k) /* In Lower portion */
               {
                  rc = A[i+i+k*lda2];
                  ic = A[1+i+i+k*lda2];
               }
               else
               {
                  rc = A[k+k+i*lda2];
                  if (flg&8) /* Hermitian, not symmetric! */
                     ic = -A[1+k+k+i*lda2];
                  else
                     ic = A[1+k+k+i*lda2];
               }
            }
            if (i == k && (flg&8)) /* hermitian matrix has diag imag zero */
               ic = ATL_rzero;
            diff = *da - rc;
            idiff = da[1] - ic;
            if (diff < 0.0)
               diff = -diff;
            if (idiff < 0.0)
               idiff = 0.0;
            if (diff > tol || idiff > tol)
            {
               if (!bv)
                  bv = ATL_NewBV(N*N);
               if (verb > 1)
                  printf("A(%d,%d)=[%e,%e];  expected=[%e,%e]\n", i, k,
                         *da, da[1], rc, ic);
               ATL_SetBitBV(bv, i*N+k);
            }
         #else
            double tst, diff;
            ATL_iptr_t idx, lda=ip->lda;
            if (flg&4) /* Upper */
               idx = (i <= k) ? i+k*lda : k+i*lda;
            else
               idx = (i >= k) ? i+k*lda : k+i*lda;
            tst = diff = IdxAwElt_ip(ip, a, i, k);
            diff -= A[idx];
            if (diff < 0.0)
               diff = -diff;

            if (diff > tol)
            {
               if (!bv)
                  bv = ATL_NewBV(N*N);
               if (verb > 1)
                  printf("A(%lu,%lu)=%e;  expected=%e, diff=%e\n",
                         (unsigned long)i, (unsigned long)k, tst, A[idx], diff);
               ATL_SetBitBV(bv, i*N+k);
            }
         #endif
         }
      }
/*
 *    Check zero padding around whole matrix
 */
      i = ip->nmuF*ip->mu - ip->mF;  /* num of extra rows that must be zero */
      if (i) /* we've got i zeros to padd to mul of mu */
      {
         ATL_iptr_t KK = ip->nfkblks*ip->kb + ip->KB0;
         const ATL_iptr_t MM = ip->nfmblks*ip->mb + ip->npmblks*ip->pmb;
         const ATL_iptr_t M = MM - i;
         for (k=0; k < KK; k++)
         {
            for (i=M; i < MM; i++)
            {
            #ifdef TCPLX
               TYPE da[2];
               IdxAwElt_ip(ip, a, i, k, da);
               if (*da != ATL_rzero || da[1] != ATL_rzero)
               {
                  fprintf(stderr, "MU PAD VIOLATION at (%lu,%lu)=(%g,%g), "
                          "not 0!\n", (unsigned long)i, (unsigned long)k,
                          *da, da[1]);
               }
            #else
               TYPE tst;
               tst = IdxAwElt_ip(ip, a, i, k);
               if (tst != ATL_rzero)
               {
                  fprintf(stderr, "MU PAD VIOLATION at (%lu,%lu)=%g, not 0!\n",
                          (unsigned long)i, (unsigned long)k, tst);
               }
            #endif
            }
         }
      }
   }
   else if (flg & 2) /* a points to gemm's B matrix in ip */
   {
      ATL_iptr_t k, j;  /* block ctr along K dim */
      const ATL_iptr_t ldb2 = ip->ldb SHIFT, kb=ip->kb, nfmblks=ip->nfmblks;
      const ATL_iptr_t nkblks = ip->nfkblks, nmblks=nfmblks+ip->npmblks;
      for (j=0; j < N; j++)
      {
         for (k=0; k < N; k++)
         {
         #ifdef TCPLX
            const TYPE *aa;
            double diff, idiff, rc, ic;
            TYPE da[2];

            IdxBwElt_ip(ip, a, k, j, da);
            if (flg&4) /* Upper matrix */
            {
               if (k <= j) /* In upper portion */
               {
                  const ATL_iptr_t IA=k+k+j*ldb2;
                  rc = A[IA];
                  if (j == k && (flg&8))  /* herm diag imag is zero */
                     ic = ATL_rzero;
                  else
                     ic = A[IA+1];
               }
               else
               {
                  const ATL_iptr_t IA=j+j+k*ldb2;
                  rc = A[IA];
                  if (flg&8) /* Hermitian, not symmetric! */
                     ic = -A[IA+1];
                  else
                     ic = A[IA+1];
               }
            }
            else /* Lower matrix */
            {
               if (k >= j) /* In Lower portion */
               {
                  const ATL_iptr_t IA=k+k+j*ldb2;
                  rc = A[IA];
                  if (j == k && (flg&8))  /* herm diag imag is zero */
                     ic = ATL_rzero;
                  else
                     ic = A[IA+1];
               }
               else
               {
                  const ATL_iptr_t IA=j+j+k*ldb2;
                  rc = A[IA];
                  if (flg&8) /* Hermitian, not symmetric! */
                     ic = -A[IA+1];
                  else
                     ic = A[IA+1];
               }
            }
            if (k == j && (flg&8)) /* hermitian matrix has diag imag zero */
               ic = ATL_rzero;
            diff = *da - rc;
            idiff = da[1] - ic;
            if (diff < 0.0)
               diff = -diff;
            if (idiff < 0.0)
               idiff = 0.0;
            if (diff > tol || idiff > tol)
            {
               if (!bv)
                  bv = ATL_NewBV(N*N);
               if (verb > 1)
                  printf("B(%d,%d)=[%e,%e];  expected=[%e,%e]\n", k, j,
                         *da, da[1], rc, ic);
               ATL_SetBitBV(bv, k*N+j);
            }
         #else
            double tst, diff;
            ATL_iptr_t idx, ldb=ip->ldb;
            if (flg&4) /* Upper */
               idx = (j >= k) ? k+j*ldb : j+k*ldb;
            else
               idx = (j <= k) ? k+j*ldb : j+k*ldb;
            tst = diff = IdxBwElt_ip(ip, a, k, j);
            diff -= A[idx];
            if (diff < 0.0)
               diff = -diff;

            if (diff > tol)
            {
               if (!bv)
                  bv = ATL_NewBV(N*N);
               if (verb > 1)
                  printf("A(%lu,%lu)=%e;  expected=%e, diff=%e\n",
                         (unsigned long)k, (unsigned long)j, tst, A[idx], diff);
               ATL_SetBitBV(bv, k*N+j);
            }
         #endif
         }
      }
/*
 *    Check zero padding around whole matrix
 */
      j = ip->nnuF*ip->nu - ip->nF;  /* num of extra rows that must be zero */
      if (j) /* we've got i zeros to padd to mul of nu */
      {
         ATL_iptr_t KK = ip->nfkblks*ip->kb + ip->KB0;
         const ATL_iptr_t NN = ip->nfnblks*ip->nb + ip->npnblks*ip->pnb;
         const ATL_iptr_t N = NN - j;
         for (j=N; j < NN; j++)
         {
            for (k=0; k < KK; k++)
            {
            #ifdef TCPLX
               TYPE da[2];
               IdxBwElt_ip(ip, a, k, j, da);
               if (*da != ATL_rzero || da[1] != ATL_rzero)
               {
                  fprintf(stderr, "NU PAD VIOLATION at (%lu,%lu)=(%g,%g), "
                          "not 0!\n", (unsigned long)k, (unsigned long)j,
                          *da, da[1]);
               }
            #else
               TYPE tst;
               tst = IdxBwElt_ip(ip, a, k, j);
               if (tst != ATL_rzero)
               {
                  fprintf(stderr, "MU PAD VIOLATION at (%lu,%lu)=%g, not 0!\n",
                          (unsigned long)k, (unsigned long)j, tst);
               }
            #endif
            }
         }
      }
   }
   return(bv);
}
