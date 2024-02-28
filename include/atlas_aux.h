/*
 * Header file for ATLAS's auxiliary routines
 */
#ifndef ATLAS_AUX_H
#define ATLAS_AUX_H
#include "atlas_misc.h"

void ATL_xerbla(int p, char *rout, char *form, ...);
int ATL_lcm(const int M, const int N);
double ATL_flushcache(long long size);
double ATL_walltime(void);
double ATL_cputime(void);

/*
 * Auxiliary routines that come in all four types
 */
void ATL_sgeset(ATL_CINT M, ATL_CINT N, const float alpha,
                const float beta, float *A, ATL_CINT lda);
void ATL_strsetL(ATL_CINT M, ATL_CINT N, const float alpha,
                 const float beta, float *A, ATL_CINT lda);
void ATL_strsetU(ATL_CINT M, ATL_CINT N, const float alpha,
                 const float beta, float *A, ATL_CINT lda);
float ATL_sgemaxnrm(ATL_CINT M, ATL_CINT N, float *A, ATL_CINT lda);
void ATL_sgeput1T_bX(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float beta, float *C, ATL_CSZT ldc);
void ATL_samL2skLNB_bX
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float beta, float *C, ATL_iptr_t ldc);
void ATL_sgeput1T_bN(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float beta, float *C, ATL_CSZT ldc);
void ATL_samL2skLNB_bN
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float beta, float *C, ATL_iptr_t ldc);
void ATL_sgeput1T_b1(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float beta, float *C, ATL_CSZT ldc);
void ATL_samL2skLNB_b1
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float beta, float *C, ATL_iptr_t ldc);
void ATL_sgeput1T_b0(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float beta, float *C, ATL_CSZT ldc);
void ATL_samL2skLNB_b0
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float beta, float *C, ATL_iptr_t ldc);
void ATL_sgeadd(const int M, const int N, const float alpha,
                const float *A, const int lda, const float beta,
                float *C, const int ldc);
void ATL_sgemove(const int M, const int N, const float alpha,
                 const float *A, const int lda, float *C, const int ldc);
void ATL_sgemoveT(const int N, const int M, const float alpha,
                  const float *A, const int lda, float *C, const int ldc);
void ATL_ssyreflect(const enum ATLAS_UPLO Uplo, const int N,
                    float *C, const int ldc);
void ATL_sgecopy(const int M, const int N, const float *A, const int lda,
                 float *C, const int ldc);

void ATL_sgescal(const int M, const int N, const float beta,
                 float *C, const int ldc);
void ATL_stradd
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float beta, float *C, ATL_CINT ldc);
void ATL_strscal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const float alpha,
    float *A, const int lda);
void ATL_shescal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const float alpha,
    float *A, const int lda);

void ATL_sgezero(const int M, const int N, float *C, const int ldc);
void ATL_ssyApAt_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float beta, float *C, ATL_CINT ldc);
void ATL_ssyApAt
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float beta, float *C, ATL_CINT ldc);
void ATL_sgeApBt_NB
   (ATL_CINT M, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *B, ATL_CINT ldb, const float beta, float *C, ATL_CINT ldc);
void ATL_sgeswapT(ATL_CINT M, ATL_CINT N, float *A, ATL_CINT lda,
                  float *B, ATL_CINT ldb);
void ATL_ssqtrans(ATL_CINT N, float *C, ATL_CINT ldc);

void ATL_szero(const int N, float *X, const int incX);
void ATL_sset(const int N, const float alpha, float *X, const int incX);
void ATL_sscal(const int N, const float alpha, float *X, const int incX);
void ATL_scopy(const int N, const float *X, const int incX,
               float *Y, const int incY);
void ATL_scpsc(const int N, const float alpha, const float *X,
               const int incX, float *Y, const int incY);
void ATL_saxpy(const int N, const float alpha, const float *X,
               const int incX, float *Y, const int incY);
void ATL_saxpy_x1_y1(const int N, const float alpha, const float *X,
                     const int incX, float *Y, const int incY);
void ATL_saxpby(const int N, const float alpha, const float *X,
                const int incX, const float beta, float *Y, const int incY);

void ATL_sgeadd_a1_b1
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a1_b1
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_a0_b1
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a0_b1
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_aX_b1
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_aX_b1
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_a1_b0
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a1_b0
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_a0_b0
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a0_b0
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_aX_b0
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_aX_b0
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_a1_bX
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a1_bX
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_a0_bX
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_a0_bX
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);
void ATL_sgeadd_aX_bX
   (ATL_CINT  M, ATL_CINT  N, const float alpha, const float *A,
    ATL_CINT  lda, const float beta, float *C, ATL_CINT  ldc);
void ATL_saxpby_aX_bX
   (ATL_CINT  N, const float alpha, const float *X, ATL_CINT  incX,
    const float beta, float *Y, ATL_CINT  incY);

void ATL_sgemove_a1
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);
void ATL_sgemove_a0
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);
void ATL_sgemove_aX
   (ATL_CINT M, ATL_CINT N, const float alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);

void ATL_sgescal_b1
   (ATL_CINT  M, ATL_CINT  N, const float beta, float *C, ATL_CINT  ldc);
void ATL_sgescal_b0
   (ATL_CINT  M, ATL_CINT  N, const float beta, float *C, ATL_CINT  ldc);
void ATL_sgescal_bX
   (ATL_CINT  M, ATL_CINT  N, const float beta, float *C, ATL_CINT  ldc);

void ATL_dgeset(ATL_CINT M, ATL_CINT N, const double alpha,
                const double beta, double *A, ATL_CINT lda);
void ATL_dtrsetL(ATL_CINT M, ATL_CINT N, const double alpha,
                 const double beta, double *A, ATL_CINT lda);
void ATL_dtrsetU(ATL_CINT M, ATL_CINT N, const double alpha,
                 const double beta, double *A, ATL_CINT lda);
double ATL_dgemaxnrm(ATL_CINT M, ATL_CINT N, double *A, ATL_CINT lda);
void ATL_dgeput1T_bX(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double beta, double *C, ATL_CSZT ldc);
void ATL_damL2skLNB_bX
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double beta, double *C, ATL_iptr_t ldc);
void ATL_dgeput1T_bN(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double beta, double *C, ATL_CSZT ldc);
void ATL_damL2skLNB_bN
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double beta, double *C, ATL_iptr_t ldc);
void ATL_dgeput1T_b1(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double beta, double *C, ATL_CSZT ldc);
void ATL_damL2skLNB_b1
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double beta, double *C, ATL_iptr_t ldc);
void ATL_dgeput1T_b0(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double beta, double *C, ATL_CSZT ldc);
void ATL_damL2skLNB_b0
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double beta, double *C, ATL_iptr_t ldc);
void ATL_dgeadd(const int M, const int N, const double alpha,
                const double *A, const int lda, const double beta,
                double *C, const int ldc);
void ATL_dgemove(const int M, const int N, const double alpha,
                 const double *A, const int lda, double *C, const int ldc);
void ATL_dgemoveT(const int N, const int M, const double alpha,
                  const double *A, const int lda, double *C, const int ldc);
void ATL_dsyreflect(const enum ATLAS_UPLO Uplo, const int N,
                    double *C, const int ldc);
void ATL_dgecopy(const int M, const int N, const double *A, const int lda,
                 double *C, const int ldc);

void ATL_dgescal(const int M, const int N, const double beta,
                 double *C, const int ldc);
void ATL_dtradd
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double beta, double *C, ATL_CINT ldc);
void ATL_dtrscal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const double alpha,
    double *A, const int lda);
void ATL_dhescal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const double alpha,
    double *A, const int lda);

void ATL_dgezero(const int M, const int N, double *C, const int ldc);
void ATL_dsyApAt_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double beta, double *C, ATL_CINT ldc);
void ATL_dsyApAt
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double beta, double *C, ATL_CINT ldc);
void ATL_dgeApBt_NB
   (ATL_CINT M, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *B, ATL_CINT ldb, const double beta, double *C, ATL_CINT ldc);
void ATL_dgeswapT(ATL_CINT M, ATL_CINT N, double *A, ATL_CINT lda,
                  double *B, ATL_CINT ldb);
void ATL_dsqtrans(ATL_CINT N, double *C, ATL_CINT ldc);

void ATL_dzero(const int N, double *X, const int incX);
void ATL_dset(const int N, const double alpha, double *X, const int incX);
void ATL_dscal(const int N, const double alpha, double *X, const int incX);
void ATL_dcopy(const int N, const double *X, const int incX,
               double *Y, const int incY);
void ATL_dcpsc(const int N, const double alpha, const double *X,
               const int incX, double *Y, const int incY);
void ATL_daxpy(const int N, const double alpha, const double *X,
               const int incX, double *Y, const int incY);
void ATL_daxpy_x1_y1(const int N, const double alpha, const double *X,
                     const int incX, double *Y, const int incY);
void ATL_daxpby(const int N, const double alpha, const double *X,
                const int incX, const double beta, double *Y, const int incY);

void ATL_dgeadd_a1_b1
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a1_b1
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_a0_b1
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a0_b1
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_aX_b1
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_aX_b1
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_a1_b0
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a1_b0
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_a0_b0
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a0_b0
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_aX_b0
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_aX_b0
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_a1_bX
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a1_bX
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_a0_bX
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_a0_bX
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);
void ATL_dgeadd_aX_bX
   (ATL_CINT  M, ATL_CINT  N, const double alpha, const double *A,
    ATL_CINT  lda, const double beta, double *C, ATL_CINT  ldc);
void ATL_daxpby_aX_bX
   (ATL_CINT  N, const double alpha, const double *X, ATL_CINT  incX,
    const double beta, double *Y, ATL_CINT  incY);

void ATL_dgemove_a1
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);
void ATL_dgemove_a0
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);
void ATL_dgemove_aX
   (ATL_CINT M, ATL_CINT N, const double alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);

void ATL_dgescal_b1
   (ATL_CINT  M, ATL_CINT  N, const double beta, double *C, ATL_CINT  ldc);
void ATL_dgescal_b0
   (ATL_CINT  M, ATL_CINT  N, const double beta, double *C, ATL_CINT  ldc);
void ATL_dgescal_bX
   (ATL_CINT  M, ATL_CINT  N, const double beta, double *C, ATL_CINT  ldc);

void ATL_cgeset(ATL_CINT M, ATL_CINT N, const float *alpha,
                const float *beta, float *A, ATL_CINT lda);
void ATL_ctrsetL(ATL_CINT M, ATL_CINT N, const float *alpha,
                 const float *beta, float *A, ATL_CINT lda);
void ATL_ctrsetU(ATL_CINT M, ATL_CINT N, const float *alpha,
                 const float *beta, float *A, ATL_CINT lda);
float ATL_cgemaxnrm(ATL_CINT M, ATL_CINT N, float *A, ATL_CINT lda);
void ATL_cgeput1T_bX(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_camL2skLNB_bX
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float *beta, float *C, ATL_iptr_t ldc);
void ATL_cgeput1T_bN(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_camL2skLNB_bN
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float *beta, float *C, ATL_iptr_t ldc);
void ATL_cgeput1T_b1(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_camL2skLNB_b1
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float *beta, float *C, ATL_iptr_t ldc);
void ATL_cgeput1T_b0(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_camL2skLNB_b0
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const float *W,  const float *beta, float *C, ATL_iptr_t ldc);
void ATL_cgeadd(const int M, const int N, const float *alpha,
                const float *A, const int lda, const float *beta,
                float *C, const int ldc);
void ATL_cgemove(const int M, const int N, const float *alpha,
                 const float *A, const int lda, float *C, const int ldc);
void ATL_cgemoveT(const int N, const int M, const float *alpha,
                  const float *A, const int lda, float *C, const int ldc);
void ATL_csyreflect(const enum ATLAS_UPLO Uplo, const int N,
                    float *C, const int ldc);
void ATL_cgecopy(const int M, const int N, const float *A, const int lda,
                 float *C, const int ldc);

void ATL_cgescal(const int M, const int N, const float *beta,
                 float *C, const int ldc);
void ATL_ctradd
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *beta, float *C, ATL_CINT ldc);
void ATL_ctrscal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const float *alpha,
    float *A, const int lda);
void ATL_chescal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const float alpha,
    float *A, const int lda);

void ATL_cgezero(const int M, const int N, float *C, const int ldc);
void ATL_csyApAt_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *beta, float *C, ATL_CINT ldc);
void ATL_csyApAt
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *beta, float *C, ATL_CINT ldc);
void ATL_cgeApBt_NB
   (ATL_CINT M, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *B, ATL_CINT ldb, const float *beta, float *C, ATL_CINT ldc);
void ATL_cgeswapT(ATL_CINT M, ATL_CINT N, float *A, ATL_CINT lda,
                  float *B, ATL_CINT ldb);
void ATL_csqtrans(ATL_CINT N, float *C, ATL_CINT ldc);

void ATL_czero(const int N, float *X, const int incX);
void ATL_cset(const int N, const float *alpha, float *X, const int incX);
void ATL_cscal(const int N, const float *alpha, float *X, const int incX);
void ATL_ccopy(const int N, const float *X, const int incX,
               float *Y, const int incY);
void ATL_ccpsc(const int N, const float *alpha, const float *X,
               const int incX, float *Y, const int incY);
void ATL_caxpy(const int N, const float *alpha, const float *X,
               const int incX, float *Y, const int incY);
void ATL_caxpy_x1_y1(const int N, const float *alpha, const float *X,
                     const int incX, float *Y, const int incY);
void ATL_caxpby(const int N, const float *alpha, const float *X,
                const int incX, const float *beta, float *Y, const int incY);

void ATL_cgeadd_a1_b1
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a1_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_a0_b1
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a0_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_aX_b1
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_aX_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_a1_b0
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a1_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_a0_b0
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a0_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_aX_b0
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_aX_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_a1_bX
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a1_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_a0_bX
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_a0_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgeadd_aX_bX
   (ATL_CINT  M, ATL_CINT  N, const float *alpha, const float *A,
    ATL_CINT  lda, const float *beta, float *C, ATL_CINT  ldc);
void ATL_caxpby_aX_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);

void ATL_cgemove_a1
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);
void ATL_cgemove_a0
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);
void ATL_cgemove_aX
   (ATL_CINT M, ATL_CINT N, const float *alpha, const float *A,
    const int lda, float *C, ATL_CINT ldc);

void ATL_cgescal_b1
   (ATL_CINT  M, ATL_CINT  N, const float *beta, float *C, ATL_CINT  ldc);
void ATL_cgescal_b0
   (ATL_CINT  M, ATL_CINT  N, const float *beta, float *C, ATL_CINT  ldc);
void ATL_cgescal_bX
   (ATL_CINT  M, ATL_CINT  N, const float *beta, float *C, ATL_CINT  ldc);

void ATL_zgeset(ATL_CINT M, ATL_CINT N, const double *alpha,
                const double *beta, double *A, ATL_CINT lda);
void ATL_ztrsetL(ATL_CINT M, ATL_CINT N, const double *alpha,
                 const double *beta, double *A, ATL_CINT lda);
void ATL_ztrsetU(ATL_CINT M, ATL_CINT N, const double *alpha,
                 const double *beta, double *A, ATL_CINT lda);
double ATL_zgemaxnrm(ATL_CINT M, ATL_CINT N, double *A, ATL_CINT lda);
void ATL_zgeput1T_bX(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zamL2skLNB_bX
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double *beta, double *C, ATL_iptr_t ldc);
void ATL_zgeput1T_bN(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zamL2skLNB_bN
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double *beta, double *C, ATL_iptr_t ldc);
void ATL_zgeput1T_b1(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zamL2skLNB_b1
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double *beta, double *C, ATL_iptr_t ldc);
void ATL_zgeput1T_b0(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zamL2skLNB_b0
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const double *W,  const double *beta, double *C, ATL_iptr_t ldc);
void ATL_zgeadd(const int M, const int N, const double *alpha,
                const double *A, const int lda, const double *beta,
                double *C, const int ldc);
void ATL_zgemove(const int M, const int N, const double *alpha,
                 const double *A, const int lda, double *C, const int ldc);
void ATL_zgemoveT(const int N, const int M, const double *alpha,
                  const double *A, const int lda, double *C, const int ldc);
void ATL_zsyreflect(const enum ATLAS_UPLO Uplo, const int N,
                    double *C, const int ldc);
void ATL_zgecopy(const int M, const int N, const double *A, const int lda,
                 double *C, const int ldc);

void ATL_zgescal(const int M, const int N, const double *beta,
                 double *C, const int ldc);
void ATL_ztradd
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *beta, double *C, ATL_CINT ldc);
void ATL_ztrscal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const double *alpha,
    double *A, const int lda);
void ATL_zhescal
   (const enum ATLAS_UPLO Uplo, const int M, const int N, const double alpha,
    double *A, const int lda);

void ATL_zgezero(const int M, const int N, double *C, const int ldc);
void ATL_zsyApAt_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *beta, double *C, ATL_CINT ldc);
void ATL_zsyApAt
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *beta, double *C, ATL_CINT ldc);
void ATL_zgeApBt_NB
   (ATL_CINT M, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *B, ATL_CINT ldb, const double *beta, double *C, ATL_CINT ldc);
void ATL_zgeswapT(ATL_CINT M, ATL_CINT N, double *A, ATL_CINT lda,
                  double *B, ATL_CINT ldb);
void ATL_zsqtrans(ATL_CINT N, double *C, ATL_CINT ldc);

void ATL_zzero(const int N, double *X, const int incX);
void ATL_zset(const int N, const double *alpha, double *X, const int incX);
void ATL_zscal(const int N, const double *alpha, double *X, const int incX);
void ATL_zcopy(const int N, const double *X, const int incX,
               double *Y, const int incY);
void ATL_zcpsc(const int N, const double *alpha, const double *X,
               const int incX, double *Y, const int incY);
void ATL_zaxpy(const int N, const double *alpha, const double *X,
               const int incX, double *Y, const int incY);
void ATL_zaxpy_x1_y1(const int N, const double *alpha, const double *X,
                     const int incX, double *Y, const int incY);
void ATL_zaxpby(const int N, const double *alpha, const double *X,
                const int incX, const double *beta, double *Y, const int incY);

void ATL_zgeadd_a1_b1
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a1_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_a0_b1
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a0_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_aX_b1
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_aX_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_a1_b0
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a1_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_a0_b0
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a0_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_aX_b0
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_aX_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_a1_bX
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a1_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_a0_bX
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_a0_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgeadd_aX_bX
   (ATL_CINT  M, ATL_CINT  N, const double *alpha, const double *A,
    ATL_CINT  lda, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zaxpby_aX_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);

void ATL_zgemove_a1
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);
void ATL_zgemove_a0
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);
void ATL_zgemove_aX
   (ATL_CINT M, ATL_CINT N, const double *alpha, const double *A,
    const int lda, double *C, ATL_CINT ldc);

void ATL_zgescal_b1
   (ATL_CINT  M, ATL_CINT  N, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zgescal_b0
   (ATL_CINT  M, ATL_CINT  N, const double *beta, double *C, ATL_CINT  ldc);
void ATL_zgescal_bX
   (ATL_CINT  M, ATL_CINT  N, const double *beta, double *C, ATL_CINT  ldc);


/*
 * Real-only aux routines
 */
void ATL_strcpypad4L(enum ATLAS_DIAG Diag, ATL_CINT N, const float *A,
                     ATL_CINT lda, float *c, ATL_CINT ldc);
void ATL_strcpypad4U(enum ATLAS_DIAG Diag, ATL_CINT N, const float *A,
                     ATL_CINT lda, float *c, ATL_CINT ldc);
void ATL_ssycpyUNB_aX(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ssycpyUNB_an(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ssycpyUNB_a1(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ssycpyLNB_aX(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ssycpyLNB_an(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ssycpyLNB_a1(ATL_CSZT N, const float alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_dtrcpypad4L(enum ATLAS_DIAG Diag, ATL_CINT N, const double *A,
                     ATL_CINT lda, double *c, ATL_CINT ldc);
void ATL_dtrcpypad4U(enum ATLAS_DIAG Diag, ATL_CINT N, const double *A,
                     ATL_CINT lda, double *c, ATL_CINT ldc);
void ATL_dsycpyUNB_aX(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dsycpyUNB_an(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dsycpyUNB_a1(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dsycpyLNB_aX(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dsycpyLNB_an(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dsycpyLNB_a1(ATL_CSZT N, const double alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);

/*
 * Specialized complex auxiliary routines
 */

void ATL_checpyUNB_aX(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyUNB_ar(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyUNB_an(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyUNB_a1(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyLNB_aX(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyLNB_ar(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyLNB_an(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_checpyLNB_a1(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyUNB_aX(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyUNB_ar(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyUNB_an(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyUNB_a1(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyLNB_aX(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyLNB_ar(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyLNB_an(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_csycpyLNB_a1(ATL_CSZT N, const float *alpha, const float *A,
                      ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_cgeput1H_bX(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_cgeput1H_bN(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_cgeput1H_b1(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_cgeput1H_b0(ATL_CSZT M, ATL_CSZT N, const float *A, ATL_CSZT lda,
                     const float *beta, float *C, ATL_CSZT ldc);
void ATL_cswap_cplx2real(ATL_CUINT N, float *C, ATL_CINT incc,
                         float *R, ATL_CINT incR, float *I, ATL_CINT incI);
void ATL_cheApAc_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *beta, float *C, ATL_CINT ldc);
void ATL_cheApAc
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *beta, float *C, ATL_CINT ldc);
void ATL_cgeApBc_NB
   (ATL_CINT M, ATL_CINT N, const float *A, ATL_CINT lda,
    const float *B, ATL_CINT ldb, const float *beta, float *C, ATL_CINT ldc);
void ATL_ccplxdivide(ATL_CINT N, const float *b, float *X, ATL_CINT incX,
                     float *Y, ATL_CINT incY);
void ATL_ccplxinvert
   (const int N, float *X, const int incX, float *Y, const int incY);

void ATL_chereflect(const enum ATLAS_UPLO Uplo, const int N,
                    float *C, const int ldc);
void ATL_cscalConj
   (const int N, const float *alpha, float *X, const int incX);
void ATL_ccopyConj
   (const int N, const float *X, const int incX, float *Y, const int incY);
void ATL_cmoveConj
   (const int N, const float *alpha, const float *X, const int incX,
    float *Y, const int incY);
void ATL_caxpyConj
   (const int N, const float *alpha, const float *X, const int incX,
    float *Y, const int incY);
void ATL_caxpyConj_x1_y1(const int N, const float *alpha, const float *X,
                         const int incX, float *Y, const int incY);
void ATL_caxpbyConj
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_caxpbyConj_aXi0_bXi0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a1_bXi0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a0_bXi0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aX_bXi0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aXi0_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a1_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a0_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aX_b1
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aXi0_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a1_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a0_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aX_b0
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aXi0_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a1_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_a0_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_caxpbyConj_aX_bX
   (ATL_CINT  N, const float *alpha, const float *X, ATL_CINT  incX,
    const float *beta, float *Y, ATL_CINT  incY);
void ATL_cgemoveC(const int N, const int M, const float *alpha,
                  const float *A, const int lda, float *C, const int ldc);

void ATL_cgeaddConj_aXi0_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a1_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a0_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aX_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a1_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a0_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aX_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a1_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a0_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aX_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a1_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_a0_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aXi0_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_cgeaddConj_aX_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_aXi0_b1
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_caxpby_aXi0_b1
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_aXi0_b1
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_aXi0_b0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_caxpby_aXi0_b0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_aXi0_b0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_aXi0_bXi0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_caxpby_aXi0_bXi0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_aXi0_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_aXi0_bX
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_caxpby_aXi0_bX
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_aXi0_bX
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_a1_bXi0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_a1_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_a0_bXi0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_a0_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);
void ATL_caxpby_aX_bXi0
   (const int N, const float *alpha, const float *X, const int incX,
    const float *beta, float *Y, const int incY);
void ATL_cgeadd_aX_bXi0
   (const int M, const int N, const float *alpha, const float *A,
    const int lda, const float *beta, float *C, const int ldc);

void ATL_cgemove_aXi0
   (const int M, const int N, const float *alpha0, const float *A,
    const int lda, float *C, const int ldc);

void ATL_cgescal_bXi0
   (const int M, const int N, const float *beta, float *C, const int ldc);

void ATL_zhecpyUNB_aX(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyUNB_ar(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyUNB_an(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyUNB_a1(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyLNB_aX(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyLNB_ar(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyLNB_an(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zhecpyLNB_a1(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyUNB_aX(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyUNB_ar(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyUNB_an(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyUNB_a1(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyLNB_aX(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyLNB_ar(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyLNB_an(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zsycpyLNB_a1(ATL_CSZT N, const double *alpha, const double *A,
                      ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_zgeput1H_bX(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zgeput1H_bN(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zgeput1H_b1(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zgeput1H_b0(ATL_CSZT M, ATL_CSZT N, const double *A, ATL_CSZT lda,
                     const double *beta, double *C, ATL_CSZT ldc);
void ATL_zswap_cplx2real(ATL_CUINT N, double *C, ATL_CINT incc,
                         double *R, ATL_CINT incR, double *I, ATL_CINT incI);
void ATL_zheApAc_NB
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *beta, double *C, ATL_CINT ldc);
void ATL_zheApAc
   (const enum ATLAS_UPLO Uplo, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *beta, double *C, ATL_CINT ldc);
void ATL_zgeApBc_NB
   (ATL_CINT M, ATL_CINT N, const double *A, ATL_CINT lda,
    const double *B, ATL_CINT ldb, const double *beta, double *C, ATL_CINT ldc);
void ATL_zcplxdivide(ATL_CINT N, const double *b, double *X, ATL_CINT incX,
                     double *Y, ATL_CINT incY);
void ATL_zcplxinvert
   (const int N, double *X, const int incX, double *Y, const int incY);

void ATL_zhereflect(const enum ATLAS_UPLO Uplo, const int N,
                    double *C, const int ldc);
void ATL_zscalConj
   (const int N, const double *alpha, double *X, const int incX);
void ATL_zcopyConj
   (const int N, const double *X, const int incX, double *Y, const int incY);
void ATL_zmoveConj
   (const int N, const double *alpha, const double *X, const int incX,
    double *Y, const int incY);
void ATL_zaxpyConj
   (const int N, const double *alpha, const double *X, const int incX,
    double *Y, const int incY);
void ATL_zaxpyConj_x1_y1(const int N, const double *alpha, const double *X,
                         const int incX, double *Y, const int incY);
void ATL_zaxpbyConj
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zaxpbyConj_aXi0_bXi0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a1_bXi0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a0_bXi0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aX_bXi0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aXi0_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a1_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a0_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aX_b1
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aXi0_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a1_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a0_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aX_b0
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aXi0_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a1_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_a0_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zaxpbyConj_aX_bX
   (ATL_CINT  N, const double *alpha, const double *X, ATL_CINT  incX,
    const double *beta, double *Y, ATL_CINT  incY);
void ATL_zgemoveC(const int N, const int M, const double *alpha,
                  const double *A, const int lda, double *C, const int ldc);

void ATL_zgeaddConj_aXi0_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a1_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a0_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aX_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a1_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a0_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aX_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a1_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a0_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aX_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a1_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_a0_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aXi0_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zgeaddConj_aX_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_aXi0_b1
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zaxpby_aXi0_b1
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_aXi0_b1
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_aXi0_b0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zaxpby_aXi0_b0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_aXi0_b0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_aXi0_bXi0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zaxpby_aXi0_bXi0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_aXi0_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_aXi0_bX
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zaxpby_aXi0_bX
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_aXi0_bX
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_a1_bXi0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_a1_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_a0_bXi0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_a0_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);
void ATL_zaxpby_aX_bXi0
   (const int N, const double *alpha, const double *X, const int incX,
    const double *beta, double *Y, const int incY);
void ATL_zgeadd_aX_bXi0
   (const int M, const int N, const double *alpha, const double *A,
    const int lda, const double *beta, double *C, const int ldc);

void ATL_zgemove_aXi0
   (const int M, const int N, const double *alpha0, const double *A,
    const int lda, double *C, const int ldc);

void ATL_zgescal_bXi0
   (const int M, const int N, const double *beta, double *C, const int ldc);


void ATL_qdgecollapse(const int M, const int N, ATL_QTYPE *C,
                      const int dldc, const int sldc);
void ATL_dsgecollapse(const int M, const int N, double *C,
                      const int dldc, const int sldc);
void ATL_ezgecollapse(const int M, const int N, ATL_QTYPE *C,
                      const int dldc, const int sldc);
void ATL_zcgecollapse(const int M, const int N, double *C,
                      const int dldc, const int sldc);
void ATL_qdtrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, ATL_QTYPE *C, const int dldc,const int sldc);
void ATL_dstrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, double *C, const int dldc, const int sldc);
void ATL_zctrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, double *C, const int dldc, const int sldc);
void ATL_eztrcollapse(const enum ATLAS_UPLO Uplo, const enum ATLAS_DIAG Diag,
                      const int N, ATL_QTYPE *C, const int dldc,const int sldc);

/*
 * This is the general LRU-based flush
 */
#if defined(ATL_USEPTHREADS) && !defined(ATL_flushcache)
   #include "atlas_pthreads.h"
   #define ATL_flushcache ATL_ptflushcache
   #define ATL_flushCacheByAddr ATL_ptflushCacheByAddr
   #define ATL_PTCACHEMUL * ATL_NTHREADS
#else
   #define ATL_PTCACHEMUL
#endif
double ATL_flushcache(long long size);
#ifdef ATL_USEPTHREADS
   void ATL_NumaTouchSpread(size_t N, void *buff);
#endif
/*
 * If we have it, use assembly-based explicit cache-line flush algorithm
 */
#if defined(ATL_ARCH_PPCG5) || defined(ATL_ARCH_PPCG4) || \
    defined(ATL_GAS_PPC) || defined(ATL_SSE2) || \
    defined(ATL_ARCH_IA64Itan) || defined(ATL_ARCH_IA64Itan2) || \
    defined(ATL_ARM64)

   #define ATL_LINEFLUSH 1
   typedef struct flStruct FLSTRUCT;
   struct flStruct
   {
      void *p;
      size_t length;
      FLSTRUCT *next;
   };
   FLSTRUCT *ATL_GetFlushStruct(void *p, size_t length, FLSTRUCT *next);
   void ATL_KillAllFlushStructs(FLSTRUCT *p);
   void ATL_flushCacheByAddr(size_t N, void *vp);
   void ATL_FlushAreasByCL(FLSTRUCT *fp);
   #if defined(ATL_USEPTHREADS) && !defined(ATL_FlushAreasByCL)
       void ATL_ptFlushAreasByCL(FLSTRUCT *fp);
       #define ATL_FlushAreasByCL ATL_ptFlushAreasByCL
   #endif

#else
    #define ATL_LINEFLUSH 0
#endif

#endif
