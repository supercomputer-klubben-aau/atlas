/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */
#ifndef ATLAS_LEVEL3_H
   #define ATLAS_LEVEL3_H

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
int ATL_sGetNB(void);
int ATL_sGetNCNB(void);
void ATL_sammm(enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
               ATL_CSZT, const float , const float*, ATL_CSZT,
               const float*,ATL_CSZT, const float , float*,ATL_CSZT);
void ATL_strsm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT M, ATL_CSZT N, const float alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT ldb);
void ATL_sgemm(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
               ATL_CSZT  M, ATL_CSZT  N, ATL_CSZT  K, const float alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float beta, float *C, ATL_CSZT  ldc);
void ATL_ssymm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float beta, float *C, ATL_CSZT  ldc);
void ATL_ssymm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float beta, float *C, ATL_CSZT  ldc);
void ATL_ssyrk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const float alp, const float *A,
               ATL_CSZT lda, const float bet, float *C, ATL_CSZT ldc);
void ATL_ssyrk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT  N, ATL_CSZT  K, const float alpha,
               const float *A, ATL_CSZT  lda, const float beta,
               float *C, ATL_CSZT  ldc);
void ATL_ssyr2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
                ATL_CSZT  N, ATL_CSZT  K, const float alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float beta, float *C, ATL_CSZT  ldc);
void ATL_ssyr2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const float alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float beta, float *C, ATL_CSZT  ldc);
void ATL_strmm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const float alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT  ldb);
void ATL_strmm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const float alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT  ldb);
int ATL_strsmR_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_strsmL_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_siptrsmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_siptrsmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
void ATL_strsm(const enum ATLAS_SIDE, const enum ATLAS_UPLO,
               const enum ATLAS_TRANS, const enum ATLAS_DIAG, ATL_CSZT,
               ATL_CSZT, const float , const float*, ATL_CSZT, float*,
               ATL_CSZT);
int ATL_strmmR_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_strmmL_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_siptrmmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_siptrmmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float ,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_dGetNB(void);
int ATL_dGetNCNB(void);
void ATL_dammm(enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
               ATL_CSZT, const double , const double*, ATL_CSZT,
               const double*,ATL_CSZT, const double , double*,ATL_CSZT);
void ATL_dtrsm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT M, ATL_CSZT N, const double alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT ldb);
void ATL_dgemm(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
               ATL_CSZT  M, ATL_CSZT  N, ATL_CSZT  K, const double alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double beta, double *C, ATL_CSZT  ldc);
void ATL_dsymm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double beta, double *C, ATL_CSZT  ldc);
void ATL_dsymm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double beta, double *C, ATL_CSZT  ldc);
void ATL_dsyrk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const double alp, const double *A,
               ATL_CSZT lda, const double bet, double *C, ATL_CSZT ldc);
void ATL_dsyrk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT  N, ATL_CSZT  K, const double alpha,
               const double *A, ATL_CSZT  lda, const double beta,
               double *C, ATL_CSZT  ldc);
void ATL_dsyr2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
                ATL_CSZT  N, ATL_CSZT  K, const double alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double beta, double *C, ATL_CSZT  ldc);
void ATL_dsyr2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const double alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double beta, double *C, ATL_CSZT  ldc);
void ATL_dtrmm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const double alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT  ldb);
void ATL_dtrmm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const double alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT  ldb);
int ATL_dtrsmR_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_dtrsmL_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_diptrsmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_diptrsmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
void ATL_dtrsm(const enum ATLAS_SIDE, const enum ATLAS_UPLO,
               const enum ATLAS_TRANS, const enum ATLAS_DIAG, ATL_CSZT,
               ATL_CSZT, const double , const double*, ATL_CSZT, double*,
               ATL_CSZT);
int ATL_dtrmmR_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_dtrmmL_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_diptrmmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_diptrmmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double ,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_cGetNB(void);
int ATL_cGetNCNB(void);
void ATL_cammm(enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
               ATL_CSZT, const float *, const float*, ATL_CSZT,
               const float*,ATL_CSZT, const float *, float*,ATL_CSZT);
void ATL_ctrsm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT M, ATL_CSZT N, const float *alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT ldb);
void ATL_cgemm(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
               ATL_CSZT  M, ATL_CSZT  N, ATL_CSZT  K, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float *beta, float *C, ATL_CSZT  ldc);
void ATL_csymm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float *beta, float *C, ATL_CSZT  ldc);
void ATL_csymm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float *beta, float *C, ATL_CSZT  ldc);
void ATL_csyrk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const float *alp, const float *A,
               ATL_CSZT lda, const float *bet, float *C, ATL_CSZT ldc);
void ATL_csyrk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT  N, ATL_CSZT  K, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *beta,
               float *C, ATL_CSZT  ldc);
void ATL_csyr2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
                ATL_CSZT  N, ATL_CSZT  K, const float *alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float *beta, float *C, ATL_CSZT  ldc);
void ATL_csyr2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const float *alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float *beta, float *C, ATL_CSZT  ldc);
void ATL_ctrmm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT  ldb);
void ATL_ctrmm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, float *B, ATL_CSZT  ldb);
int ATL_ctrsmR_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ctrsmL_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ciptrsmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ciptrsmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
void ATL_ctrsm(const enum ATLAS_SIDE, const enum ATLAS_UPLO,
               const enum ATLAS_TRANS, const enum ATLAS_DIAG, ATL_CSZT,
               ATL_CSZT, const float *, const float*, ATL_CSZT, float*,
               ATL_CSZT);
int ATL_ctrmmR_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ctrmmL_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ciptrmmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_ciptrmmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const float *,
               const float*, ATL_CSZT, float*, ATL_CSZT);
int ATL_zGetNB(void);
int ATL_zGetNCNB(void);
void ATL_zammm(enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
               ATL_CSZT, const double *, const double*, ATL_CSZT,
               const double*,ATL_CSZT, const double *, double*,ATL_CSZT);
void ATL_ztrsm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT M, ATL_CSZT N, const double *alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT ldb);
void ATL_zgemm(const enum ATLAS_TRANS, const enum ATLAS_TRANS,
               ATL_CSZT  M, ATL_CSZT  N, ATL_CSZT  K, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double *beta, double *C, ATL_CSZT  ldc);
void ATL_zsymm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double *beta, double *C, ATL_CSZT  ldc);
void ATL_zsymm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double *beta, double *C, ATL_CSZT  ldc);
void ATL_zsyrk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const double *alp, const double *A,
               ATL_CSZT lda, const double *bet, double *C, ATL_CSZT ldc);
void ATL_zsyrk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT  N, ATL_CSZT  K, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *beta,
               double *C, ATL_CSZT  ldc);
void ATL_zsyr2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
                ATL_CSZT  N, ATL_CSZT  K, const double *alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double *beta, double *C, ATL_CSZT  ldc);
void ATL_zsyr2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const double *alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double *beta, double *C, ATL_CSZT  ldc);
void ATL_ztrmm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT  ldb);
void ATL_ztrmm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               const enum ATLAS_TRANS TransA, const enum ATLAS_DIAG Diag,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, double *B, ATL_CSZT  ldb);
int ATL_ztrsmR_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ztrsmL_IP(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ziptrsmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ziptrsmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
void ATL_ztrsm(const enum ATLAS_SIDE, const enum ATLAS_UPLO,
               const enum ATLAS_TRANS, const enum ATLAS_DIAG, ATL_CSZT,
               ATL_CSZT, const double *, const double*, ATL_CSZT, double*,
               ATL_CSZT);
int ATL_ztrmmR_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ztrmmL_ST(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ziptrmmR(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);
int ATL_ziptrmmL(ATL_UINT bv, ATL_CSZT, ATL_CSZT, const double *,
               const double*, ATL_CSZT, double*, ATL_CSZT);

/*
 * Routines with prefixes C and Z only
 */
void ATL_cherk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const float alp, const float *A,
               ATL_CSZT lda, const float bet, float *C, ATL_CSZT ldc);
void ATL_cher2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT N, ATL_CSZT K, const float *alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float beta, float *C, ATL_CSZT  ldc);
void ATL_chemm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float *beta, float *C, ATL_CSZT  ldc);
void ATL_cherk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT N, ATL_CSZT K, const float alpha,
               const float *A, ATL_CSZT  lda, const float beta,
               float *C, ATL_CSZT ldc);
void ATL_cher2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const float *alpha,
                const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
                const float beta, float *C, ATL_CSZT  ldc);
void ATL_chemm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const float *alpha,
               const float *A, ATL_CSZT  lda, const float *B, ATL_CSZT  ldb,
               const float *beta, float *C, ATL_CSZT  ldc);

void ATL_zherk(const enum ATLAS_UPLO, const enum ATLAS_TRANS,
               ATL_CSZT N, ATL_CSZT K, const double alp, const double *A,
               ATL_CSZT lda, const double bet, double *C, ATL_CSZT ldc);
void ATL_zher2k(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT N, ATL_CSZT K, const double *alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double beta, double *C, ATL_CSZT  ldc);
void ATL_zhemm(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double *beta, double *C, ATL_CSZT  ldc);
void ATL_zherk_APR(const enum ATLAS_UPLO Uplo,const enum ATLAS_TRANS Trans,
               ATL_CSZT N, ATL_CSZT K, const double alpha,
               const double *A, ATL_CSZT  lda, const double beta,
               double *C, ATL_CSZT ldc);
void ATL_zher2k_APR(const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
                ATL_CSZT  N, ATL_CSZT  K, const double *alpha,
                const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
                const double beta, double *C, ATL_CSZT  ldc);
void ATL_zhemm_APR(const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
               ATL_CSZT  M, ATL_CSZT  N, const double *alpha,
               const double *A, ATL_CSZT  lda, const double *B, ATL_CSZT  ldb,
               const double *beta, double *C, ATL_CSZT  ldc);


#endif
