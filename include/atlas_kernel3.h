/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#ifndef ATLAS_KERNEL3_H
#define ATLAS_KERNEL3_H

/*
 * Real level 3 kernels
 */
void ATL_ssymmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_ssymmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_ssymmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_ssymmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_strsmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strsmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_strmmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
int ATL_ssyrkLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyrkUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyrkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyrkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyr2kLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyr2kUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyr2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_ssyr2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
void ATL_dsymmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_dsymmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_dsymmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_dsymmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_dtrsmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrsmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_dtrmmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
int ATL_dsyrkLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyrkUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyrkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyrkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyr2kLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyr2kUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyr2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_dsyr2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);

/*
 * Complex level 3 kernels
 */
void ATL_chemmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_chemmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_chemmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_chemmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_csymmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_csymmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_csymmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_csymmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_ctrsmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmLUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmLUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrsmRUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ctrmmRUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
int ATL_cherkLC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cherkUC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cherkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cherkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyrkLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyrkUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyrkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyrkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cher2kLC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cher2kUC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cher2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_cher2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyr2kLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyr2kUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyr2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_csyr2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
void ATL_zhemmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zhemmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zhemmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zhemmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zsymmRU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zsymmLU
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zsymmRL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_zsymmLL
   (ATL_CSZT  M, ATL_CSZT  N, const void *alpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *beta, void *C, ATL_CSZT ldc);
void ATL_ztrsmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmLUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmLUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRLCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUCN
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrsmRUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
void ATL_ztrmmRUCU
   (ATL_CSZT  M, ATL_CSZT  N, const void *valpha, const void *A, ATL_CSZT lda,
    void *C, ATL_CSZT ldc);
int ATL_zherkLC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zherkUC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zherkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zherkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyrkLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyrkUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyrkLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyrkUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zher2kLC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zher2kUC
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zher2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zher2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyr2kLT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyr2kUT
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyr2kLN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);
int ATL_zsyr2kUN
   (ATL_CSZT  N, ATL_CSZT  K, const void *valpha, const void *A, ATL_CSZT lda,
    const void *B, ATL_CSZT ldb, const void *vbeta, void *C, ATL_CSZT ldc);

/*
 * Real level 3 kernel auxiliaries
 */
void ATL_ssycopyU_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_ssycopyL_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_N_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_U_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_N_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_U_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_N_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_U_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_N_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_U_a0
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_ssycopyU_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_ssycopyL_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_N_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_U_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_N_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_U_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_N_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_U_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_N_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_U_a1
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_ssycopyU_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_ssycopyL_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_N_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2L_U_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_N_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyU2U_U_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_N_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2L_U_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_N_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strcopyL2U_U_aX
   (ATL_CSZT  N, const float alpha, const float *A, ATL_CSZT lda, float *C);
void ATL_strinvertUU(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_strinvertLU(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_strinvertUN(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_strinvertLN(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putU_bX
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putL_bX
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputU_bX
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputL_bX
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putU_b1
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putL_b1
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputU_b1
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputL_b1
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putU_b0
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_ssyr2k_putL_b0
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputU_b0
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strputL_b0
   (ATL_CSZT  N, const float *v, const float beta, float *A, ATL_CSZT lda);
void ATL_strsmKLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_strsmKRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const float alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_dsycopyU_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dsycopyL_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_N_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_U_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_N_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_U_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_N_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_U_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_N_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_U_a0
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dsycopyU_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dsycopyL_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_N_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_U_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_N_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_U_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_N_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_U_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_N_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_U_a1
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dsycopyU_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dsycopyL_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_N_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2L_U_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_N_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyU2U_U_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_N_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2L_U_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_N_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrcopyL2U_U_aX
   (ATL_CSZT  N, const double alpha, const double *A, ATL_CSZT lda, double *C);
void ATL_dtrinvertUU(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_dtrinvertLU(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_dtrinvertUN(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_dtrinvertLN(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putU_bX
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putL_bX
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputU_bX
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputL_bX
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putU_b1
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putL_b1
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputU_b1
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputL_b1
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putU_b0
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dsyr2k_putL_b0
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputU_b0
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrputL_b0
   (ATL_CSZT  N, const double *v, const double beta, double *A, ATL_CSZT lda);
void ATL_dtrsmKLLTN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLLTU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLLNN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLLNU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLUTN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLUTU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLUNN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKLUNU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRLTN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRLTU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRLNN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRLNU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRUTN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRUTU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRUNN
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_dtrsmKRUNU
   (ATL_CSZT  M, ATL_CSZT  N, const double alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);

/*
 * Complex level 3 kernel auxiliaries
 */
void ATL_cCtrsmKL
   (enum ATLAS_UPLO Uplo, enum ATLAS_TRANS Trans, enum ATLAS_DIAG Diag,
    ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *B, ATL_CSZT ldb);
void ATL_checopy
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_csycopy
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2L_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2Lc_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2L_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2Lc_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2U_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2Uc_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2U_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyU2Uc_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2L_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2Lc_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2L_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2Lc_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2U_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2Uc_N
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2U_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrcopyL2Uc_U
   (ATL_CSZT  N, const float *A, ATL_CSZT lda, float *C);
void ATL_ctrmv_scalLNU_an1
   (ATL_CSZT  N, const float *alpha, const float *A, ATL_CSZT lda, float *X);
void ATL_ctrmv_scalLNN_aX
   (ATL_CSZT  N, const float *alpha, const float *A, ATL_CSZT lda, float *X);
void ATL_ctrmv_scalUNU_an1
   (ATL_CSZT  N, const float *alpha, const float *A, ATL_CSZT lda, float *X);
void ATL_ctrmv_scalUNN_aX
   (ATL_CSZT  N, const float *alpha, const float *A, ATL_CSZT lda, float *X);
void ATL_ctrinvertUU(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_ctrinvertLU(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_ctrinvertUN(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_ctrinvertLN(ATL_CSZT  N, float *A, ATL_CSZT lda);
void ATL_ctrputU_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputL_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putU_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putL_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputU_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputL_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putU_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putL_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputU_bX
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputL_bX
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putU_bX
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putL_bX
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputU_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputL_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putU_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putL_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputU_bn1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrputL_bn1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putU_bn1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_csyr2k_putL_bn1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putU_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putL_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputU_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputL_b0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putU_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putL_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputU_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputL_b1
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putU_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cher2k_putL_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputU_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_cheputL_bXi0
   (ATL_CSZT  N, const float *v, const float *beta, float *A, ATL_CSZT lda);
void ATL_ctrsm0LLTN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LLTU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LLNN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LLNU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LLCN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LLCU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUTN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUTU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUNN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUNU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUCN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0LUCU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLTN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLTU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLNN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLNU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLCN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RLCU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUTN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUTU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUNN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUNU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUCN
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_ctrsm0RUCU
   (ATL_CSZT  M, ATL_CSZT  N, const float *alpha, const float *A,
    ATL_CSZT lda, float *C, ATL_CSZT ldc);
void ATL_zCtrsmKL
   (enum ATLAS_UPLO Uplo, enum ATLAS_TRANS Trans, enum ATLAS_DIAG Diag,
    ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *B, ATL_CSZT ldb);
void ATL_zhecopy
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_zsycopy
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2L_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2Lc_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2L_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2Lc_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2U_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2Uc_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2U_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyU2Uc_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2L_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2Lc_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2L_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2Lc_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2U_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2Uc_N
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2U_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrcopyL2Uc_U
   (ATL_CSZT  N, const double *A, ATL_CSZT lda, double *C);
void ATL_ztrmv_scalLNU_an1
   (ATL_CSZT  N, const double *alpha, const double *A, ATL_CSZT lda, double *X);
void ATL_ztrmv_scalLNN_aX
   (ATL_CSZT  N, const double *alpha, const double *A, ATL_CSZT lda, double *X);
void ATL_ztrmv_scalUNU_an1
   (ATL_CSZT  N, const double *alpha, const double *A, ATL_CSZT lda, double *X);
void ATL_ztrmv_scalUNN_aX
   (ATL_CSZT  N, const double *alpha, const double *A, ATL_CSZT lda, double *X);
void ATL_ztrinvertUU(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_ztrinvertLU(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_ztrinvertUN(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_ztrinvertLN(ATL_CSZT  N, double *A, ATL_CSZT lda);
void ATL_ztrputU_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputL_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putU_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putL_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputU_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputL_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putU_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putL_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputU_bX
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputL_bX
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putU_bX
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putL_bX
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputU_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputL_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putU_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putL_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputU_bn1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrputL_bn1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putU_bn1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zsyr2k_putL_bn1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putU_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putL_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputU_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputL_b0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putU_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putL_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputU_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputL_b1
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putU_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zher2k_putL_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputU_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_zheputL_bXi0
   (ATL_CSZT  N, const double *v, const double *beta, double *A, ATL_CSZT lda);
void ATL_ztrsm0LLTN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LLTU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LLNN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LLNU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LLCN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LLCU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUTN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUTU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUNN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUNU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUCN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0LUCU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLTN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLTU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLNN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLNU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLCN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RLCU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUTN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUTU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUNN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUNU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUCN
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);
void ATL_ztrsm0RUCU
   (ATL_CSZT  M, ATL_CSZT  N, const double *alpha, const double *A,
    ATL_CSZT lda, double *C, ATL_CSZT ldc);

#endif
