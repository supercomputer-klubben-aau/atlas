/*
 * Prototypes trsm kernel routines used in gemm or amm-based BLAS routines
 */
#ifndef ATLAS_TRSMK_H
   #define ATLAS_TRSMK_H 1
#include "atlas_misc.h"
#include "atlas_amm.h"


int ATL_strsmKL_rk4
   (enum ATLAS_SIDE Side, enum ATLAS_UPLO Uplo, enum ATLAS_TRANS TA,
    enum ATLAS_DIAG Diag, ATL_CINT M, ATL_CINT N, const float alpha,
    const float *A, ATL_CINT lda, float *B, ATL_CINT ldb);
void ATL_sktrsmLLN_rk4  /* rank-4 based trsm kernel */
(
   ATL_CINT  M,   /* size of orig triangular matrix A */
   ATL_CINT  N,   /* number of RHS in B */
   const float alpha,  /* scale factor for B */
   const float *A, /* MxM lower matrix A, diag has inverse of original diag */
   TYPE *B,       /* on input, B, on output X, of A x = b */
   ATL_CINT  ldb, /* leading dim of B */
   float *W        /* Mx4 workspace with good alignment */
);

int ATL_strsmKL_amm  /* amm-based trsm kernel primitive */
(
   amminfo_t *mmp,      /* chosen amm kernel info */
   enum ATLAS_DIAG Diag,
   ATL_CUINT mb0,       /* sizeof first triangular block */
   ATL_CUINT nmblks,    /* CEIL(M/mb); M = (nmblks-1)*b+mb0 */
   ATL_CUINT nnblks,    /* FLOOR(N/nb) */
   ATL_CUINT nr,        /* mod(N/nb) */
   const SCALAR alpha,  /* scale factor for rhs */
   const float *L,       /* L, already copied to row-panel storage as above */
   size_t incL,         /* distance between blocks in L */
   float *R,             /* ptr to col-major rhs (B on input, X on output */
   size_t ldr,          /* leading dim of RHS matrix R */
   ATL_CUINT nmu,       /* nmu = mb/mu, assume mod(mb/mu) == 0 */
   ATL_CUINT nnu,       /* nnu = NB/nu, assume mod(NB/nu) == 0 */
   TYPE *c,             /* bxNB with trailing readable space for amm */
   TYPE *W,             /* aligned workspace for amm's B blocks */
   size_t incW          /* stride between b*NB blocks */
);

int ATL_strsmL_amm      /* general Left amm-based trsm */
(                       /* RETURNS: 0 on success, non-0 if op not done */
   enum ATLAS_SIDE Side,
   enum ATLAS_UPLO Uplo,
   enum ATLAS_TRANS TA,
   enum ATLAS_DIAG Diag,
   ATL_CINT  M,         /* size of triangular matrix A */
   ATL_CINT  N,         /* number of RHS in B */
   const SCALAR alpha,  /* scale factor for B */
   const float *T,     /* MxM triangular matrix A */
   ATL_CINT  ldt,       /* leading dim of T */
   float *R,           /* on input, B, on output X, of A x = b */
   ATL_CINT  ldr        /* leading dim of R */
);

void Mjoin(PATL,trsmK_L2blk)   /* takes data in col-maj Lower format */
(                              /* puts in format required by trsmL_amm */
   enum ATLAS_DIAG Diag,
   const int mb0,   /* size of first triangular block */
   const int b,     /* size of all other blocks */
   int nmblks,      /* number of blocks of size b (not incl mb0) */
   const TYPE *L,   /* lower triang mat of x = inv(L) b */
   size_t ldl,      /* leading dim of L */
   TYPE *W,         /* wrkspc to copy blocks to */
   size_t incW,     /* stride between blocks (incW >= MAX(b,mb0)^2) */
   cm2am_t a2blk    /* copy rect blks to amm storage */
);

int ATL_dtrsmKL_rk4
   (enum ATLAS_SIDE Side, enum ATLAS_UPLO Uplo, enum ATLAS_TRANS TA,
    enum ATLAS_DIAG Diag, ATL_CINT M, ATL_CINT N, const double alpha,
    const double *A, ATL_CINT lda, double *B, ATL_CINT ldb);
void ATL_dktrsmLLN_rk4  /* rank-4 based trsm kernel */
(
   ATL_CINT  M,   /* size of orig triangular matrix A */
   ATL_CINT  N,   /* number of RHS in B */
   const double alpha,  /* scale factor for B */
   const double *A, /* MxM lower matrix A, diag has inverse of original diag */
   TYPE *B,       /* on input, B, on output X, of A x = b */
   ATL_CINT  ldb, /* leading dim of B */
   double *W        /* Mx4 workspace with good alignment */
);

int ATL_dtrsmKL_amm  /* amm-based trsm kernel primitive */
(
   amminfo_t *mmp,      /* chosen amm kernel info */
   enum ATLAS_DIAG Diag,
   ATL_CUINT mb0,       /* sizeof first triangular block */
   ATL_CUINT nmblks,    /* CEIL(M/mb); M = (nmblks-1)*b+mb0 */
   ATL_CUINT nnblks,    /* FLOOR(N/nb) */
   ATL_CUINT nr,        /* mod(N/nb) */
   const SCALAR alpha,  /* scale factor for rhs */
   const double *L,       /* L, already copied to row-panel storage as above */
   size_t incL,         /* distance between blocks in L */
   double *R,             /* ptr to col-major rhs (B on input, X on output */
   size_t ldr,          /* leading dim of RHS matrix R */
   ATL_CUINT nmu,       /* nmu = mb/mu, assume mod(mb/mu) == 0 */
   ATL_CUINT nnu,       /* nnu = NB/nu, assume mod(NB/nu) == 0 */
   TYPE *c,             /* bxNB with trailing readable space for amm */
   TYPE *W,             /* aligned workspace for amm's B blocks */
   size_t incW          /* stride between b*NB blocks */
);

int ATL_dtrsmL_amm      /* general Left amm-based trsm */
(                       /* RETURNS: 0 on success, non-0 if op not done */
   enum ATLAS_SIDE Side,
   enum ATLAS_UPLO Uplo,
   enum ATLAS_TRANS TA,
   enum ATLAS_DIAG Diag,
   ATL_CINT  M,         /* size of triangular matrix A */
   ATL_CINT  N,         /* number of RHS in B */
   const SCALAR alpha,  /* scale factor for B */
   const double *T,     /* MxM triangular matrix A */
   ATL_CINT  ldt,       /* leading dim of T */
   double *R,           /* on input, B, on output X, of A x = b */
   ATL_CINT  ldr        /* leading dim of R */
);

void Mjoin(PATL,trsmK_L2blk)   /* takes data in col-maj Lower format */
(                              /* puts in format required by trsmL_amm */
   enum ATLAS_DIAG Diag,
   const int mb0,   /* size of first triangular block */
   const int b,     /* size of all other blocks */
   int nmblks,      /* number of blocks of size b (not incl mb0) */
   const TYPE *L,   /* lower triang mat of x = inv(L) b */
   size_t ldl,      /* leading dim of L */
   TYPE *W,         /* wrkspc to copy blocks to */
   size_t incW,     /* stride between blocks (incW >= MAX(b,mb0)^2) */
   cm2am_t a2blk    /* copy rect blks to amm storage */
);
#endif
