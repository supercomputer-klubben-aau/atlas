#ifndef ATLAS_AMM_H
   #define ATLAS_AMM_H

#ifndef ATL_MaxMalloc
   #include "atlas_maxmalloc.h"
#endif
#ifndef ATL_MaxMalloc
   #define ATL_MaxMalloc 268435456UL
#endif
#include "atlas_misc.h"

#define ATL_AMK_NBIT 5
#define ATL_AMK_KVEC  1
#define ATL_AMK_KRUN  2
#define ATL_AMK_KMIN  4
#define ATL_AMK_KMAX  8
#define ATL_AMK_KALL 16

#ifdef TREAL
   typedef void (*cm2am_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, const size_t, TYPE*);
   typedef void (*tcm2am_t)(const size_t, const SCALAR, const TYPE*,
                           const size_t, TYPE*);
   typedef void (*am2cm_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, TYPE*, const size_t);
   typedef void (*ablk2cmat_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const SCALAR, TYPE *, const size_t);
   typedef void (*cmat2ablk_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const size_t, const SCALAR,TYPE*);
   typedef void (*ammswp_t)(ATL_CINT, TYPE*,ATL_CSZT,TYPE*);
#else
   typedef void (*cm2am_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, const size_t, TYPE*, TYPE*);
   typedef void (*tcm2am_t)(const size_t, const SCALAR, const TYPE*,
                           const size_t, TYPE*, TYPE*);
   typedef void (*am2cm_t)(const size_t, const size_t, const SCALAR,
                           const TYPE*, const TYPE*, TYPE*, const size_t);
   typedef void (*ablk2cmat_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const TYPE*, const SCALAR,
                               TYPE *, const size_t);
   typedef void (*cmat2ablk_t)(const size_t, const size_t, const SCALAR,
                               const TYPE*, const size_t, const SCALAR,
                               TYPE*,TYPE*);
   typedef void (*ammswp_t)(ATL_CINT, TYPE*,ATL_CSZT,TYPE*,TYPE*);
#endif
typedef void (*ammkern_t)(ATL_CSZT,ATL_CSZT,ATL_CSZT,const TYPE*,const TYPE*,
                          TYPE*, const TYPE*, const TYPE*, const TYPE*);
#define ushort unsigned short
#define uint unsigned int
#define uchar unsigned char

typedef struct sminfo sminfo_t;  /* struct for trsm microkernel */
struct sminfo
{
   void *utrsm;                  /* function pointer to trsm microkernel */
   cm2am_t a2blk, b2blk;         /* A/B copy funcs */
   ammkern_t amm_b0, amm_b1;     /* func ptrs for beta=0,1 gemm */
   #ifdef TCPLX
      ammkern_t amm_bn;          /* func ptr for beta=-1 gemm */
   #endif
   size_t incA;                  /* (1|lda)SHIFT */
   ushort mmflg;                 /* flag bitvec from amm master list */
   ushort kb, rb;                /* triangular & RHS blocking */
   uchar mu, nu, ku;             /* unrollings of underlying amm */
   uchar bv;/*0:Right?, 1:Upper?, 2:TransA?, 3:Conj?, 4:NonUnit? 5: ALLT? */
};
typedef struct tminfo tminfo_t;  /* struct for trmm microkernel */
struct tminfo
{
   ushort flg;                   /* 0 bit: same kernel as gemm? */
   tcm2am_t t2blk;               /* copy routine for triangular matrix */
   cm2am_t r2blk;                /* copy routines for regular matrix */
   ablk2cmat_t blk2c;            /* cpy from acc-maj blk to col-maj C */
   ammkern_t amm_b0, amm_b1;     /* func ptrs for beta=0,1 trmm */
   #ifdef TCPLX
      ammkern_t amm_bn;          /* func ptr for beta=-1 trmm */
   #endif
   size_t incA;                  /* (1|lda) */
   ushort kb;                    /* triangular blocking */
   ushort mu, nu, ku;            /* unrollings of underlying amm */
   ushort vlen;                  /* vector len for the kernels */
   ushort kvec;                  /* kvec kernel ? */
};

typedef struct amminfo amminfo_t;
struct amminfo
{
   cm2am_t a2blk, b2blk;
   ablk2cmat_t Cblk2cm, Cblk2cm_b1;
   cmat2ablk_t cm2Cblk;
   ammkern_t amm_b0, amm_b1, amm_bn, amm_k1_b0, amm_k1_b1, amm_k1_bn;
   ushort IDX, mb, nb, kb, kbmin, kbmax;
   uchar flag, mu, nu, ku, vlen;
};

typedef struct ipinfo ipinfo_t;  /* struct for inner-product based gemm */
struct ipinfo
{
   ammkern_t ammK1_b0, ammK1_b1; /* amm safe to use for KB0 calc */
   ammkern_t amm_b0, amm_b1;     /* amm to use after first K peel */
   cm2am_t a2blk, b2blk;         /* A/B copy funcs */
   ablk2cmat_t blk2c, blk2c_b1;  /* cpy from acc-maj blk to col-maj C */
   #ifdef TCPLX
      ammkern_t ammK1_bn;
      ammkern_t amm_bn;
      const TYPE *alpA, *alpB, *alpC, *ONE;
   #endif
   size_t nfnblks, npnblks; /* # of of full and partial blks along N */
   size_t nfmblks, npmblks; /* # of of full and partial blks along M */
   size_t nfkblks;          /* # of full kblks = (K-kb0)/kb */
   size_t incAk, incBk;     /* block-kb increment for A/B ptrs */
   size_t pincAm, incAm;    /* partial & full mb-block increment */
   size_t pincBn, incBn;    /* partial & full nb-block increment */
   size_t lda, ldb, ldc;    /* leading dim for all 3 matrices */
   #ifndef TCPLX
      TYPE alpA, alpB, alpC;
   #endif
   uint szC;                /* size of real portion of C workspace */
   uint pszA, szA;          /* size of real portion of part & full blks of A */
   uint pszB, szB;          /* size of real portion of part & full blks of B */
   uint idx;                /* index in view of this amm case */
   uint idxK;               /* index in atlas_?amm_kern.h */
   ushort exsz;             /* extra bytes to alloc past end */
   ushort mF, nF;           /* size of final part of mat before copy */
   ushort mb, pmb, nb, pnb; /* full & partial blk factors for M & N */
   ushort kb;               /* sz of all K blks except first */
   ushort kb0;              /* K%kb or if zero, kb */
   ushort KB0;              /* if amm k-vec, ((kb0+ku-1)/ku)*ku, else kb0 */
   ushort nmu, nnu;         /* # of unrollings along each dim */
   ushort pnmu, pnnu;       /* # of unrollings in partial blocks */
   ushort nmuF, nnuF;       /* # of unrolling for final (remainder) block */
   uchar mu, nu, ku;        /* unrolling for each dim used by amm */
   uchar vlen;              /* vector length (needed for KVEC kerns) */
};

typedef struct opinfo opinfo_t;
struct opinfo
{
   cm2am_t a2blk, b2blk;    /* A/B copy routs */
   ablk2cmat_t blk2C;       /* cpy from acc-maj blk to col-maj C */
   ammkern_t amm_b0;        /* beta=0 amm kern */
   #ifdef TCPLX
      ammkern_t amm_b1;     /* beta=1 amm kern (needed for complex) */
      ammkern_t amm_bn;     /* beta=n amm kern (needed for complex) */
      const TYPE *alpA, *alpB, *beta, *ONE;
   #endif
   size_t nfnblks, npnblks; /* # of of full and partial blks along N */
   size_t nfmblks, npmblks; /* # of of full and partial blks along M */
   size_t pincAm, incAm;    /* partial & full mb-block increment */
   size_t pincBn, incBn;    /* partial & full nb-block increment */
   size_t lda, ldb, ldc;    /* stride between row elts for each mat */
   #ifndef TCPLX
      TYPE alpA, alpB, beta;
   #endif
   uint szC;                /* size of real C blk in workspace */
   uint pszA, szA;          /* size of real portion of part & full blks of A */
   uint pszB, szB;          /* size of real portion of part & full blks of B */
   ushort exsz;             /* extra bytes to alloc past end */
   ushort idx;              /* kern array index */
   ushort mF, nF;           /* size of final part of mat before copy */
   ushort mb, pmb, nb, pnb; /* full & partial blk factors for M & N */
   ushort kb;               /* K from original problem */
   ushort KB;               /* ((kb+ku-1)/vl)*vl if kvec, else kb */
   ushort nmuF, nnuF;       /* # of unrollings in final M/N blocks */
   ushort nmu, nnu;         /* # of unrollings along each dim */
   ushort pnmu, pnnu;       /* # of unrollings in partial blocks */
   uchar mu, nu, ku;        /* unrolling for each dim used by amm */
   uchar vlen;              /* vector length (needed for KVEC kerns) */
};
#undef ushort
#undef uchar
#undef uint

enum ATL_AMMALG {ATL_amm1b, ATL_ammrkK, ATL_NMK};

void Mjoin(PATL,opinfo) /* populate opinfo_t using atlas_opgen_view K-3 idx */
   (opinfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);

void Mjoin(PATL,iptrsmInfo)
   (ipinfo_t *ip, enum ATLAS_SIDE SD, enum ATLAS_UPLO Uplo, enum ATLAS_TRANS TA,
    size_t N, size_t K, size_t lda, size_t ldb, const SCALAR alpha);
void Mjoin(PATL,ipgenInfo)
   (ipinfo_t *ip, int flg, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
void Mjoin(PATL,ipmenInfo)
   (ipinfo_t *ip, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
void Mjoin(PATL,ipmenUMInfo)
   (ipinfo_t *ip, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
void Mjoin(PATL,ipmekInfo)
   (ipinfo_t *ip, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
void Mjoin(PATL,ipnekInfo)
   (ipinfo_t *ip, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
double Mjoin(PATL,ipsyrkInfo) /* returns predicted time */
   (ipinfo_t*, int, enum ATLAS_TRANS,size_t,size_t,size_t,size_t,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,opsyrkInfo)
   (opinfo_t *out, int flag, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K,
    ATL_CSZT lda, ATL_CSZT ldc, const SCALAR alpha, const SCALAR beta);
ablk2cmat_t Mjoin(PATL,opsyr2kInfo)
   (opinfo_t *out, int flag, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K,
    ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
void Mjoin(PATL,ipnekInfoPop)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta,
    size_t nfmblks, size_t npmblks, ATL_UINT mb, ATL_UINT pmb,
    size_t nfnblks, size_t npnblks, ATL_UINT nb, ATL_UINT pnb);

void Mjoin(PATL,ipmekInfoPop)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta,
    size_t nfmblks, size_t npmblks, ATL_UINT mb, ATL_UINT pmb,
    size_t nfnblks, size_t npnblks, ATL_UINT nb, ATL_UINT pnb);

void Mjoin(PATL,ipmenInfoPop)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta,
    size_t nfmblks, size_t npmblks, ATL_UINT mb, ATL_UINT pmb,
    size_t nfnblks, size_t npnblks, ATL_UINT nb, ATL_UINT pnb);

void Mjoin(PATL,ipgenInfoPop)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta,
    size_t nfmblks, size_t npmblks, ATL_UINT mb, ATL_UINT pmb,
    size_t nfnblks, size_t npnblks, ATL_UINT nb, ATL_UINT pnb);

static int Mjoin(PATL,geGetAmmmIndx)(size_t M, size_t N, size_t K)
{
   ATL_assert(0);
}

static int Mjoin(PATL,tGetIPInfo_tMN)
   (ipinfo_t *ip, int P, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, const SCALAR alpha, size_t lda, size_t ldb,
    const SCALAR beta, size_t ldc)
{
   ATL_assert(0);
}
static void Mjoin(PATL,sqComputeIPInfo)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ATL_assert(0);
}
static void Mjoin(PATL,geComputeIPInfo)
   (ipinfo_t *ip, int idx, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    size_t M, size_t N, size_t K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ATL_assert(0);
}
void Mjoin(PATL,doAmmBlk)
   (ipinfo_t*,ATL_CUINT,ATL_iptr_t,ATL_iptr_t,ATL_iptr_t,
    const TYPE*, const TYPE*, TYPE*, TYPE*, TYPE*, TYPE*, TYPE*,
    const SCALAR, const ablk2cmat_t);

void Mjoin(PATL,iploopsNMK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsNMK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsMK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsNK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsMNK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
void Mjoin(PATL,iploopsMNK)
   (ipinfo_t *ip, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    const int MVS, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC,
    const SCALAR beta, const ablk2cmat_t blk2c);
static int Mjoin(PATL,rkGetAmmInfoInt)(char wh, int idx)
{
   ATL_assert(0);
}
static void *Mjoin(PATL,rkGetAmmInfoPtr)(int idx, int what, int alp, int bet)
{
   ATL_assert(0);
}
static int Mjoin(PATL,sqGetAmmInfoInt)(char wh, int idx)
{
   ATL_assert(0);
}
static void *Mjoin(PATL,sqGetAmmInfoPtr)(int idx, int what, int alp, int bet)
{
   ATL_assert(0);
}
static int Mjoin(PATL,geGetAmmInfoInt)(char wh, int idx)
{
   ATL_assert(0);
}
static void *Mjoin(PATL,geGetAmmInfoPtr)(int idx, int what, int alp, int bet)
{
   ATL_assert(0);
}
static int Mjoin(PATL,tGetRKInfo_tK)
   (opinfo_t*a0, int a1, enum ATLAS_TRANS a2, enum ATLAS_TRANS a3,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ATL_assert(0);
}
static int Mjoin(PATL,tGetRKInfo_tNK)
   (opinfo_t*a0, int a1, enum ATLAS_TRANS a2, enum ATLAS_TRANS a3,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ATL_assert(0);
}
int Mjoin(PATL,GetInfo_1b_opnk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_1b_opmk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_1b_opsq)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_1b_oprk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_opnk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_opmk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_opsq)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetInfo_oprk)
   (opinfo_t*, enum ATLAS_TRANS, enum ATLAS_TRANS, ATL_CSZT M, ATL_CSZT N,
    ATL_CSZT K, ATL_CSZT lda, ATL_CSZT ldb, ATL_CSZT ldc,
    const SCALAR alpha, const SCALAR beta);
static void Mjoin(PATL,GetRKInfo)
   (opinfo_t*a0, int a1, enum ATLAS_TRANS a2, enum ATLAS_TRANS a3,
    ATL_CSZT M, ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldb, size_t ldc,
    const SCALAR alpha, const SCALAR beta)
{
   ATL_assert(0);
}
int Mjoin(PATL,GetSyrkIdx)(unsigned int flg, ATL_CSZT N, ATL_CSZT K, double);
double Mjoin(PATL,sSyrkTimeEst)
   (int id, unsigned int flg, size_t N, size_t K, double symul);
int Mjoin(PATL,GetSyrkInfo)
   (amminfo_t*,int, enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,int);
void Mjoin(PATL,GetSyrkOP)
   (opinfo_t *out, int flag, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
    ATL_CSZT N, ATL_CSZT K, size_t lda, size_t ldc,
    const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetRankKInfo)
   (amminfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB, ATL_CSZT M,
    ATL_CSZT N, ATL_CSZT K, const SCALAR alpha, const SCALAR beta);
int Mjoin(PATL,GetAmmmInfo)
   (amminfo_t *out, enum ATLAS_TRANS TA, enum ATLAS_TRANS TB, ATL_CSZT M,
    ATL_CSZT N, ATL_CSZT K, const SCALAR alpha, const SCALAR beta);
#ifdef TREAL
int Mjoin(PATL,GetTrsmInfo)
   (amminfo_t *out, int ialp, enum ATLAS_TRANS TA, ATL_CSZT M, ATL_CSZT N,
    const SCALAR beta);
#endif
int Mjoin(PATL,tGetAmmmInfo)
   (amminfo_t *out, const unsigned int P, enum ATLAS_TRANS TA,
    enum ATLAS_TRANS TB, ATL_CSZT M, ATL_CSZT N, ATL_CSZT K,
    const SCALAR alpha, const SCALAR beta);
ablk2cmat_t Mjoin(PATL,tGetSyammInfo)
   (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K,
    const SCALAR alpha, const SCALAR beta);
ablk2cmat_t Mjoin(PATL,tGetSyammInfo_K)
   (amminfo_t *out, const int P, enum ATLAS_TRANS TA, ATL_CSZT N, ATL_CSZT K);

void Mjoin(PATL,oploopsM)
   (opinfo_t *rkinf, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    int MV, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC);
void Mjoin(PATL,oploopsN)
   (opinfo_t *rkinf, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    int MV, TYPE *a, TYPE *b, TYPE *rC, TYPE *iC);
void Mjoin(PATL,opblk)
   (opinfo_t *op, size_t i, size_t j, const TYPE *A, const TYPE *B, TYPE *C,
    TYPE *pA, TYPE *pAn, TYPE *pB, TYPE *pBn, TYPE *rC, TYPE *iC);

void Mjoin(PATL,ammmK)
   (amminfo_t*, const int mb, const int nmu, const int nb, const int nnu,
    ATL_CINT nfkblks, const int kb, const int kb0, const int KB0, const TYPE *A,
    const size_t lda, const size_t incAk, const TYPE*B, const size_t ldb,
    const size_t incBk, const ablk2cmat_t blkc2c, TYPE*, const size_t ldc,
    TYPE *a, ATL_CINT inca, TYPE *b, ATL_CINT incb, TYPE *rC, TYPE *iC,
    const SCALAR alpA, const SCALAR alpB, const SCALAR alpC, const SCALAR beta);

void Mjoin(PATL,gemm)
   (const enum ATLAS_TRANS,const enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT,
    const SCALAR, const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,
    TYPE*,ATL_CSZT);
void Mjoin(PATL,ammm)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
#ifdef TREAL
int Mjoin(PATL,rk4n4)(enum ATLAS_TRANS,ATL_CSZT,const SCALAR,const TYPE*,
                      ATL_CSZT,const TYPE*,ATL_CSZT,TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,ammm_rk2)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmREC)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmMNK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_aliased_rkK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_tN)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_IP)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_rkK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammm_1b)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmKMNK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmKNMK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmNKM)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT, const SCALAR,
    const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,TYPE*,ATL_CSZT);
int Mjoin(PATL,ammmNMK)
   (enum ATLAS_TRANS,enum ATLAS_TRANS,ATL_CSZT,ATL_CSZT,ATL_CSZT,
    const SCALAR,const TYPE*,ATL_CSZT,const TYPE*,ATL_CSZT,const SCALAR,
    TYPE*,ATL_CSZT);
#ifdef TCPLX
void Mjoin(PATL,heRL2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D);
void Mjoin(PATL,heLL2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D, ATL_iptr_t ldd);
void Mjoin(PATL,heRU2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D);
void Mjoin(PATL,heLU2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D, ATL_iptr_t ldd);
#endif
void Mjoin(PATL,syRL2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D);
void Mjoin(PATL,syLL2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D, ATL_iptr_t ldd);
void Mjoin(PATL,syRU2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D);
void Mjoin(PATL,syLU2ipBlk)
   (ipinfo_t *ip, ATL_CUINT bv, cm2am_t cpN, cm2am_t cpT, ATL_iptr_t IBLK,
    ATL_iptr_t KBLK, const TYPE *A, TYPE *W, TYPE *D, ATL_iptr_t ldd);
cm2am_t Mjoin(PATL,ipsyGetCopyA)(ipinfo_t*, ATL_CUINT flg, cm2am_t *CPT);
cm2am_t Mjoin(PATL,ipsyGetCopyB)(ipinfo_t*, ATL_CUINT flg, cm2am_t *CPT);

int Mjoin(PATL,ipsymmL)
   (ATL_CUINT, ATL_CSZT, ATL_CSZT, const SCALAR, const TYPE*,
    ATL_CSZT, const TYPE*, ATL_CSZT, const SCALAR, TYPE*, ATL_CSZT);
int Mjoin(PATL,ipsymmR)
   (ATL_CUINT, ATL_CSZT, ATL_CSZT, const SCALAR, const TYPE*,
    ATL_CSZT, const TYPE*, ATL_CSZT, const SCALAR, TYPE*, ATL_CSZT);
int Mjoin(PATL,ipsyr2k)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const TYPE*, ATL_CSZT,
    const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,ipher2k)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const TYPE*, ATL_CSZT,
    const TYPE, TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,opsyr2k)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const TYPE*, ATL_CSZT,
    const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,opher2k)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const TYPE*, ATL_CSZT,
    const TYPE, TYPE*,ATL_CSZT);
#endif
void Mjoin(PATL,ipsyrk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
void Mjoin(PATL,ipherk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const TYPE, const TYPE*, ATL_CSZT, const TYPE, TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,opsyrk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,opherk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const TYPE, const TYPE*, ATL_CSZT, const TYPE, TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,syrk_amm)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,herk_amm)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const TYPE, const TYPE*, ATL_CSZT, const TYPE, TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,umsyrk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,umherk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const TYPE, const TYPE*, ATL_CSZT, const TYPE, TYPE*,ATL_CSZT);
#endif
int Mjoin(PATL,sqsyrk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const SCALAR, const TYPE*, ATL_CSZT, const SCALAR, TYPE*,ATL_CSZT);
#ifdef TCPLX
int Mjoin(PATL,sqherk)
   (const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CSZT, ATL_CSZT,
    const TYPE, const TYPE*, ATL_CSZT, const TYPE, TYPE*,ATL_CSZT);
#endif
/*
 * Functions for testing Info return flag
 */
#define ATL_AMMFLG_KRUNTIME(flg_) ((flg_) & 1)
#define ATL_AMMFLG_KMAJOR(flg_) ((flg_) & 2)
/*
 * Helper functions we'd like to inline
 */
#ifdef ATL_GLOBELT
   #if defined(SCPLX) || defined(SREAL)
      #include "atlas_samm_kern.h"
   #else
      #if !defined(DREAL) && !defined(DCPLX)
         #error "TYPE must be defined for GLOBELT!"
      #endif
      #include "atlas_damm_kern.h"
   #endif
   #define ATL_GLOBIDX 1
#endif
#ifdef ATL_GLOBIDX
static TYPE INLINE *IdxAw_rkK(opinfo_t *op, TYPE *A, size_t i)
{  /* index A workspace for ith block for rank-K (only 1 kb) */
   size_t m;

   m = op->nfmblks;
   m = Mmin(m, i);
   i -= m;
   A += (m * op->szA + i * op->pszA)SHIFT;
   return(A);
}

static TYPE INLINE *IdxAw_ip(ipinfo_t *ip, TYPE *A, size_t i, size_t k)
{  /* index A workspace for inner-product gemm */
   const size_t nfmblks = ip->nfmblks, szA = (i < nfmblks) ? ip->szA:ip->pszA;
   size_t m;

   m = Mmin(nfmblks, i);
   i -= m;
   m = (m * ip->szA + i * szA) * (ip->nfkblks+1);  /* offset to kpan */
   m = (m+k*szA)SHIFT;
   return(A+m);
}

static TYPE INLINE *IdxBw_ip(ipinfo_t *ip, TYPE *B, size_t k, size_t j)
{  /* index B workspace for inner-product gemm */
   const size_t nfnblks = ip->nfnblks, szB = (j < nfnblks) ? ip->szB:ip->pszB;
   size_t n;

   n = Mmin(nfnblks, j);
   j -= n;
   n = (n * ip->szB + j * szB) * (ip->nfkblks+1);  /* offset to kpan */
   n = (n+k*szB)SHIFT;
   return(B+n);
}

#ifdef ATL_GLOBELT
/*
 * Get (i,j) elt out of block-major A storage.  Used for debugging,
 * so not optimized.  This code presently not debugged for KVEC.
 */
#ifdef TCPLX
   static void IdxAwElt_ip
      (ipinfo_t *ip, const TYPE *A, ATL_iptr_t i,ATL_iptr_t k, TYPE *elt)
#else
   static TYPE IdxAwElt_ip
      (ipinfo_t *ip, const TYPE *A, ATL_iptr_t i,ATL_iptr_t k)
#endif
{ /* get element (i,k) out of block storage */
   ATL_iptr_t iblk = i / ip->mb, kblk = k / ip->kb;
   ATL_UINT flg;
   const TYPE *a;
   #ifdef TCPLX
      ATL_iptr_t szA=ip->szA;
   #endif
   const unsigned int mu = ip->mu;

/*
 * Compute block number we are on (iblk), and i within block (i)
 */
   if (iblk >= ip->nfmblks)
   {
      ATL_iptr_t ii;
      i -= ip->mb * ip->nfmblks;
      ii = i / ip->pmb;
      iblk = ip->nfmblks + ii;
      i -= ii* ip->pmb;
      #ifdef TCPLX
         szA = ip->pszA;
      #endif
   }
   else
      i -= ip->mb * iblk;
   k -= kblk*ip->kb;           /* k now K index within block */
   a = IdxAw_ip(ip, (TYPE*)A, iblk, kblk); /* a pts to kb*mb block */
   flg = ATL_AMM_GetFLAG(ip->idxK);
/*
 * blk has is ku-major: mu kvec length blocks
 */
   if (ATL_AMM_KMAJOR(flg))
   {
      unsigned int ks = k/ip->vlen, is = i/mu;
      unsigned int bl=ip->vlen*mu, pl = ((ip->kb+ip->vlen-1)/ip->vlen)*bl;

      a += is*pl;        /* a pts to correct K-panel */
      i -= is*mu;
      a += ks*bl;        /* a pts to correct sub-block */
      k -= ks*ip->vlen;
      a += k*ip->vlen+i; /* a pts at correct element */
   }
/*
 * blk has standard mu-entry blks, repeated K times (vectorized along M)
 */
   else
   {
      ATL_CUINT is=i/mu, kb=(kblk < ip->nfkblks)?ip->kb:ip->KB0, pl=kb*mu;
      a += is*pl;         /* go to correct kpanel */
      i -= is*mu;         /* leave only i%mu remainder */
      a += k*mu;          /* go to correct mu block */
      a += i;             /* go to corect elt within mu blk */
   }
   #ifdef TCPLX
      elt[1] = *a;
      elt[0] = a[szA];
   #else
      return(*a);
   #endif
}

/*
 * Get (k,j) elt out of block-major B storage.  Used for debugging,
 * so not optimized.  This code presently not debugged for KVEC.
 */
#ifdef TCPLX
   static void IdxBwElt_ip
      (ipinfo_t *ip, const TYPE *B, ATL_iptr_t k, ATL_iptr_t j, TYPE *elt)
#else
   static TYPE IdxBwElt_ip
      (ipinfo_t *ip, const TYPE *B, ATL_iptr_t k, ATL_iptr_t j)
#endif
{ /* get element (k,j) out of access-major block storage */
   ATL_iptr_t jblk = j / ip->nb, kblk = k / ip->kb;
   ATL_UINT flg;
   const TYPE *b;
   #ifdef TCPLX
      ATL_iptr_t szB=ip->szB;
   #endif
   const unsigned int nu = ip->nu;

/*
 * Compute block number we are on (jblk), and j within block (j)
 */
   if (jblk >= ip->nfnblks)
   {
      ATL_iptr_t jj;
      j -= ip->nb * ip->nfnblks;
      jj = j / ip->pnb;
      jblk = ip->nfnblks + jj;
      j -= jj* ip->pnb;
      #ifdef TCPLX
         szB = ip->pszB;
      #endif
   }
   else
      j -= ip->nb * jblk;
   k -= kblk*ip->kb;           /* k now K index within block */
   b = IdxBw_ip(ip, (TYPE*)B, kblk, jblk); /* b pts to kb*mb block */
   flg = ATL_AMM_GetFLAG(ip->idxK);
/*
 * blk has is ku-major: nu kvec length blocks
 */
   if (ATL_AMM_KMAJOR(flg))
   {
      unsigned int ks = k/ip->vlen, js = j/nu;
      unsigned int bl=ip->vlen*nu, pl = ((ip->kb+ip->vlen-1)/ip->vlen)*bl;

      b += js*pl;        /* b pts to correct K-panel */
      j -= js*nu;
      b += ks*bl;        /* b pts to correct sub-block */
      k -= ks*ip->vlen;
      b += k*ip->vlen+j; /* b pts at correct element */
   }
/*
 * blk has standard nu-entry blks, repeated K times (vectorized along M)
 */
   else
   {
      ATL_CUINT js=j/nu, kb=(kblk < ip->nfkblks)?ip->kb:ip->KB0, pl=kb*nu;
      b += js*pl;         /* go to correct kpanel */
      j -= js*nu;         /* leave only i%nu remainder */
      b += k*nu;          /* go to correct nu block */
      b += j;             /* go to corect elt within nu blk */
   }
   #ifdef TCPLX
      elt[1] = *b;
      elt[0] = b[szB];
   #else
      return(*b);
   #endif
}
#endif

static const TYPE INLINE *IdxA_ip
   (ipinfo_t *ip, const TYPE *A, size_t i, size_t k)
{  /* index global A base-pointer to find (i,k) block for inner-product gemm */
   size_t m;

   m = ip->nfmblks;
   m = Mmin(m, i);
   i -= m;
   A += m * ip->incAm + i * ip->pincAm;
   A += k * ip->incAk;
   return(A);
}

static const INLINE TYPE *IdxB_ip
   (ipinfo_t *ip, const TYPE *B, size_t k, size_t j)
{  /* index global A base-pointer to find (i,k) block for inner-product gemm */
   size_t n;

   n = ip->nfnblks;
   n = Mmin(n, j);
   j -= n;
   B += ip->incBn * n + ip->pincBn * j;
   B += k * ip->incBk;
   return(B);
}
static INLINE TYPE *IdxC_ip(ipinfo_t *ip, TYPE *C, size_t i, size_t j)
{
   size_t n;

   n = ip->nfmblks;
   n = Mmin(n, i);
   i -= n;
   C += (ip->mb * n + ip->pmb * i)SHIFT;
   n = ip->nfnblks;
   n = Mmin(n, j);
   j -= n;
   C += (ip->nb * n + ip->pnb * j)*((ip->ldc) SHIFT);
   return(C);
}
#endif

#endif  /* end include file guard */
