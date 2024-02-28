/*   #define DEBUG 1 */
/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_trsmK.h"
#ifdef DEBUG
   #include "atlas_tst.h"
#endif
/*
 * Handles Left version of trsm for Lower,NoTrans or Upper,Trans
 *
 * On input, L has already be copied to a special format.  If nmblks=CEIL(M/b),
 * then L has nmblks row-panels, but since L is Lower triangular, the first
 * row panel has only one block (triangular), the next row panel contains
 * 2 blocks (first square, second triangular), and the last row panel
 * has nmblks-1 square blocks, followed by one triangular blocks.
 * Each row panel grabs blocks from L, and ends in the diagonal/triangular blk.
 * Square blocks are in the access-major storage required by the provided amm,
 * and have been scaled by -1 during the copy (so that they can be used to
 * subtract off contributions from the solved coefficients in prior triangular
 * blocks).
 * Triangular blocks are presently in block-major storage, with the diagonal
 * entries already inverted for NonUnit (and ignored for Unit), so that we
 * can call trsKL_rk4 directly.  Later, we need to make sure these guys
 * are multiples of 4, and then avoid the malloc/copy in trsmK_rk4.
 *
 * Workspace (w) must have room to store all blocks of B necessary to compute
 * one col-panel of R and one block of c.  Each block is size incw.
 * incw >= b*NB, where extra length can be added to preserve alignment.
 * So, total size is: (nmblks-1)*incw, incw >= b*NB.
 */
void Mjoin(PATL,ktrsmLLN_rk4)
   (ATL_CINT M, ATL_CINT N, const SCALAR alp, const TYPE *A, TYPE *B,
    ATL_CINT ldb, TYPE *W);
#define trsmK(Di_, m_, n_, al_, a_, b_, ldb_, w_) \
   Mjoin(PATL,ktrsmLLN_rk4)(m_, n_, al_, a_, b_, ldb_, w_)
int Mjoin(PATL,trsmKL_amm)
(
   amminfo_t *mmp,      /* chosen amm kernel info */
   enum ATLAS_DIAG Diag,
   ATL_CUINT mb0,       /* sizeof first triangular block */
   ATL_CUINT nmblks,    /* CEIL(M/mb); M = (nmblks-1)*b+mb0 */
   ATL_CUINT nnblks,    /* FLOOR(N/nb) */
   ATL_CUINT nr,        /* mod(N/nb) */
   const SCALAR alpha,  /* scale factor for rhs */
   const TYPE *L,       /* L, already copied to row-panel storage as above */
   size_t incL,         /* distance between blocks in L */
   TYPE *R,             /* ptr to col-major rhs (B on input, X on output */
   size_t ldr,          /* leading dim of RHS matrix R */
   ATL_CUINT nmu,       /* nmu = mb/mu, assume mod(mb/mu) == 0 */
   ATL_CUINT nnu,       /* nnu = NB/nu, assume mod(NB/nu) == 0 */
   TYPE *c,             /* bxNB with trailing readable space for amm */
   TYPE *W,             /* aligned workspace for amm's B blocks */
   size_t incW          /* stride between b*NB blocks */
)
{
   const int b=mmp->mb, NB=mmp->nb;
   const size_t incR = ldr*NB - (mb0+(nmblks-1)*b);
   ATL_CUINT ku=mmp->ku, ldl0=(mb0+3)&(~3);
   ATL_UINT kb0;
   ATL_UINT r;
   TYPE *pR = R;
   cm2am_t b2blk=mmp->b2blk;       /* col2blk for B, wt alpha=1 */
   ablk2cmat_t blk2c=mmp->Cblk2cm; /* blk2col for C, wt beta of alpha above */
   ammkern_t amm_b0=mmp->amm_b0, amm_b1=mmp->amm_b1;

   if (mb0 != b)
   {
      if (!(mmp->flag & 1))
         amm_b0 = mmp->amm_k1_b0;
      else if (mb0 < mmp->kbmin || (mb0/ku)*ku != mb0)
         amm_b0 = mmp->amm_k1_b0;
   }
/*
 * Loop over the full column-panels of RHS
 */
   for (r=0; r < nnblks; r++, pR += incR)
   {
      ATL_INT i;
      const TYPE *pLb = L+incL;  /* ptr to L's block, loops along rowpan */
/*
 *    First row-panel has no gemm updates
 */
      #ifdef DEBUG
         i = (mb0+3)&(~3);
         Mjoin(PATL,geprint)("L0", i, i, L, i);
      #endif
      trsmK(Diag, mb0, NB, alpha, L, pR, ldr, c);
      #ifdef DEBUG
         Mjoin(PATL,geprint)("B~", mb0, NB, pR, ldr);
      #endif
      if (nmblks > 1)
         b2blk(mb0, NB, ATL_rone, pR, ldr, W); /* put B0 in amm's B storage */
      pR += mb0;
      for (i=1; i < nmblks; i++)     /* loop over rowpans of L and the */
      {                              /* blks in the colpan of B they update */
         TYPE *w = W;
         ammkern_t amm=amm_b0;
         ATL_INT j, kb=mb0;

         for (j=0; j != i; j++)      /* apply gemm updates to C */
         {
            const TYPE *LN = pLb + incL;
            TYPE *wN = w + incW;
            amm(nmu, nnu, kb, pLb, w, c, LN, wN, c);
            pLb = LN;
            w = wN;
            amm = amm_b1;
            kb = b;
         }
/*
 *       Subtract all solved co-efficients from RHSs using accumulated GEMM
 *       updates.  We apply alpha here, so the triangular solve will use alpha=1
 */
         blk2c(b, NB, ATL_rone, c, alpha, pR, ldr);
/*
 *       Now complete solve for X on this block of RHS with TRSM
 */
         #ifdef DEBUG
            Mjoin(PATL,geprint)("Lx", b, b, pLb, b);
         #endif
         trsmK(Diag, b, NB, ATL_rone, pLb, pR, ldr, c);
         pLb += incL;
/*
 *       If this isn't the last row-panel, I'll need the final answer copied
 *       to amm's B storage for subtracting these equations' contributions out
 */
         if (i != nmblks-1)
            b2blk(b, NB, ATL_rone, pR, ldr, w); /* put Bi in amm's B storage */
         pR += b;
      }
   }
/*
 * Is there a partial NB block at end?  This block must handle padding for both
 * B and C.
 */
   if (nr)
   {
      const int nu=mmp->nu, nnu=(nr+nu-1)/nu, NB=nnu*nu;
      ATL_INT i;
      const TYPE *pLb = L+incL;  /* ptr to L's block, loops along rowpan */
/*
 *    First row-panel has no gemm updates
 */
      #ifdef DEBUG
         i = (mb0+3)&(~3);
         Mjoin(PATL,geprint)("L0", i, i, L, i);
      #endif
      trsmK(Diag, mb0, nr, alpha, L, pR, ldr, c);
      #ifdef DEBUG
         Mjoin(PATL,geprint)("B~", mb0, NB, pR, ldr);
      #endif
      if (nmblks > 1)
         b2blk(mb0, nr, ATL_rone, pR, ldr, W); /* put B0 in amm's B storage */
      pR += mb0;
      for (i=1; i < nmblks; i++)     /* loop over rowpans of L and the */
      {                              /* blks in the colpan of B they update */
         TYPE *w = W;
         ammkern_t amm=amm_b0;
         ATL_INT j, kb=mb0;

         for (j=0; j != i; j++)      /* apply gemm updates to C */
         {
            const TYPE *LN = pLb + incL;
            TYPE *wN = w + incW;
            amm(nmu, nnu, kb, pLb, w, c, LN, wN, c);
            pLb = LN;
            w = wN;
            amm = amm_b1;
            kb = b;
         }
/*
 *       Subtract all solved co-efficients from RHSs using accumulated GEMM
 *       updates.  We apply alpha here, so the triangular solve will use alpha=1
 */
         blk2c(b, nr, ATL_rone, c, alpha, pR, ldr);
/*
 *       Now complete solve for X on this block of RHS with TRSM
 */
         #ifdef DEBUG
            Mjoin(PATL,geprint)("Lx", b, b, pLb, b);
         #endif
         trsmK(Diag, b, nr, ATL_rone, pLb, pR, ldr, c);
         pLb += incL;
/*
 *       If this isn't the last row-panel, I'll need the final answer copied
 *       to amm's B storage for subtracting these equations' contributions out
 */
         if (i != nmblks-1)
            b2blk(b, nr, ATL_rone, pR, ldr, w); /* put Bi in amm's B storage */
         pR += b;
      }
   }
   return(0);
}

#define L2blk Mjoin(PATL,trsmK_L2blk)
int Mjoin(PATL,trsmL_amm)
(
   enum ATLAS_SIDE Side,
   enum ATLAS_UPLO Uplo,
   enum ATLAS_TRANS TA,
   enum ATLAS_DIAG Diag,
   ATL_CINT  M,        /* size of triangular matrix A */
   ATL_CINT  N,        /* number of RHS in B */
   const SCALAR alpha, /* scale factor for B */
   const TYPE *T,      /* MxM triangular matrix A */
   ATL_CINT  ldt,
   TYPE *R,            /* on input, B, on output X, of A x = b */
   ATL_CINT  ldr       /* leading dim of R */
)
{
   amminfo_t mminf;
   ATL_UINT nmblks, nsqblks, nnblks, nb, nmu, nnu, mu, nu, nr;
   size_t wrksz, szL, szB, szC, incL, incW;
   void *vp=NULL;
   TYPE *pC, *W, *L;
   ATL_INT i;
   int b, mb0;
   #ifdef DEBUG
      int *eBV=NULL;
      TYPE *myB;

      myB = malloc(M*N*sizeof(TYPE));
      ATL_assert(myB);
      Mjoin(PATL,geprint)("Lg", M, M, T, ldt);
      Mjoin(PATL,gecopy)(M, N, R, ldr, myB, M);
      if (Diag == AtlasUnit)
         Mjoin(PATL,reftrsmLLNU)(M, N, alpha, T, ldt, myB, M);
      else
         Mjoin(PATL,reftrsmLLNN)(M, N, alpha, T, ldt, myB, M);
      Mjoin(PATL,geprint)("Xg", M, N, myB, M);
   #endif
   if (Side == AtlasRight)
      return(1);
   if (Uplo == AtlasUpper && TA == AtlasNoTrans)
      return(2);
   if (Uplo == AtlasLower && TA == AtlasTrans)
      return(3);
   #if 0
   if (M < 24)   /* tiny triangles should just use trsmKL_rk4 directly */
      return(4); /* since our B can't be much bigger than 4 */
   if (N < 4)
      return(5);
   #endif
   if (alpha == 0.0)
      i = 0;
   else if (alpha == 1.0)
      i = 1;
   else
      i = (alpha == -1.0) ? -1 : 2;
   b = Mjoin(PATL,GetTrsmInfo)(&mminf, i, TA, M, N, alpha);
   if (b < 4)
      return(Mjoin(PATL,trsmKL_rk4)(Side, Uplo, TA, Diag, M, N, alpha,
                                    T, ldt, R, ldr));

   nb = mminf.nb;
   mu = mminf.mu;
   nu = mminf.nu;
   nmblks = ((M+b-1)/b);
   nnblks = N / nb;
   nr = N - nnblks*nb;
   nsqblks = ((nmblks-1)*nmblks)>>1; /* # of square blks below diag */
   incL = (b+3)&(~3);                /* triang blks must be multiples of 4 */
   incL *= incL;
   incW = b * nb;
   szL = incL * (nsqblks + nmblks);  /* don't store upper blks of matrix */
   szB = (nmblks-1) * incW;          /* last blk not copied to B form */
   szC = incW;                       /* C just 1 block in size */
   wrksz = 3*ATL_Cachelen + ATL_MulBySize(szL+szB+szC);
   if (wrksz < ATL_MaxMalloc)
      vp = malloc(wrksz);
   if (!vp)
      return(4);
   pC = ATL_AlignPtr(vp);
   W = pC + szC;            /* W is workspace for amm's B */
   W = ATL_AlignPtr(W);
   L = W + szB;             /* L is workspace for triangular matrix */
   L = ATL_AlignPtr(L);     /* put it last so last diag blk is amm buffer */
   i = nmblks*b;
   if (i == M)
      mb0 = b;
   else
      mb0 = M - i + b;
   nmu = b / mu;
   nnu = nb / nu;
   L2blk(Diag, mb0, b, nmblks, T, ldt, L, incL, mminf.a2blk);
   ATL_assert(!Mjoin(PATL,trsmKL_amm)(&mminf, Diag, mb0, nmblks, nnblks,
                                      nr, alpha, L, incL, R, ldr, nmu, nnu,
                                      pC, W, incW));
   #ifdef DEBUG
      ATL_assert(nmu*mu == b);
      ATL_assert(b*(nmblks-1)+mb0 == M);
      i = (M+3)&(~3);
      eBV = Mjoin(PATL,cmpmatBV)(0, (1.0*M)*M*1e-15, M, N, myB, M, R, ldr);
      Mjoin(PATL,geprint)("Xc", M, N, R, ldr);
      Mjoin(PATL,gediff)(M, N, R, ldr, myB, M, myB, M);
      ATL_print2dBV(M, N, eBV);
      free(eBV);
      free(myB);
   #endif
   free(vp);
   return(0);
}
