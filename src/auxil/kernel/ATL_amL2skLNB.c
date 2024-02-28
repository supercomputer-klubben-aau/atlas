/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2018 R. Clint Whaley
 */
#include "atlas_aux.h"
/*
 * Takes access-major 'Lower' C block W, and translates to column-major
 * lower symm C storage C.  This function operates on 1 block at most.
 * shVL is power of vector length.
 */
#ifdef BETA1
   #define BNM _b1
#elif defined(BETA0)
   #define BNM _b0
#elif defined(BETAN) || defined(BETAN1)
   #define BNM _bN
#else
   #define BNM _bX
#endif
void Mjoin(Mjoin(PATL,amL2skLNB),BNM)
   (ATL_iptr_t N, ATL_CUINT MU, ATL_CUINT NU, ATL_CUINT shVL,
    const TYPE *W, const SCALAR beta, TYPE *C, ATL_iptr_t ldc)
{
#if 0
   ATL_CUINT incW = ((MU*NU+(1<<shVL)-1)>>shVL)<<shVL;
   ATL_UINT i, in;

   for (i=0; i < N; i = in)
   {
      ATL_UINT j, jn, mu = N-i;
      mu = (mu >= MU) ? MU : mu;
      in = i + mu;
      for (j=0; j < in; j = jn, W += incW)
      {
         ATL_UINT nu = N-j;
         TYPE *c = C+i+j*ldc;
         const TYPE *w = W;
         nu = (nu >= NU) ? NU : nu;
         jn = j + nu;
         if (jn <= i) /* block does not cross diagonal */
         {
            ATL_UINT ii, jj;
            for (jj=0; jj < nu; jj++, c += ldc, w += MU)
               for (ii=0; ii < mu; ii++)
                  c[ii] = w[ii];
         }
         else  /* diagonal crosses block */
         {
            ATL_UINT ii, jj;
            for (jj=0; jj < nu; jj++, c += ldc, w += MU)
            {
               ATL_CUINT J=j+jj;
               for (ii=0; ii < mu; ii++)
               {
                  ATL_CUINT I=i+ii;
                  if (I >= J)
                     c[ii] = w[ii];
               }
            }
         }
      }
   }
#else
   ATL_CUINT blksz = ((MU*NU+(1<<shVL)-1)>>shVL)<<shVL;
   ATL_UINT j, jn, cinc;
   const TYPE *rW=W, *rw;
   const ATL_iptr_t incC = ldc*NU;
   for (cinc=j=0; j < N; j = jn, C += incC, cinc += blksz)
   {
      ATL_UINT i, in, nu = N-j, nd;
      nu = (nu >= NU) ? NU : nu;
      jn = j + nu;
      i = (j/MU)*MU;
/*
 *    Diagonal presently at (j,j), compute how many blocks of i diag crosses,
 *    crossing should happen at (jn,jn)
 */
      rw = rW;
      nd = i + ((jn-i+MU-1)/MU)*MU;
      for(nd=Mmin(nd,N); i < nd; i = in)
      {
         ATL_UINT mu=N-i, ii, jj, nblks;
         const TYPE *w = rw + cinc;
         TYPE *c = C + i;
         mu = (mu >= MU) ? MU : mu;
         in = i + mu;
         for (jj=0; jj < nu; jj++, c += ldc, w += MU)
         {
            ATL_CUINT J=j+jj;
            for (ii=0; ii < mu; ii++)
            {
               ATL_CUINT I=i+ii;
               if (I >= J)
                  c[ii] = w[ii];
            }
         }
         nblks = (in+NU-2) / NU;  /* # of row blocks in this panel */
         rw += nblks*blksz; /* finished all blocks in rowpanel */
/*
 *       Diagonal blocks may finish rows, if so increment rW
 */
         rW = (in <= jn) ? rw : rW;
      }
/*
 *    Diagonal above any remaining blocks
 */
      for (; i < N; i = in)
      {
         ATL_UINT mu=N-i, ii, jj, nblks;
         const TYPE *w = rw + cinc;
         TYPE *c = C + i;
         mu = (mu >= MU) ? MU : mu;
         in = i + mu;
         for (jj=0; jj < nu; jj++, c += ldc, w += MU)
            for (ii=0; ii < mu; ii++)
               c[ii] = w[ii];
         nblks = (in+NU-2) / NU;  /* # of row blocks in this panel */
         rw += nblks*blksz;
      }
   }
#endif
}
