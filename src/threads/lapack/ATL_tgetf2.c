#include "atlas_cbc.h"
#include "atlas_tlapack.h"
#include "atlas_level2.h"
#ifdef TCPLX
   #define GER Mjoin(PATL,geru)
#else
   #define GER Mjoin(PATL,ger)
#endif
int Mjoin(PATL,tgetf2C)
(
   ATL_CUINT P,   /* number of threads working on this problem */
   ATL_CUINT rank,/* rank of this thread */
   void *vchk,    /* boolean sync array */
   ATL_CUINT M,   /* number of rows in A;  M >= P*N for parallelism! */
   ATL_CUINT N,   /* number of cols in A */
   TYPE *A,       /* IN: matrix to factor; OUT: LU */
   size_t lda1, /* element stride between cols of A */
   int *ipiv      /* pivot array */
)
{
   #ifdef TCPLX
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
      size_t lda2 = lda1+lda1;
   #else
      #define none ATL_rnone
      #define lda2 lda1
   #endif
   int info=0;
/*
 * Rank 0 owns the diagonal, so these cores will just perform updates as
 * driven by rank 0
 */
   ATL_INT j, mp, mr;
   ATL_CUINT Mp = M / P;     /* # of rows all cores get */
   ATL_CUINT Mr = M - Mp*P;  /* extra rows left over */
/*
 * Quick return if no work to do
 */
   if (M < 1 || N < 1)
      return(0);
/*
 * For parallism, the least number of local rows must be larger than N; i.e.
 * i.e., the whole diagonal block must be owned by the rank 0.
 */
   if (Mp < N)
      return(Mjoin(PATL,getf2)(M, N, A, lda1, ipiv));
/*
 * All ranks except 0 have only trailing matrix to worry about, so their
 * local row count never changes, they don't have a diagonal, etc.
 */
   if (rank)
   {
      ATL_CUINT Mprev = Mp*rank + Mmin(rank-1, Mr);
      TYPE *Ac=A+(Mprev SHIFT);
      ATL_CUINT ML = (rank <= Mr) ? Mp+1 : Mp;

      for (j=0; j < N; j++, Ac += lda2)
      {
         #ifdef TCPLX
            TYPE pv[2];
         #else
            TYPE pv;
         #endif
         int pivL, pivG;
         pivL = cblas_iamax(ML, Ac, 1);
         pivG = pivL + Mprev;
         #ifdef TCPLX
            pv[0] = Ac[pivL+pivL];
            pv[1] = Ac[pivL+pivL+1];
            Mjoin(PATL,cbc_comb_iamax_nopost0)(P, rank, &pivG, pv, vchk);
         #else
            pv = Ac[pivL];
            Mjoin(PATL,cbc_comb_iamax_nopost0)(P, rank, &pivG, &pv, vchk);
         #endif
         #ifdef TCPLX
            if (pv[0] != ATL_rzero || pv[1] != ATL_rzero)
            {
               Mjoin(PATL,cplxinvert)(1, pv, 1, pv, 1);
               cblas_scal(ML, pv, Ac, 1);
            }
         #else
            if (pv != ATL_rzero)
               cblas_scal(ML, ATL_rone/pv, Ac, 1);
         #endif
            else if (!info)
               info = j+1;
         GER(ML, N-j-1, none, Ac,1, A+(((j+1)*lda1+j)SHIFT),lda1, Ac+lda2,lda1);
      }
   }
/*
 * Rank 0 owsn the diagonal block and so it will drive the algorithm,
 * and do pivoting and solving.
 */
   else
   {
      size_t ldap1 = (lda1+1)SHIFT;
      TYPE *Ac=A, *Ad=A;
      int m = Mp;

      for (j=0; j < N; j++, Ac += lda2, Ad += ldap1)
      {
         #ifdef TCPLX
            TYPE pv[2];
         #else
            TYPE pv;
         #endif
         int pivL, pivG;
         int bval;
         pivL = cblas_iamax(m--, Ad, 1);
         pivG = j + pivL;
         #ifdef TCPLX
            pv[0] = Ad[pivL+pivL];
            pv[1] = Ad[pivL+pivL+1];
            bval = Mjoin(PATL,cbc_comb_iamax_nopost0)(P, 0, &pivG, pv, vchk);
         #else
            pv = Ad[pivL];
            bval = Mjoin(PATL,cbc_comb_iamax_nopost0)(P, 0, &pivG, &pv, vchk);
         #endif
/*         printf("pivG=%d\n", pivG); */
         ipiv[j] = pivG;
         if (pivG != j)
            cblas_swap(N, A+(j SHIFT), lda1, A+(pivG SHIFT), lda1);
         ATL_cbc_post(0, vchk);  /* allow other cores to leave iamax! */
         #ifdef TCPLX
            if (pv[0] != ATL_rzero || pv[1] != ATL_rzero)
            {
               Mjoin(PATL,cplxinvert)(1, pv, 1, pv, 1);
               cblas_scal(m, pv, Ad+2, 1);
            }
         #else
            if (pv != ATL_rzero)
               cblas_scal(m, ATL_rone/pv, Ad+1, 1);
         #endif
            else if (!info)
               info = j+1;
         GER(m, N-j-1, none, Ad+(1 SHIFT), 1, Ad+lda2, lda1, Ad+ldap1, lda1);
      }
   }
   return(info);
}
#undef GER
#ifdef none
   #undef none
#endif
#ifdef lda2
   #undef lda2
#endif
#if 0  /* untested code! */
int Mjoin(PATL,tgetf2CR)
(
   ATL_CUINT P,   /* number of threads working on this problem */
   ATL_CUINT rank,/* rank of this thread */
   void *vchk,    /* boolean sync array */
   ATL_CUINT M,   /* number of rows in A;  M >= P*N for parallelism! */
   ATL_CUINT N,   /* number of cols in A */
   TYPE *A,       /* IN: matrix to factor; OUT: LU */
   size_t lda,    /* element stride between cols of A */
   int *ipiv      /* pivot array */
)
{
   #ifdef TCPLX
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
      const TYPE one[2]  = {ATL_rone, ATL_rzero};
      size_t lda2 = lda+lda;
   #else
      #define one ATL_rone
      #define none ATL_rnone
      #define lda2 lda
   #endif
   int ierr=0;
   ATL_CUINT Nleft = (N>>2)<<1;
   if (Nleft)
   {
       ATL_CUINT Nright=N-Nleft, NRp = Nright/P, NRpr = Nright-NRp*P;
       ATL_CUINT nrp = (rank < NRpr) ? NRp+1: NRp;
       ATL_CUINT nrprev = (NRp*rank + Mmin(rank, NRpr)) SHIFT;
       ATL_CUINT MN = Mmin(M,N);
       ATL_INT MM, Mp;
       TYPE *Ac, *An;
       int i;
/*
 *     Factor left portion, the apply pivots to right, then sync so we
 *     we know it is save to start writing TRSM on pivoted matrix
 */
       ierr = Mjoin(PATL,tgetf2CR)(P, rank, vchk, M, N, A, lda, ipiv);
       Ac = A + Nleft*lda2;
       if (nrp)
          ATL_laswp(nrp, Ac+nrprev, lda, 0, Nleft, ipiv, 1);
       An = Ac + (Nleft SHIFT);
       ATL_barrier();
       ATL_membar;   /* core must 'see' all changes made above! */
/*
 *     If local NRHS < 4, reduce parellism to avoid L1BLAS perf with contention
 */
       if (NRp < 4)
       {
          ATL_CUINT p4 = (Nright > 4) ? (Nright>>2) : 1;
          if (rank < p4)
          {
             if (p4 > 1)
             {
                ATL_CUINT np = Nright/p4, nr = Nright-np*p4;
                ATL_CUINT n = (rank < nr) ? np+1: np;
                ATL_CUINT prev = (np*rank + Mmin(rank, nr)) SHIFT;
                cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                           CblasUnit, Nleft, n, one, A, lda, Ac+prev, lda);
             }
             else
                cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                           CblasUnit, Nleft, Nright, one, A, lda, Ac, lda);
          }
       }
       else
          cblas_trsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
                     CblasUnit, Nleft, nrp, one, A, lda, Ac+nrprev, lda);
/*
 *     Await completion of solve, and update with GEMM by splitting M.
 *     This algorithm really designed for panel factorization, where M
 *     dominates N.  Would not be good for problems where M is very small,
 *     and N very large!  Can add case to that affect if it becomes important.
 */
       MM = M - Nleft;
       Mp = MM / P;
       ATL_barrier();
       ATL_membar;   /* core must 'see' all changes made above! */
       if (Mp >= 4)
       {
          ATL_CUINT mr = M-Mp*P, m = (rank < mr) ? Mp+1: Mp;
          ATL_CUINT prev = (Mp*rank + Mmin(rank, mr)) SHIFT;
          cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, Nright,
                     Nleft, none, A+prev+(Nleft SHIFT), lda,
                     Ac+prev, lda, one, An, lda);
       }
       else if (rank == P-1)
          cblas_gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, MM, Nright,
                     Nleft, none, A+(Nleft SHIFT), lda, Ac, lda, one, An, lda);
/*
 *     Await completion of GEMM and then factor the right
 */
       ATL_barrier();
       ATL_membar;   /* core must 'see' all changes made above! */
       i = Mjoin(PATL,tgetf2CR)(P, rank, vchk, MM, Nright, An, lda, ipiv+Nleft);
/*
 *     Must use a special laswp rout that adds 1st arg to all piv entries.
 *     this avoids having a sync for update of ipiv, which is handled last,
 *     when we are forced to sync anyway
 * HERE HERE: still need to write this rout, and then parallelize this call!
 */

       ATL_laswp_off(Nleft, A, lda, Nleft, MN, ipiv, 1);
       if (i) if (!ierr) ierr = i + Nleft;
       if (rank == P-1)
          for (i=Nleft; i != MN; i++)
             ipiv[i] += Nleft;
   }
   else ierr = Mjoin(PATL,tgetf2C)(P, rank, vchk, M, N, A, lda, ipiv);
/*
 * Not safe to use data until everybody sees changes, so don't ret until done!
 */
   ATL_barrier();
   ATL_membar;   /* core must 'see' all changes before returning */
   return(ierr);
}
#ifdef one
   #undef one
#endif
#ifdef none
   #undef none
#endif
#ifdef lda2
   #undef lda2
#endif
#endif
#include "atlas_pca.h"

void Mjoin(PATL,DoWorkGETF2_nowrk)(ATL_LAUNCHSTRUCT_t *lp, void *vp)
{
   ATL_thread_t *tp=vp;
   ATL_TGETF2_M_t *lup=((ATL_TGETF2_M_t*)lp->opstruct)+tp->rank;
   int *ipiv = lup->ipiv;
   ATL_CINT M=lup->M, N=lup->N, lda=lup->lda, MN = Mmin(M,N);
   const int p = lup->p, rank = lup->rank;
   ATL_CINT mp = M/p, mr = M - mp*p;
   ATL_INT m, locpiv, globpiv, k, j, i;
   #ifdef TCPLX
      ATL_CINT lda2 = lda+lda;
   #else
      #define lda2 lda
      #define none ATL_rnone
   #endif
   TYPE *A, *Ac, *a, *v;
   TYPE pivval, apv, apv2;
   #ifdef TCPLX
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #endif
   volatile ATL_INT *maxindx=lup->maxindx, *stage=lup->stage;
   void (*my_ger)(const int M, const int N, const SCALAR alpha,
                  const TYPE *X, const int incX,
                  const TYPE *Y, const int incY, TYPE *A, const int lda);

   #if 1
      if (M >= N*p && p > 1)
      {
         int ierr;
         ierr = Mjoin(PATL,tgetf2C)(p, rank, NULL, M, N, lup->A, lda, ipiv);
         if (ierr)
            lup->info = ierr;
         return;
      }
   #endif
   #ifdef TCPLX
      my_ger = Mjoin(PATL,geru);
   #else
      my_ger = Mjoin(PATL,ger);
   #endif
   m = (rank) ? mp : mp+mr;
   Ac = A = lup->A;
   a = (rank) ? A + ((m*rank + mr)SHIFT) : A;
   for (j=0; j < MN; j++, Ac += lda2, a += lda2)
   {
      locpiv = cblas_iamax(m, a, 1);
/*
 *    Combine local pivot into global
 */
      if (!rank)
      {
         globpiv = j+locpiv;
         #ifdef TCPLX
            apv = Mabs(Ac[globpiv+globpiv]) + Mabs(Ac[globpiv+globpiv+1]);
         #else
            apv = Mabs(Ac[globpiv]);
         #endif
         #if 1
            Mjoin(PATL,cbc_comb_iamax_nopost0)(p, rank, &globpiv, &apv, NULL);
         #else
         for (i=1; i < p; i++)
         {
            while(stage[i] < j);
            k = maxindx[i];
            apv2 = Mabs(Ac[k SHIFT]);
            #ifdef TCPLX
               apv2 += Mabs(Ac[k+k+1]);
            #endif
            if (apv < apv2)
            {
               apv = apv2;
               globpiv = k;
            }
            maxindx[i] = -1;
         }
         #endif
         ipiv[j] = globpiv;
/*         printf("globpiv=%d\n", globpiv); */
         if (globpiv != j)
            cblas_swap(N, A+(j SHIFT), lda, A+(globpiv SHIFT), lda);
         ATL_cbc_post(rank, NULL);
         stage[0] = j;
         m--;                                           /* just finished */
         #ifdef TCPLX
            a += 2;                                     /* one row */
         #else
            a++;                                        /* one row */
         #endif
      }
      else /* all threads except 0 write their results, and await 0 */
      {
         #ifdef TCPLX
            apv = Mabs(a[locpiv+locpiv]) + Mabs(a[locpiv+locpiv+1]);
         #else
            apv = Mabs(a[locpiv]);
         #endif
         #if 1
            globpiv = locpiv+rank*mp+mr;
            Mjoin(PATL,cbc_comb_iamax_nopost0)(p, rank, &globpiv, &apv, NULL);
         #else
         maxindx[rank] = locpiv+rank*mp+mr;
         stage[rank] = j;
         while (stage[0] < j);
         #endif
      }
      #ifdef TCPLX
         if (Ac[j+j] != ATL_rzero || Ac[j+j+1] != ATL_rzero)
         {
            TYPE inv[2];
            Mjoin(PATL,cplxinvert)(1, Ac+j+j, 1, inv, 1);
            cblas_scal(m, inv, a, 1);
         }
      #else
         pivval = Ac[j];
         if (pivval != ATL_rzero)
            cblas_scal(m, ATL_rone/pivval, a, 1);
      #endif
      else /* pivot is zero, we have a singular matrix! */
         lup->info = j;   /* all threads have same info */

      #ifdef TCPLX
         my_ger(m, N-j-1, none, a, 1, Ac+((j+lda)<<1), lda, a+lda2, lda);
         my_ger = Mjoin(PATL,geru_L2);
      #else
         my_ger(m, N-j-1, ATL_rnone, a, 1, Ac+j+lda, lda, a+lda, lda);
         my_ger = Mjoin(PATL,ger_L2);
      #endif
   }
}

void Mjoin(PATL,DoWorkGETF2)(ATL_LAUNCHSTRUCT_t *lp, void *vp0)
{
   ATL_thread_t *tp=vp0;
   ATL_TGETF2_M_t *lup=((ATL_TGETF2_M_t*)lp->opstruct)+tp->rank;
   int *ipiv = lup->ipiv;
   ATL_CINT M=lup->M, N=lup->N, lda=lup->lda, MN = Mmin(M,N);
   const int p = lup->p, rank = lup->rank;
   int pivrank;
   ATL_CINT mp = M/p, mr = M - mp*p;
   ATL_INT m, locpiv, globpiv, k, j, i, ldw, ldw0, ldw1;
   void *vp;
   TYPE *a, *W, *Wc, *w, **WRKS=lup->works, *v;
   TYPE pivval, apv, apv2, pv2;
   volatile ATL_INT *maxindx=lup->maxindx, *stage=lup->stage;
   #ifdef TCPLX
      const TYPE none[2] = {ATL_rnone, ATL_rzero};
   #endif

   m = (rank) ? mp : mp+mr;
   a = (rank) ? (((TYPE*)lup->A)+((mp*rank + mr)SHIFT)) : lup->A;
/*
 * Make ldw's a multiple of 16 bytes that is not a power of 2; 0's ldw
 * is larger by mr than all other ldws (ldw1)
 */
#if defined(DREAL) || defined(SCPLX)
   ldw0 = ((mp+mr+1)>>1)<<1;
   ldw1 = ((mp+1)>>1)<<1;
   if (!(ldw0 & (ldw0-1)))
      ldw0 += 2;
   if (!(ldw1 & (ldw1-1)))
      ldw1 += 2;
#elif defined(SREAL)
   ldw0 = ((mp+mr+3)>>2)<<2;
   ldw1 = ((mp+3)>>2)<<2;
   if (!(ldw0 & (ldw0-1)))
      ldw0 += 4;
   if (!(ldw1 & (ldw1-1)))
      ldw1 += 4;
#else
   ldw0 = mp+mr;
   ldw1 = mp;
   if (!(ldw0 & (ldw0-1)))
      ldw0++;
   if (!(ldw1 & (ldw1-1)))
      ldw1++;
#endif
   ldw = (rank) ? ldw1 : ldw0;
   vp = malloc(ATL_MulBySize(ldw)*N+ATL_Cachelen);
/*
 * If anyone fails to allocate the space, free any allocated spaces and
 * call the no-copy version
 */
   j = (vp != NULL);
   if (!rank)
   {
      for (i=1; i < p; i++)
      {
         while (stage[i] != -2);
         j &= maxindx[i];
         maxindx[i] = -1;
      }
      *maxindx = j;
      stage[0] = -2;
   }
   else
   {
      maxindx[rank] = j;
      stage[rank] = -2;
      while (stage[0] != -2);
   }
   if (*maxindx == 0)
   {
      if (vp)
         free(vp);
      Mjoin(PATL,DoWorkGETF2_nowrk)(lp, vp0);
      return;
   }
   ATL_assert(vp);
   WRKS[rank] = w = W = ATL_AlignPtr(vp);
   Mjoin(PATL,gecopy)(m, N, a, lda, W, ldw);
   for (j=0; j < MN; j++, w += (ldw SHIFT))
   {
      locpiv = cblas_iamax(m, w, 1);
      apv = Mabs(w[locpiv SHIFT]);
      #ifdef TCPLX
         apv += Mabs(w[locpiv+locpiv+1]);
      #endif
/*
 *    Combine local pivot into global
 */
      if (!rank)
      {
         globpiv = j+locpiv;
         Mjoin(PATL,cbc_comb_iamax_nopost0)(p, rank, &globpiv, &apv, NULL);
         pivrank = (globpiv > mp+mr) ? (globpiv - mr) / mp : 0;
         ipiv[j] = globpiv;
         if (pivrank)
         {
            locpiv = globpiv-mr-pivrank*mp;
            cblas_swap(N, W+(j SHIFT), ldw,
                       WRKS[pivrank]+(locpiv SHIFT), ldw1);
         }
         else
         {
            if (globpiv != j)
               cblas_swap(N, W+(j SHIFT), ldw, W+(globpiv SHIFT), ldw);
         }
         ATL_cbc_post(rank, NULL);
         stage[0] = j;
         m--;                                           /* just finished */
         #ifdef TCPLX
            w += 2;                                     /* one row */
         #else
            w++;                                        /* one row */
         #endif
      }
      else /* all threads except 0 write their results, and await 0 */
      {
         globpiv = mr+rank*mp+locpiv;
         Mjoin(PATL,cbc_comb_iamax_nopost0)(p, rank, &globpiv, &apv, NULL);
      }
      #ifdef TCPLX
         v = &WRKS[0][(j*ldw0+j)SHIFT];
         if (*v != ATL_rzero || v[1] != ATL_rzero)
         {
            TYPE inv[2];
            Mjoin(PATL,cplxinvert)(1, v, 1, inv, 1);
            cblas_scal(m, inv, w, 1);
         }
      #else
         pivval = WRKS[0][j*ldw0+j];
         if (pivval != ATL_rzero)
            cblas_scal(m, ATL_rone/pivval, w, 1);
      #endif
      else /* pivot is zero, we have a singular matrix! */
         lup->info = j;   /* all threads have same info */

      #ifdef TCPLX
         Mjoin(PATL,geru_L2)(m, N-j-1, none, w, 1,
                             WRKS[0]+((j*(ldw0+1)+ldw0)SHIFT), ldw0,
                             w+ldw+ldw, ldw);
      #else
         Mjoin(PATL,ger_L2)(m, N-j-1, ATL_rnone, w, 1, WRKS[0]+j*(ldw0+1)+ldw0,
                            ldw0, w+ldw, ldw);
      #endif
   }
   stage[rank] = MN;  /* let core 0 know we are done */
/*
 * Copy answer back out of workspace and then free workspace
 */
   Mjoin(PATL,gecopy)(rank?mp:mp+mr, N, W, ldw, a, lda);
/*
 * Core 0 waits for all other cores to finish before he frees his work:
 * all non-zero cores access 0's workspace, but 0 does not access others' work
 * after iamax barrier
 */
   if (!rank)
   {
      for (i=1; i < p; i++)
         while(stage[i] != MN);
   }
   free(vp);
}

int Mjoin(PATL,StructIsInitGETF2)(void *vp)
{
   return(((ATL_TGETF2_M_t*)vp)->M);
}

int Mjoin(PATL,tgetf2_nocp)(ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT lda, int *ipiv)
{
   ATL_TGETF2_M_t lu2s[ATL_NTHREADS];
   ATL_INT maxindx[ATL_NTHREADS], stage[ATL_NTHREADS];
   TYPE *works[ATL_NTHREADS];

   ATL_CINT MN = Mmin(M,N);
   ATL_INT p = ATL_NTHREADS, m, mr, i, j;

   if (M < 1 || N < 1)
      return(0);
   m = M / ATL_NTHREADS;
   mr = M - m*ATL_NTHREADS;
/*
 * This logic is necessary since tgetf2 assumes only one processor owns entire
 * logical block.  Can remove if we rewrite tgetf2 to allow the diagonal to
 * span multiple processors
 */
   if (m+mr < N)
   {
      p = M / N;
      if (p)
         m = M / p;
   }
   if (p < 2)   /* not enough rows, call serial algorithm */
      return(Mjoin(PATL,getf2)(M, N, A, lda, ipiv));
   for (i=0; i < p; i++)
   {
      stage[i] = maxindx[i] = -1;
      lu2s[i].M = M;
      lu2s[i].N = N;
      lu2s[i].A = A;
      lu2s[i].lda = lda;
      lu2s[i].ipiv = ipiv;  /* only thread 0 will write ipiv */
      lu2s[i].info = 0;
      lu2s[i].maxindx = maxindx;
      lu2s[i].stage = stage;
      lu2s[i].p = p;
      lu2s[i].rank = i;
      lu2s[i].works = works;
   }
   for (; i < ATL_NTHREADS; i++)
      lu2s[i].M = 0;
   ATL_goparallel(p, Mjoin(PATL,DoWorkGETF2_nowrk), lu2s, NULL);
   return(lu2s[0].info);
}
int Mjoin(PATL,tgetf2)(ATL_CINT M, ATL_CINT N, TYPE *A, ATL_CINT lda, int *ipiv)
{
   ATL_TGETF2_M_t lu2s[ATL_NTHREADS];
   ATL_INT maxindx[ATL_NTHREADS], stage[ATL_NTHREADS];
   TYPE *works[ATL_NTHREADS];

   ATL_CINT MN = Mmin(M,N);
   ATL_INT p = ATL_NTHREADS, m, mr, i, j;

   if (M < 1 || N < 1)
      return(0);
   m = M / ATL_NTHREADS;
   mr = M - m*ATL_NTHREADS;
/*
 * This logic is necessary since tgetf2 assumes only one processor owns entire
 * logical block.  Can remove if we rewrite tgetf2 to allow the diagonal to
 * span multiple processors
 */
   if (m+mr < N)
   {
      p = M / N;
      if (p)
         m = M / p;
   }
   if (p < 2)   /* not enough rows, call serial algorithm */
      return(Mjoin(PATL,getf2)(M, N, A, lda, ipiv));
   for (i=0; i < p; i++)
   {
      stage[i] = maxindx[i] = -1;
      lu2s[i].M = M;
      lu2s[i].N = N;
      lu2s[i].A = A;
      lu2s[i].lda = lda;
      lu2s[i].ipiv = ipiv;  /* only thread 0 will write ipiv */
      lu2s[i].info = 0;
      lu2s[i].maxindx = maxindx;
      lu2s[i].stage = stage;
      lu2s[i].p = p;
      lu2s[i].rank = i;
      lu2s[i].works = works;
   }
   for (; i < ATL_NTHREADS; i++)
      lu2s[i].M = 0;
   ATL_goparallel(p, Mjoin(PATL,DoWorkGETF2), lu2s, NULL);
   return(lu2s[0].info);
}
#ifndef TCPLX
   #undef lda2
#endif
