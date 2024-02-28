#include "atlas_simd.h"
#ifdef SREAL
   #define ATL_vmskld(v_, p_, m_) v_ = _mm256_maskload_ps(p_, m_)
   #define ATL_vmskst(p_, v_, m_) _mm256_maskstore_ps(p_, m_, v_)
   #define TYPE float
   #define BSZ 16
#else
   #define ATL_vmskld(v_, p_, m_) v_ = _mm256_maskload_pd(p_, m_)
   #define ATL_vmskst(p_, v_, m_) _mm256_maskstore_pd(p_, m_, v_)
   #define TYPE double
   #define BSZ 12
#endif
#ifdef BETAN1
   #define BETAN
#endif
#ifdef ALPHAN1
   #define ALPHAN
#endif
/* HERE vl=4, mu=3, nu=3 bs=12 */
void ATL_USERCPMM
(
   const size_t M,      /* number of rows in A */
   const size_t N,      /* number of columns in A */
   const TYPE alpha,    /* scalar for b */
   const TYPE *b,       /* matrix stored in 3x3-major order */
   const TYPE beta,     /* scalar for C */
   TYPE *C,             /* matrix to be copied to access-major format */
   const size_t ldc     /* stride between row elements */
)
{
   const unsigned int mf = M/3, nf = N/3;
   const unsigned int m = mf*3, n = nf*3, mr = M-m, nr = N-n;
   unsigned int pansz = BSZ;
   const size_t incC0 = (ldc+1)*3;
   unsigned int i, j;
   TYPE *C0=C, *C1=C0+ldc, *C2=C1+ldc;
   #ifdef SREAL
      const __m256i m0={-1L,0xFFFFFFFFL,0L,0L},
         m1={0xFFFFFFFF00000000L,0xFFFFFFFFL,0L,0L}, m2={0L, 0xFFFFFFFFL,0L,0L};
   #else
      const __m256i m0={-1,-1,-1,0}, m1={ 0,-1,-1,0}, m2={ 0, 0,-1,0};
   #endif
   #if defined(ALPHAN) || defined(BETAN)
      TYPE NONE=(-1.0);
   #endif
   #if defined(ALPHAX) || defined(ALPHAN)
      ATL_VTYPE valpha;
   #endif
   #if defined(BETAX) || defined(BETAN)
      ATL_VTYPE vbeta;
   #endif

   #ifdef BETAN
      ATL_vbcast(vbeta, &NONE);
   #elif defined(BETAX)
      ATL_vbcast(vbeta, &beta);
  #endif
   #ifdef ALPHAN
      ATL_vbcast(valpha, &NONE);
   #elif defined(ALPHAX)
      ATL_vbcast(valpha, &alpha);
  #endif

   for (j=0; j < nf; j++)
   {
      const TYPE *p = b;
      unsigned int psz = pansz+BSZ, incC = incC0 - (mf-j)*3;
      ATL_VTYPE c0, b0;

      #ifdef BETA0
         ATL_vmskld(c0, p+0, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+0, m0);
         ATL_vmskld(b0, C0, m0);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C0, m0);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+0, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C0, c0, m0);
      C0 += 3;
      #ifdef BETA0
         ATL_vmskld(c0, p+3, m1);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+3, m1);
         ATL_vmskld(b0, C1, m1);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C1, m1);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+3, m1);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C1, c0, m1);
      C1 += 3;
      #ifdef BETA0
         ATL_vmskld(c0, p+6, m2);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+6, m2);
         ATL_vmskld(b0, C2, m2);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C2, m2);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+6, m2);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C2, c0, m2);
      C2 += 3;
      p += pansz;
      for (i=j+1; i < mf; i++, p += psz, psz += BSZ)
      {
      #ifdef BETA0
         ATL_vmskld(c0, p+0, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+0, m0);
         ATL_vmskld(b0, C0, m0);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C0, m0);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+0, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C0, c0, m0);
      C0 += 3;
      #ifdef BETA0
         ATL_vmskld(c0, p+3, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+3, m0);
         ATL_vmskld(b0, C1, m0);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C1, m0);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+3, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C1, c0, m0);
      C1 += 3;
      #ifdef BETA0
         ATL_vmskld(c0, p+6, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(c0, c0, valpha);
         #endif
      #elif defined(ALPHA1) && (defined(BETAX) || defined(BETAN))
         ATL_vmskld(c0, p+6, m0);
         ATL_vmskld(b0, C2, m0);
         ATL_vmac(c0, b0, vbeta);
      #else
         ATL_vmskld(c0, C2, m0);
         #if defined(BETAX) || defined(BETAN)
            ATL_vmul(c0, c0, vbeta);
         #endif
         ATL_vmskld(b0, p+6, m0);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmac(c0, b0, valpha);
         #else
            ATL_vadd(c0, c0, b0);
         #endif
      #endif
      ATL_vmskst(C2, c0, m0);
      C2 += 3;
      }
      switch(mr)
      {
      case 1:
            C0[0] = beta*C0[0] + alpha*p[0];
            C1[0] = beta*C1[0] + alpha*p[3];
            C2[0] = beta*C2[0] + alpha*p[6];
         break;
      case 2:
            C0[0] = beta*C0[0] + alpha*p[0];
            C0[1] = beta*C0[1] + alpha*p[1];
            C1[0] = beta*C1[0] + alpha*p[3];
            C1[1] = beta*C1[1] + alpha*p[4];
            C2[0] = beta*C2[0] + alpha*p[6];
            C2[1] = beta*C2[1] + alpha*p[7];
         break;
      default:;
      }
      C0 += incC;
      C1 += incC;
      C2 += incC;
      pansz += BSZ;
      b += pansz;
   }
   switch(nr)
   {
      const TYPE *p;
   case 1:
      p = b;
            C0[0] = beta*C0[0] + alpha*p[0];
      break;
   case 2:
      p = b;
            C0[0] = beta*C0[0] + alpha*p[0];
            C0[1] = beta*C0[1] + alpha*p[1];
            C1[1] = beta*C1[1] + alpha*p[4];
      break;
   default:;
   }
}
