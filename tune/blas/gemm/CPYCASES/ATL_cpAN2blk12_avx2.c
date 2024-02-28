#include "atlas_misc.h"
#ifdef ALPHAN1
   #define ALPHAN 1
#endif
#ifdef SREAL
   #include <immintrin.h>
   #define ATL_VTYPE __m128
   #define ATL_vgather(v_, p_, i_) v_ = _mm256_i64gather_ps(p_, i_, 4)
   #define ATL_vmskgat(v_, p_, i_, m_) v_ = _mm256_mask_i64gather_ps(z0, p_, i_, m_, 4)
   #define ATL_vuld(v_, p_) v_ = _mm_loadu_ps(p_)
   #define ATL_vmul(d_, s1_, s2_) d_ =  _mm_mul_ps(s1_, s2_)
   #define ATL_vst(p_, v_) _mm_store_ps(p_, v_)
   #define ATL_vzero(v_) v_ = _mm_setzero_ps()
#else
   #define ATL_VLEN 4
   #include "atlas_simd.h"
   #undef ATL_VLEN
   #define ATL_vgather(v_, p_, i_) v_ = _mm256_i64gather_pd(p_, i_, 8)
   #define ATL_vmskgat(v_, p_, i_, m_) \
      v_ = _mm256_mask_i64gather_pd(z0, p_, i_, m_, 8)
#endif
void ATL_USERCPMM(const size_t K, const size_t D, const SCALAR alpha,
                  const TYPE *A, const size_t lda, TYPE *b)
{

   const __m256i off = {0, lda, lda<<1, (lda<<1)+lda};
   const size_t D12=D/12, Dr=D-12*D12, lda4=lda<<2, lda12=(lda4+(lda<<1))<<1;
   const TYPE *A4=A+lda4, *A8=A4+lda4;
   size_t j;
   #ifdef ALPHAN
      const ATL_VTYPE valp = {ATL_rnone,ATL_rnone,ATL_rnone,ATL_rnone};
   #elif defined(ALPHAX)
      const ATL_VTYPE valp = {alpha,alpha,alpha,alpha};
   #endif
   for (j=0; j < D12; j++)
   {
      size_t k;
      for (k=0; k < K; k++)
      {
         ATL_VTYPE a0, a1, a2;
         ATL_vgather(a0, A+k, off);
         ATL_vgather(a1, A4+k, off);
         ATL_vgather(a2, A8+k, off);
         #if defined(ALPHAN) || defined(ALPHAX)
            ATL_vmul(a0, a0, valp);
            ATL_vmul(a1, a1, valp);
            ATL_vmul(a2, a2, valp);
         #endif
         ATL_vst(b, a0);
         ATL_vst(b+4, a1);
         ATL_vst(b+8, a2);
         b += 12;
      }
      A += lda12;
      A4 += lda12;
      A8 += lda12;
   }
   if (Dr)
   {
      const ATL_VTYPE z0={ATL_rzero,ATL_rzero,ATL_rzero,ATL_rzero};
      ATL_VTYPE vmsk;
      const TYPE *pA=A;
      const size_t D4 = Dr>>2, Drr=Dr-(D4<<2), ZPAD = 3 - ((Dr+3)>>2);
      size_t k;
      #ifdef SREAL
         int msk[4];
         #define NONE -1
      #else
         long long int msk[4];
         #define NONE -1L
      #endif
      if (Drr)
      {
         msk[0] = msk[1] = msk[2] = msk[3] = 0;
         switch(Drr)
         {
         case 3:
            msk[2] = NONE;
         case 2:
            msk[1] = NONE;
         case 1:
            msk[0] = NONE;
         default:;
         }
         ATL_vuld(vmsk, (void*)msk);
      }
      #undef NONE
      for (k=0; k < K; k++)
      {
         ATL_VTYPE a0, a1;
         switch(D4)
         {
         case 2:
            ATL_vgather(a0, A+k, off);
            ATL_vgather(a1, A4+k, off);
            #if defined(ALPHAN) || defined(ALPHAX)
               ATL_vmul(a0, a0, valp);
               ATL_vmul(a1, a1, valp);
            #endif
            ATL_vst(b, a0);
            ATL_vst(b+4, a1);
            b += 8;
            pA = A8;
            break;
         case 1:
            ATL_vgather(a0, A+k, off);
            #if defined(ALPHAN) || defined(ALPHAX)
               ATL_vmul(a0, a0, valp);
            #endif
            ATL_vst(b, a0);
            b += 4;
            pA = A4;
         default:;
         }
         if (Drr)
         {
            ATL_vmskgat(a0, pA+k, off, vmsk);
            #if defined(ALPHAN) || defined(ALPHAX)
               ATL_vmul(a0, a0, valp);
            #endif
            ATL_vst(b, a0);
            b += 4;
         }
         for (j=0; j < ZPAD; j++, b += 4)
            ATL_vst(b, z0);
      }
   }
}
