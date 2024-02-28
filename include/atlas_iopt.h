/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2016 R. Clint Whaley
 */
#ifndef ATLAS_IOPT_H
   #define ATL_IOPT_H

#define ATL_IsPow2(i_) /* 1 iff i_ is a power-of-2, 0 else */ \
   ( !((i_)&((i_)-1)) )

/*
 * ans_ undefined if i_=0 (no bits set)
 */
#if defined(__GNUC__) && defined(ATL_GAS_x8664) || defined(ATL_GAS_x8632)
   #define ATL_FAST_LSB 1  /* signal we've got a fast least set bit */
   #define ATL_iLeastSetBit(ans_, i_) \
      __asm__("bsf %1, %%eax" : "=a"(ans_):"r"(i_))
#elif defined(__GNUC__) && defined(ATL_GAS_ARM64)
   #define ATL_FAST_LSB 1
   #define ATL_iLeastSetBit(ans_, i_) \
      __asm__("rbit %0, %0; clz %0, %0" : "=r"(ans_):"0"(i_))

#else
   #define ATL_iLeastSetBit(ans_, i_) \
   { \
      register unsigned int p_, k_=(i_), d_=16, msk_=(1<<16)-1, rej_; \
      /* Lower 16 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ = rej_; \
      k_ >>= rej_; \
      d_ >>= 1; \
      msk_ >>= d_; \
      /* Lower 8 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ += rej_; \
      k_ >>= rej_; \
      d_ >>= 1; \
      msk_ >>= d_; \
      /* Lower 4 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ += rej_; \
      k_ >>= rej_; \
      d_ >>= 1; \
      msk_ >>= d_; \
      /* Lower 2 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ += rej_; \
      k_ >>= rej_; \
      ans_ = p_ + 1 - (k_&1); \
   }
#endif
   #define ATL_iLeastSetBit8(ans_, i_) /* set bit known to be 1st 8 bits */ \
   { \
      register unsigned int p_, k_=(i_), d_=4, msk_=(1<<4)-1, rej_; \
      /* Lower 4 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ = rej_; \
      k_ >>= rej_; \
      d_ >>= 1; \
      msk_ >>= d_; \
      /* Lower 2 */ \
      rej_ = (k_&msk_) ? 0 : d_; /* # of lower bits to reject */ \
      p_ += rej_; \
      k_ >>= rej_; \
      d_ >>= 1; \
      msk_ >>= d_; \
      ans_ = p_ + 1 - (k_&1); \
   }
   /* set bit known to be in 1st 8 bits, or behavior undefined */
   #define ATL_iLeastSetBit4(ans_, i_) \
   { \
      register unsigned int p_, k_=(i_); \
      /* Lower 2 */ \
      p_ = (k_&0x3) ? 0 : 2; /* # of lower bits to reject */ \
      k_ >>= p_; \
      ans_ = p_ + 1 - (k_&1); \
   }
#endif /* end multiple include guard */

#ifdef ATL_WANT_ILCM
   #define ATL_WANT_IGCD
#endif

#if defined(ATL_WANT_IGCD) && !defined ATL_GCD_H
   #define ATL_GCD_H

#ifndef INLINE
   #if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
      #ifdef __GNUC__
         #define INLINE __inline__
      #else
         #define INLINE inline
      #endif
   #else
      #define INLINE
   #endif
#endif
#ifndef Mmax
   #define Mmax(i_, j_) ((i_ >= j_) ? i_:j_)
#endif
#ifndef Mmin
   #define Mmin(i_, j_) ((i_ >= j_) ? j_:i_)
#endif

static INLINE unsigned int ATL_iGCD(const unsigned int M, const unsigned int N)
{
/*
 * Compute the Greatest Common Divisor using Stein's Algorithm.  Has been
 * modified to exploit bit-level operations, with inline assembly for x86
 * and ARM64.  On these systems, seems to be faster than normal Stein's
 * always.  If you don't have a FAST_LSB, however, sometimes its faster than
 * normal Stein's, and sometimes not.  On systems where its slower, should
 * be easy to tweak ATL_iLeastSetBit to at least tie.
 */
   unsigned int gcdsh=0;
   if (M != N)
   {
      unsigned int tmp, min, max;
      max = Mmax(M,N);
      min = Mmin(M,N);
      if (ATL_IsPow2(max))                /* power of 2s in min are whole GCD */
      {
         ATL_iLeastSetBit(tmp, min);      /* compute shift to make min odd */
         return(1<<tmp);
      }
      if (ATL_IsPow2(min))               /* Min pwr2 between max&min are GCD */
      {
         unsigned int sh;
         ATL_iLeastSetBit(tmp, max);      /* compute shift to make max odd */
         ATL_iLeastSetBit(sh,  min);      /* compute shift to make min odd */
         gcdsh = Mmin(tmp, sh);
         return(1<<gcdsh);
      }
      tmp = min;
      if (!(max&1))                        /* if max isn't odd, make it so */
      {                                    /* so we match end-of-loop state */
         unsigned int sh;
         int shA, shI;                     /* shift mAx, shift mIn */

         ATL_iLeastSetBit(shI, min);       /* compute shift to make min odd */
         min >>= shI;                      /* all 2s in min removed, now odd */
         ATL_iLeastSetBit(shA, max);       /* compute shift to make max odd */
         sh = Mmin(shA, shI);
         max >>= sh;                       /* max is divided by common pwr 2 */
         gcdsh += sh;
         if (min == 1)                     /* easily predicted not true */
            goto MIN_POW2;                 /* hint comp put ret outside blk */
         do
         {
            max -= (max&1) ? min:0;        /* make max even */
            max >>= 1;                     /* max now even or odd */
         }
         while (max >= min);
         tmp = max;
         max = min;                        /* max=min, so odd */
         min = tmp;
      }
      if (tmp)
      {
         do                                /* max>=min, max is known odd */
         {
            int shI;                  /* shift mIn */
/*
 *          Problem with branches on even/odd is that they are unpredictable.
 *          If we assume random scattering of 1s, wrong 50% of time.  However,
 *          for random, chance of 4 0s in a row is 6.25%, so this is a
 *          predictable branch for non-sparse min.  We will remove all
 *          contiguous 0s at once, sparse min will rapidly be reduced.
 */
            #ifndef ATL_FAST_LSB
            if (min&0xF)
               ATL_iLeastSetBit4(shI, min) /* compute shift to make min odd */
            else
            #endif
               ATL_iLeastSetBit(shI, min); /* compute shift to make min odd */
            min >>= shI;                   /* all 2s in min removed, now odd */
            if (min == 1)                  /* easily predicted not true */
               goto MIN_POW2;              /* hint comp put ret outside blk */
            do
            {
               max -= (max&1) ? min:0;     /* make max even */
               max >>= 1;                  /* max now even or odd */
            }
            while (max >= min);
            tmp = max;
            max = min;                     /* max=min, so odd */
            min = tmp;
         }
         while(tmp);                       /* max is known to be odd */
      }
      return(max<<gcdsh);
   }
   return(M);
MIN_POW2:
   return(1<<gcdsh);
}
#endif /* end GCD multiple include guard */

#if defined(ATL_WANT_ILCM) && !defined ATL_LCM_H
   #define ATL_LCM_H
static INLINE unsigned int ATL_iLCM(unsigned int M, unsigned int N)
/*
 * Comnputes LCM = M*N / GCD;
 */
{
   unsigned long long mul;
   int gcd;
   gcd = ATL_iGCD(M,N);
   if (gcd == ((M>=N) ? M:N))
      return(gcd);
   mul = M * N;
   #ifdef ATL_FAST_LSB  /* avoid division if LSB is one instruction */
   if (ATL_IsPow2(gcd))
   {
      ATL_iLeastSetBit(gcd, gcd);
      gcd = (mul >> gcd);
   }
   else
   #endif
      gcd = mul / gcd;
   return(gcd);
}
   #undef INLINE
#endif /* end LCM multiple include guard */

