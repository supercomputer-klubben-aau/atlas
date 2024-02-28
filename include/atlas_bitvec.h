#ifndef ATLAS_BITVEC_H
   #define ATLAS_BITVEC_H
   #include "atlas_misc.h"
   #include "atlas_type.h"

/*
 * Macros for number of bits in array entry, shift value to div/mul by nbits,
 * mask to do modulo(bpiBV), and all bits set
 */
#if ATL_LSIZE == 8
   #define bpiBV 64  /* # bits per integral type */
   #define shBV 6
   #define modmskBV 0x3FL
   #define allsetBV 0xFFFFFFFFFFFFFFFFL
#elif ATL_LSIZE == 4
   #define bpiBV 32  /* # bits per integral type */
   #define shBV 5
   #define modmskBV 0x1F
   #define allsetBV 0xFFFFFFFF
#else
   #error "LONG is neither 8 or 4 bytes!"
#endif
#define ATL_BV_t unsigned long
/*
 * A bitvector consists of an integral array.  The first element of the
 * array indicates the number bits stored in the bitvec array.
 */
#define ATL_FreeBV(bv_) free(bv_)  /* bitvec just an integral array */
ATL_BV_t *ATL_NewBV(unsigned long nbits);
void ATL_InitBV(unsigned long nbits, ATL_BV_t *bv, ATL_BV_t msk);
ATL_BV_t *ATL_ExpandBV(ATL_BV_t *bv, unsigned long newbits);
unsigned long ATL_GetTotBitsBV(ATL_BV_t *bv);
ATL_BV_t ATL_SetEltBV
   (ATL_BV_t *bv, unsigned long ielt, ATL_BV_t val);
ATL_BV_t ATL_GetEltBV(ATL_BV_t *bv, unsigned long ielt);
int ATL_SetBitBV(ATL_BV_t *bv, unsigned long pos);
int ATL_UnsetBitBV(ATL_BV_t *bv, unsigned long pos);
int ATL_IsBitSetBV(ATL_BV_t *bv, unsigned long bpos);
long ATL_FindFirstSetBitBV(ATL_BV_t *bv, unsigned long bs);
long ATL_FindFirstUnsetBitBV(ATL_BV_t *bv, unsigned long bs);
void ATL_ReverseAllBitsBV(ATL_BV_t *bv);
void ATL_SetAllBitsBV(ATL_BV_t *bv);
void ATL_UnsetAllBitsBV(ATL_BV_t *bv);
int ATL_IsBitRangeSetBV(ATL_BV_t *bv, unsigned long b0, unsigned long bN);
int ATL_IncorpBV(ATL_BV_t *dst, ATL_BV_t *src, unsigned long pos);

#endif
