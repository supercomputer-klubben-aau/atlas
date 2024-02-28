#ifndef ATLAS_TBITVEC_H
   #define ATLAS_TBITVEC_H
   #include "atlas_bitvec.h"
   #include "atlas_cbc.h"
   #include "atlas_tprim.h"

/*
 * Global bitvecs distribute the BV amongst at most P areas that are all
 * accessed with separate locks, minimizing contention during normal access.
 * The accessor functions all take a rank, which will be used to select the
 * native BV to scope first.  Only when the native BV is empty will we look
 * at other's BVs.  This gives us work stealing while minimizing contention.
 * The main source of contention is from false sharing, which is hugely
 * magnificed with bitvecs.  To prevent this, we expand all independent accessor
 * areas to multiples of a safe linesize (currently 128bytes, set in
 * atlas_tprim.h).  The data structure then looks like, each rounded to 128:
 * 1. summary area:
 *     0     1        2       3      4        5          6       7
 *    <P> <gnbits> <sumSz> <lckSz> <bvSz> <nLrgBlks> <nSmlBlks> <b>
 * 2. local lock storage area : *each lock* rounded up to 128 bytes
 * 3. locBV area: <nunset> [serial BV] -> <nleft> <nbits> <BV> => 128 rounded
 * -> <b> is the number of bits in a small block, big block as b+1 bits
 */
#define ATL_TBV_P      0
#define ATL_TBV_GNBITS 1
#define ATL_TBV_SUMSZ  2
#define ATL_TBV_LCKSZ  3
#define ATL_TBV_BVSZ   4
#define ATL_TBV_NLGB   5
#define ATL_TBV_NSMB   6
#define ATL_TBV_B      7
#define ATL_TBV_NUNSET 8
#ifndef ATL_UL_t
   #define ATL_UL_t unsigned long
#endif
#ifndef ATL_CUL_t
   #define ATL_CUL_t const unsigned long
#endif
size_t ATL_tSizeofBV(unsigned long nbits, unsigned int P);
void ATL_tInitBV(void *vp, unsigned long nbits, unsigned int P);
void *ATL_tNewBV(unsigned long nbits, unsigned int P);
void ATL_tFreeBV(void *bv);
unsigned long ATL_tGetTotBitsBV(void *bv); /* global number bits in bv */
long ATL_tInfoBV(const void *vp, long what); /* what is TBV_* above */
void ATL_tPrintBV(FILE *fpout, char *nm, void *bv);
/*
 * if (flg&1): no need to set locks before changing
 */
void ATL_tReverseAllBitsBV(void *vp, unsigned int flg);
/*
 * sets *lastbit to the last global bit "owned" by thread rank, and
 * RETURNS: first global bit owned by thread rank.
 */
unsigned long ATL_tGetLocalBoundsBV(void *bv, unsigned int rank,
                                    unsigned long *lastbit);
/*
 * Read and optionally change bit at position pos, return old value
 */
long ATL_tScopeBitBV(void *vp, unsigned long bit, unsigned int flg);
#define ATL_tIsBitSetBV(bv_, p_) ATL_tScopeBitBV(bv_, p_, 513)
#define ATL_tSetBitBV(bv_, p_) ATL_tScopeBitBV(bv_, p_, 2)
#define ATL_tUnsetBitBV(bv_, p_) ATL_tScopeBitBV(bv_, p_, 4)
/*
 * Combines all P individual global bitvecs into one local bitvec.
 * Mutexes are not locked, so info may be wrong.
 * RETURNS: number of unset bits in new local bitvec
 */
long ATL_tGlb2locBV(ATL_BV_t *lbv, void *gbv, unsigned long pos);
/*
 * On input *nbits is the max number of bits in bv to set to the bit pattern
 * in the least significant bits of mask.  On exit, the number of bits actually
 * set (may be less if some bits owned by another mutex; only one lock will
 * be made).
 * RETURNED: original bitpattern, *nbits set to the number of bits
 */
unsigned long ATL_tSetRangeBV(void *bv, unsigned int *nbits, unsigned long pos,
                              unsigned long mask);
/*
 * Find first unset bit (starting at rank's block), w/o locking (unsafe)
 */
long ATL_tFindUnsetBitBV(void *vp, unsigned int rank);
/*
 * These functions find a [set,unset] bit, and complement it.
 * They return the global bit position that was changed.  They first look
 * at the bits "owned" by the given rank, and if those bits are all [unset,set],
 * they look at other thread's local BV, and try to change those bits.
 * Thus, this method can be used for work stealing, for instance.  These
 * functions may need to await locks, and so they can be delayed in returning
 * while they wait on locks held by other threads.
 * A return of -1 indicates that all bits are [unset,set].
 */
long ATL_tSetUnsetBitBV(void *bv, unsigned int rank);


#endif
