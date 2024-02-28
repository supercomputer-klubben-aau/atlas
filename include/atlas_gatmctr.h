#ifndef ATLAS_GATMCTR_H
   #define ATLAS_GATMCTR_H

#include "atlas_bitvec.h"
#include "atlas_threads.h"
#include "atlas_type.h"
/*
 * Define bit positions in the flag bitvec of boolean flag
 * NOTE: POLL implies LAST1
 */
#define ATL_GAC_MOV   1  /* set: move numbers priv->pub */
#define ATL_GAC_NOPUB 2  /* set: make all data private (static scheduling) */
#define ATL_GAC_NOPRV 4  /* set: make all data public */
#ifndef  GAC_FLAG_SET
   #define GAC_FLAG_SET(field_, bit_) ( ((field_) & (bit_)) != 0 )
#endif
#define ATL_GAC_PUB ATL_GAC_NOPRV
#define ATL_GAC_PRV ATL_GAC_NOPUB
#define ATL_GAC_MIX ATL_GAC_MOV
/*
 * Define offsets for initial SAFELS of global part of data structure
 */
#define ATL_GAC_P   0     /* these defines used to index */
#define ATL_GAC_FLG 1     /* data describing the gatmctr */
#define ATL_GAC_INC 2     /* only used in mixed case */
#define ATL_GAC_N   3     /* only used in mixed case */

/*
 * Global atomic counters return number between [1,N], with 0 meaning no
 * more count.  Count range is split amongst rank-linked counters, and so
 * is affectively random.  However, NOPRV always returns 1 as last job #.
 * Below, off is offset (negative) to get to the unaligned ptr to free.
 * (1) NOPRV: Perfect work stealing, but high overhead.  Last call returns 1.
 *            Use for large jobs that dominate mutex cost, and perfect stealing
 *            is very helpful.
 * -> [<P><flg><incLow><incCnt><off>][<low1>...<lowP>][<act1>]...[<actP>][actL]
 * (2) NOPUB: Pure static scheduling, very low overhead.  Use for very
 *            small jobs, where limited dynamic recovery enabled by having
 *            threads finish their portion, and re-enroll as a new rank.
 * -> data: [<P><flg><off>][<cnt1><low1>]....[<cntP><lowP>]
 * (3) else : mixed case with private and public jobs.  Low average overhead,
 *            but some threads may return 0 before the job is finished.
 *            Useful only when you've got another batch of jobs to do after
 *            this one.  See below for how work stealing works, and struct.
 * NOTE: only (1) above guarantees that the last non-zero return is 1.
 *       (2) & (3) must assume random count return.
 *
 * MIXED PUBLIC AND PRIVATE COUNTERS:
 * The total range is [low,N], low >= 1, N > low, so totSpan = N-low+1.
 * Public jobs are in range [lowPub,N] and count down (PubSpan=cntPub-lowPub+1,
 * where lowPub initialized to N).  Private jobs are in range [low,lowPub-1]
 * and count up (PrvSpan=lowPub-cntPrv, where cntPrv initialized to low).
 * This look like:
 *     ..Private count..|..Public Count..
 *      ---------------> <---------------
 *     [cntPrv,LowPub-1]|[lowPub, cntPub]
 *
 * -> private jobs consumed by cntPrv++, public jobs by cntPub--,
 *    with cntPrv initialized to low, and cntPub to N.
 * -> Work is moved from private to public by subtracting from LowPub, which
 *    increases the public span while decreasing the private span.
 * -> If public work taken since last access (lastPubcnt != PubCnt), move
 *    private range to public
 * -> Start wt P jobs in public, rest private, except for small problems.
 * -> Once cntPub starts to change, private count will be moved to public
 *    to keep number of public jobs constant.
 *
 * SMALL PROBLEMS:
 * ===============
 * For very small problems, dynamic scheculing has too much overhead, so we
 * set the public count to 0 (cntPub=0).  This gives us counters that
 * implement pure static scheduling for P threads.  If the last core wakes
 * up terribly late, this can cause horrible delays, but this can be
 * ameliorated by having finished threads rejoin the work with a new virtual
 * rank.  Thus, in the best case, small problems have no scheduling overhead,
 * and in worst (where scheduling would die), we have only the cost of
 * re-enrolling, so something like: (job cost) + (mutex cost).
 *
 * DATA STRUCTURE EXPLANATION
 * ==========================
 * Keeping independent items on separate cache lines prevents false sharing,
 * which is critical for parallel performance (false sharing on locks or
 * data can cause safe access times to increase by orders of magnitude).
 * ATLAS therefore defines a macro ATL_SAFELS which is assumed >= actual
 * line size.  Presently it is set at 128 bytes. The data structure is then
 * manually laid out to minimize false sharing, where each [] is SAFELS region.
 * The global description is:
 *   [<P><flag><inc><N><off>]
 * The per-core information is:
 *    P*([<cntPrv><lastPubCnt>][lowPub][<cntPub>][<lock>]) :
 *  totSz = SAFELS + P*inc, inc=3*SAFELS+CEIL(sizeof(lock)/SAFELS),
 */

size_t ATL_gatmctr_sizeof(unsigned int P, long flg);
void *ATL_gatmctr_init(void *ac, unsigned int P, long cnt, long flg);
void *ATL_gatmctr_alloc(unsigned int nctr, unsigned int P, long cnt, long flg,
                        long *inc);
void ATL_gatmctr_initByRank(void *ac, unsigned int rank, long lowPrv,
                            long npriv, long npub);
void *ATL_gatmctr_new(unsigned int P, long N, long flg);
void ATL_gatmctr_free(void *ac);
void ATL_gatmctr_print(FILE *fpout, void *ac);
long ATL_gatmctr_dec(void *ac, unsigned int rank);

#endif  /* end ifndef ATLAS_GATMCTR_H */
