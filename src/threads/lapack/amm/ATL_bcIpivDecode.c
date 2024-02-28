/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * Copyright (C) 2015 R. Clint Whaley
 */
#include "atlas_misc.h"
#include "atlas_pca.h"
#include "atlas_bcamm.h"
/*
 * Takes ipiv split in (rank,lbn, sbn,sbr) bitpattern, and translates it back
 * to standard getrf global row indices.  See ATL_tipivEncode for more details.
 */
void ATL_bcIpivDecode
(
   ATL_bcpiv_t *bp,/* block-cyclic pivot ptr returned by init */
   int n,          /* number of pivot entries to decode */
   int I           /* ipiv index to start the decoding at */
)
{
   ATL_UINT nprior, k;
   ATL_UINT *ipiv = bp->ipiv + I;
   ATL_CUINT P=bp->R, MU=bp->MU, B=bp->B;
   ATL_CUINT nbSBR=bp->nbSBR, nbRNK=bp->nbRNK, nbSBN=bp->nbSBN;
   ATL_CUINT mskSBR=(1<<nbSBR)-1, mskRNK=(1<<nbRNK)-1, mskSBN=(1<<nbSBN)-1;

   for (k=0; k < n; k++)
   {
      ATL_UINT rank, sbn, sbr, lbn = ipiv[k];
      rank = lbn & mskRNK;
      lbn >>= nbRNK;
      sbr = lbn & mskSBR;
      lbn >>= nbSBR;
      sbn = lbn & mskSBN;
      lbn >>= nbSBN;
/*
 *    Scalapack's indxl2g with srcproc=0
 */
      ipiv[k] = lbn*P*B + rank*B + sbn*MU+sbr;
   }
}
