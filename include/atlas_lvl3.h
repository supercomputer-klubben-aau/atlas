/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */

#ifndef ATLAS_LVL3_H
#define ATLAS_LVL3_H

#include "atlas_misc.h"
#include "atlas_amm.h"
#include "atlas_f77.h"
#include "atlas_level3.h"
void *Mjoin(PATL,utrsmR_alloc)
   (sminfo_t *ip, int N, TYPE **Diag, TYPE **L, TYPE **R, TYPE **w);
void *Mjoin(PATL,utrsmL_alloc)
   (sminfo_t*ip, int N, TYPE **Diag, TYPE **L, TYPE **R, TYPE **w);
void* Mjoin(PATL,trsmR_alloc) (const ipinfo_t *gip, sminfo_t *tip,
      ATL_CSZT N, ATL_CSZT *Tsz, ATL_CSZT *Rsz, TYPE **Diag, TYPE **L,
      TYPE **R, TYPE **w);
void* Mjoin(PATL,trmmR_alloc) (const ipinfo_t *gip, tminfo_t *tip,
      ATL_CSZT N, ATL_CSZT *Tsz, ATL_CSZT *Rsz, ATL_CSZT *Csz, TYPE **bw,
      TYPE **L, TYPE **R, TYPE **w);
void* Mjoin(PATL,trsmL_alloc) (const ipinfo_t *gip, sminfo_t *tip,
      ATL_CSZT N, ATL_CSZT *Tsz, ATL_CSZT *Rsz, TYPE **Diag, TYPE **L,
      TYPE **R, TYPE **w);
void* Mjoin(PATL,trmmL_alloc) (const ipinfo_t *gip, tminfo_t *tip,
      ATL_CSZT N, ATL_CSZT *Tsz, ATL_CSZT *Rsz, ATL_CSZT *Csz, TYPE **bw,
      TYPE **L, TYPE **R, TYPE **w);
void Mjoin(PATL,tminfoR_UT) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoR_UN) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoR_LT) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoR_LN) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
#ifdef TCPLX
void Mjoin(PATL,tminfoR_UC) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoR_LC) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
#endif
   int Mjoin(PATL,trsmR_LTUN)(ipinfo_t *gip, sminfo_t *tip,
         const enum ATLAS_DIAG Diag, ATL_CSZT R, ATL_CSZT N,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *X, ATL_CSZT ldx,
         ATL_CSZT Tsz, ATL_CSZT Rsz, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w );
   int Mjoin(PATL,trmmR_LTUN)(ipinfo_t *gip, tminfo_t *tip,
         ATL_CSZT M, ATL_CSZT N, const SCALAR alpha, const TYPE *A,
         ATL_CSZT lda, TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz,
         ATL_CSZT Csz, TYPE *tbw, TYPE *L, TYPE *R, TYPE *w );
   int Mjoin(PATL,trsmR_LNUT)(ipinfo_t *gip, sminfo_t *tip,
         const enum ATLAS_DIAG Diag, ATL_CSZT R, ATL_CSZT N,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *X, ATL_CSZT ldx,
         ATL_CSZT Tsz, ATL_CSZT Rsz, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w );
   int Mjoin(PATL,trmmR_LNUT)(ipinfo_t *gip, tminfo_t *tip,
         ATL_CSZT M, ATL_CSZT N, const SCALAR alpha, const TYPE *A,
         ATL_CSZT lda, TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz,
         ATL_CSZT Csz, TYPE *tbw, TYPE *L, TYPE *R, TYPE *w );
void Mjoin(PATL,tminfoL_UT) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoL_UN) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoL_LT) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoL_LN) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
#ifdef TCPLX
void Mjoin(PATL,tminfoL_UC) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
void Mjoin(PATL,tminfoL_LC) (tminfo_t *tip, ipinfo_t *gip,
      const enum ATLAS_DIAG Diag,
      ATL_CSZT N, ATL_CSZT R, const SCALAR alpha, ATL_CSZT lda, ATL_CSZT ldb);
#endif
   int Mjoin(PATL,trsmL_LTUN)(ipinfo_t *gip, sminfo_t *tip,
         const enum ATLAS_DIAG Diag, ATL_CSZT R, ATL_CSZT N,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *X, ATL_CSZT ldx,
         ATL_CSZT Tsz, ATL_CSZT Rsz, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w );
   int Mjoin(PATL,trmmL_LTUN)(ipinfo_t *gip, tminfo_t *tip,
         ATL_CSZT M, ATL_CSZT N, const SCALAR alpha, const TYPE *A,
         ATL_CSZT lda, TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz,
         ATL_CSZT Csz, TYPE *tbw, TYPE *L, TYPE *R, TYPE *w );
   int Mjoin(PATL,trsmL_LNUT)(ipinfo_t *gip, sminfo_t *tip,
         const enum ATLAS_DIAG Diag, ATL_CSZT R, ATL_CSZT N,
         const SCALAR alpha, const TYPE *A, ATL_CSZT lda, TYPE *X, ATL_CSZT ldx,
         ATL_CSZT Tsz, ATL_CSZT Rsz, TYPE *diag, TYPE *L, TYPE *RW, TYPE *w );
   int Mjoin(PATL,trmmL_LNUT)(ipinfo_t *gip, tminfo_t *tip,
         ATL_CSZT M, ATL_CSZT N, const SCALAR alpha, const TYPE *A,
         ATL_CSZT lda, TYPE *X, ATL_CSZT ldx, ATL_CSZT Tsz, ATL_CSZT Rsz,
         ATL_CSZT Csz, TYPE *tbw, TYPE *L, TYPE *R, TYPE *w );


#endif
