/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1999 R. Clint Whaley
 */
#ifndef ATLAS_KERN3_H
#define ATLAS_KERN3_H

#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_kernel3.h"
#include "atlas_reflevel3.h"
#include "atlas_amm.h"
#include Mstr(Mjoin(ATLAS_PRE,opgen_view.h))
/*
 * Gemm entry points
 */
#ifndef ATL_gemm
   #define ATL_gemm Mjoin(PATL,gemm)
#endif
#define ATL_almm Mjoin(PATL,ammm_aliased_rkK)

#ifdef Left_
   #define Side_ AtlasLeft
   #define SideNM L
#elif defined(Right_)
   #define Side_ AtlasRight
   #define SideNM R
#endif

#ifdef Upper_
   #define Uplo_ AtlasUpper
   #define UploNM U
#elif defined(Lower_)
   #define Uplo_ AtlasLower
   #define UploNM L
#endif

#ifdef UnitDiag_
   #define Unit_ AtlasUnit
   #define UnitNM U
#elif defined(NonUnitDiag_)
   #define Unit_ AtlasNonUnit
   #define UnitNM N
#endif

#ifdef Transpose_
   #define Trans_ AtlasTrans
   #define TransNM T
#elif defined(Notranspose_)
   #define Trans_ AtlasNoTrans
   #define TransNM N
#elif defined(ConjTrans_)
   #define Trans_ AtlasConjTrans
   #define TransNM C
#endif

#if ATL_VWopgen_MAX_KB < 48
   #define ATL_XOVER ATL_VWopgen_MAX_KB
#else
   #define ATL_XOVER 48
#endif
#ifndef TRSM_Xover
   #define TRSM_Xover ATL_XOVER
#endif
#ifndef TRMM_Xover
   #define TRMM_Xover ATL_XOVER
#endif
#ifndef HER2K_Xover
   #define HER2K_Xover ATL_XOVER
#endif
#ifndef SYR2K_Xover
   #define SYR2K_Xover ATL_XOVER
#endif
#ifndef HERK_Xover
   #define HERK_Xover ATL_XOVER
#endif
#ifndef SYRK_Xover
   #define SYRK_Xover ATL_XOVER
#endif
#ifndef HEMM_Xover
   #define HEMM_Xover ATL_XOVER
#endif
#ifndef SYMM_Xover
   #define SYMM_Xover ATL_XOVER
#endif

#endif
