/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#ifndef ATLAS_ENUM_H
   #define ATLAS_ENUM_H

   #define CBLAS_ENUM_ONLY
   #include "cblas.h"
   #undef CBLAS_ENUM_ONLY

   #define ATLAS_ORDER CBLAS_ORDER
      #define AtlasRowMajor CblasRowMajor
      #define AtlasColMajor CblasColMajor
   #define ATLAS_TRANS CBLAS_TRANSPOSE
      #define AtlasNoTrans CblasNoTrans
      #define AtlasTrans CblasTrans
      #define AtlasConjTrans CblasConjTrans
   #define ATLAS_UPLO CBLAS_UPLO
      #define AtlasUpper CblasUpper
      #define AtlasLower CblasLower
   #define ATLAS_DIAG CBLAS_DIAG
      #define AtlasNonUnit CblasNonUnit
      #define AtlasUnit CblasUnit
   #define ATLAS_SIDE CBLAS_SIDE
      #define AtlasLeft  CblasLeft
      #define AtlasRight CblasRight

#endif

