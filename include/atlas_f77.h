/*
 * Automatically Tuned Linear Algebra Software v3.11.41
 * (C) Copyright 1997 R. Clint Whaley
 */
#ifndef ATLAS_F77_H
#define ATLAS_F77_H

   #ifndef ATL_F77_SUBROUTINE
      #define ATL_F77_SUBROUTINE void
   #endif
   #ifndef F77_INTEGER
      #define F77_INTEGER int
   #endif
   #if defined(CRAY)
      #define UseTransChar 1
      #include <fortran.h>
      #define F77_CHAR _fcd
      #define ATL_F2C_TransChar(c) (*(_fcdtocp(c) ))
      #define ATL_C2F_TransChar(c) (_cptofcd(&(c), 1))
   #elif defined(StringStructVal)
      typedef struct {char *cp; F77_INTEGER len;} F77_CHAR;
      #define ATL_F2C_TransChar(c) (*(c.cp))
      #define UseTransChar 2
   #elif defined(StringStructPtr)
      typedef struct {char *cp; F77_INTEGER len;} F77_CHAR;
      #define ATL_F2C_TransChar(c) (*(c->cp))
      #define UseTransChar 3
   #else
      #define ATL_DeclareSlens
      #define F77_CHAR char *
      #define ATL_F2C_TransChar(c) (*(c))
      #define ATL_C2F_TransChar(c) (&(c))
      #define ATL_STRLEN_1 ,F77_INTEGER ATL_Slen1
      #define ATL_STRLEN_2 ,F77_INTEGER ATL_Slen1, F77_INTEGER ATL_Slen2
      #define ATL_STRLEN_3 ,F77_INTEGER ATL_Slen1, F77_INTEGER ATL_Slen2, \
                           F77_INTEGER ATL_Slen3
      #define ATL_STRLEN_4 ,F77_INTEGER ATL_Slen1, F77_INTEGER ATL_Slen2, \
                           F77_INTEGER ATL_Slen3, F77_INTEGER ATL_Slen4
      #define ATL_STRLEN_1_para ,ATL_Slen1
      #define ATL_STRLEN_2_para ,ATL_Slen1, ATL_Slen2
      #define ATL_STRLEN_3_para ,ATL_Slen1, ATL_Slen2, ATL_Slen3
      #define ATL_STRLEN_4_para ,ATL_Slen1, ATL_Slen2, ATL_Slen3, ATL_Slen4
   #endif

   #ifndef ATL_STRLEN_1
      #define ATL_STRLEN_1
      #define ATL_STRLEN_2
      #define ATL_STRLEN_3
      #define ATL_STRLEN_4
      #define ATL_STRLEN_1_para
      #define ATL_STRLEN_2_para
      #define ATL_STRLEN_3_para
      #define ATL_STRLEN_4_para
   #endif

#endif
