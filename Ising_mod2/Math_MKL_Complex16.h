//
//  Math_MKL_Complex16.h
//  Ising_mod2
//
//  Created by Pablo BalRes on 26/01/2020.
//  Copyright © 2020 Pablo BalRes. All rights reserved.
//

#ifndef Math_MKL_Complex16_h
#define Math_MKL_Complex16_h

#include <stdio.h>
#include <math.h>
#include "mkl.h"


MKL_Complex16 Conju_MKL(MKL_Complex16 V);

MKL_Complex16 multipl_conj_MKL(MKL_Complex16 V,MKL_Complex16 W);

MKL_Complex16 Sum_MKL(MKL_Complex16 V, MKL_Complex16 W);

//Diagonalize a general matrix
void Diag_Matrix(MKL_Complex16 **Ham, int NumTot,double *E);

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda );

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

#endif /* Math_MKL_Complex16_h */




//#include "mkl_types.h"

//#include "mkl.h"
/*
 #ifndef lapack_complex_double
 #define lapack_complex_double   MKL_Complex16
 #endif
 
 #define MKL_Complex16 _Complex
 */


//typedef struct{ double re; double im; } MKL_Complex16;
//#define MKL_Complex16 std::complex
//#define MKL_Complex16 std::<double>complex
//#define MKL_Complex16 <double>complex
//typedef MKL_Complex16 complex  ;
// MKL_Complex16 ---> double complex
//typedef double<complex> Dcomplex;
//typedef _Complex double DCOMPLEX;
//#define MKL_Complex16 DCOMPLEX
//#define MKL_INT size_t
//#include "mkl.h" //MKL_Complex16
//#include "mkl_types.h"
