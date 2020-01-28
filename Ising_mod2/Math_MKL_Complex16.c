//
//  Math_MKL_Complex16.c
//  Ising_mod2
//
//  Created by Pablo BalRes on 26/01/2020.
//  Copyright © 2020 Pablo BalRes. All rights reserved.
//

#include "Math_MKL_Complex16.h"

MKL_Complex16 Conju_MKL(MKL_Complex16 V){
    V.imag = -V.imag;
    return V;
}
MKL_Complex16 multipl_conj_MKL(MKL_Complex16 V, MKL_Complex16 W){
    MKL_Complex16 Z;
    Z.real = V.real * W.real - V.imag * W.imag;
    Z.imag = V.imag * W.real + V.real * W.imag;
    return Z;
}
MKL_Complex16 Sum_MKL(MKL_Complex16 V, MKL_Complex16 W){
    MKL_Complex16 Z;
    Z.real = V.real + W.real;
    Z.imag = V.imag + W.imag;
    return Z;
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda ) {
    MKL_INT i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ )
            printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
        printf( "\n" );
    }
}

/* Auxiliary routine: printing a real matrix */
void  print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
    MKL_INT i, j;
    printf( "\n %s\n", desc );
    for( i = 0; i < m; i++ ) {
        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
        printf( "\n" );
    }
}

void Diag_Matrix(MKL_Complex16 **Ham, int NumTot,double *E){
    //int *Dimention, Lwork, LRWork;
    int Info;
    //MKL_Complex16 *Work;
    //double *Rwork;
    //double *W;
    
    //Dimention = &NumTot;
    
    //double Ener[*Dimention];
    //MKL_Complex16 Matr[*Dimention][*Dimention];
    //W = (double *) malloc(NumTot * sizeof(double));
    //E = (double *) malloc(NumTotStates*sizeof(double));

    /*
    Lwork = 22 * (2 * *Dimention - 1);//define the dimention of the Work array
    LRWork = 3 * *Dimention - 2; //define the dimention of the Rwork array
    if(LRWork - 2 > 1){
        Rwork = (double *) malloc(LRWork * sizeof(double));
        
    }else {
        Rwork = (double *) malloc(1 * sizeof(double));
        
    }
    //MKL_Complex16
    Work = (MKL_Complex16 *) malloc(Lwork * sizeof(MKL_Complex16));
    
    Work = 0;
    Rwork = 0;
    */
    /*
     ZHEEV The routine computes all the eigenvalues and, optionally, the eigenvectors of a square complex Hermitian matrix A, in our case the Hamiltonian H.
     For more info check the page http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html
     
     JOBZ  -> 'N':  Compute eigenvalues only or 'V':  Compute eigenvalues and eigenvectors
     UPLO  -> 'U':  Upper triangle of A is stored or 'L':  Lower triangle of A is stored.
     N     -> The order of the matrix A.  N >= 0.
     A     -> COMPLEX*16 array, dimension (LDA, N)
     On entry, the Hermitian matrix A.  If UPLO = 'U', the
     leading N-by-N upper triangular part of A contains the
     upper triangular part of the matrix A.  If UPLO = 'L',
     the leading N-by-N lower triangular part of A contains
     the lower triangular part of the matrix A.
     On exit, if JOBZ = 'V', then if INFO = 0, A contains the
     orthonormal eigenvectors of the matrix A.
     If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
     or the upper triangle (if UPLO='U') of A, including the
     diagonal, is destroyed.
     LDA   -> The leading dimension of the array A.  LDA >= max(1,N)
     W     -> (output)Array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
     WORK  -> (workspace/output) array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     LWORK -> The length of the array WORK.  LWORK >= max(1,2*N-1).
     For optimal efficiency, LWORK >= (NB+1)*N,
     where NB is the blocksize for ZHETRD returned by ILAENV. (NOTE: if LWORK = -1, see the link)
     RWORK -> DOUBLE PRECISION array, dimension (max(1, 3*N-2))
     INFO  -> (output)INTEGER
     = 0:  successful exit
     < 0:  if INFO = -i, the i-th argument had an illegal value
     > 0:  if INFO = i, the algorithm failed to converge; i
     off-diagonal elements of an intermediate tridiagonal
     form did not converge to zero.
     */
    /*
     char uu='U';
     const char *U = &uu;
     char vv='V';
     const char *V = &vv;
     LAPACKE_zheev(V, U, Dimention, *Ham, Dimention, E, Work, &Lwork, Rwork, &Info);
     
     //Check if the process was complished succesfully
     if (Info !=0) {
     
     }
     */
    MKL_Complex16 *a = (MKL_Complex16 *) malloc( NumTot * NumTot *sizeof(MKL_Complex16));
    
    for (int i = 0; i < NumTot; i++) {
        for (int j = 0; j < NumTot; j++) {
            a[i*NumTot+j] = Ham[i][j];
        }
    }
    
    
    
    Info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', 'U', NumTot, a, NumTot, E);
    
    if (Info != 0) {
        printf("There are problems in diagonalization: ");
        (Info > 0)?printf("the algorithm failed to converge"):printf("the i-th parameter had an illegal value.");
        exit(-1);
    }
    
    
    for (int i = 0; i < NumTot; i++) {
        for (int j = 0; j < NumTot; j++) {
            Ham[i][j] = a[i*NumTot+j];
        }
    }
    
    /* Print eigenvalues */
    print_rmatrix( "Eigenvalues", 1, NumTot, E, 1 );
    /* Print eigenvectors */
    print_matrix( "Eigenvectors (stored columnwise)", NumTot, NumTot, a, NumTot );
    
    free(a);
     
}

