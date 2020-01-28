//
//  pruebas.c
//  Ising_mod2
//
//  Created by Pablo BalRes on 20/01/2020.
//  Copyright © 2020 Pablo BalRes. All rights reserved.
//

#include "pruebas.h"
#include <stdio.h>
#include <stdlib.h>  // malloc function
#include <stdbool.h> //library to include boolean algebra
//#include <complex.h> //include complex numbers
#include <math.h>
#include <complex.h> //include complex numbers
#include "mkl.h"

bool BoundCond;

bool DiagType = false; //In case we want to implement other packs of diagonalization

//static variables useful to fix the parameters of the problem
static int NumTotSpin; //total number of spin of the system
//static int NumTotStates; // total number of states 2^NumTotSpin, ==> the dimention of the hilber space ( =dim(H) )
//int **Spin=NULL;
static double hfield; //value of the longitudinal magnetic field h
static double gfield; //value of the trasversal magnetic field g
//static double hfield_copy; //copy of h field, hfieldcopy used only in DAVIDSON algort (see Diag_Hamiltonian() ).



void Initialize_Data()
{
    FILE *Condition_file;
    int temp1;
    MKL_INT a,b;
    int ii=4;
    a=1;
    b=a;
    //MKL_Complex16 z;
    MKL_Complex16 z[ii][ii];
    
    for (int i=0; i<ii; i++) {
        for (int j=0; j<ii; j++) {
            z[i][j].real=1;
            //z[i][j].imag=0.1;
            printf("(%lf+%lf)",z[i][j].real,z[i][j].imag);
        }
        printf("\n");
    }
    z[0][0].real=1;
    
    /*
    z.real=1;
    z.imag=8;
    z.real=z.real+1+ii+b;
    printf("LOOOOK    a+b=%d or %lf + i %lf\n\n\n", a + b, z.real,z.imag);
    */
    
    //open the file to read the conditions
    Condition_file = fopen("Condition.in", "r");
    // check if there is a file
    if (Condition_file == NULL)
    {

        perror("Error while opening the file:\n");
        printf("%s [ERROR]\n","Condition.in");
        exit(EXIT_FAILURE);
    }
    if(!feof(Condition_file)) {
        fscanf(Condition_file, "%d[^\n]", &NumTotSpin);
        fscanf(Condition_file, "%*[^\n]");
        
        fscanf(Condition_file, "%lf[^\n]", &gfield);
        fscanf(Condition_file, "%*[^\n]");

        fscanf(Condition_file, "%lf[^\n]", &hfield);
        fscanf(Condition_file, "%*[^\n]");

        fscanf(Condition_file, "%d[^\n]",&temp1);
        BoundCond = temp1;
        fscanf(Condition_file, "%*[^\n]");

        fscanf(Condition_file, "%d[^\n]", &temp1);
        DiagType = temp1 ;
        printf("\nInside Condition.in file,\n");
    }
        
    

        
        
    
    fclose(Condition_file);
}


int main(){
    BoundCond = false;
    
    Initialize_Data();
    
    //printf("que picha %d\n",BoundCond);
    /*
    MKL_Complex16 c,a[2],b[2];
    
    int n =2;
    a[0].real = 0;
    b[0].real = 0;
    a[1].real = 0;
    b[1].real = 0;
    
    a[0].imag = 0;
    b[0].imag = 0;
    a[1].imag = 0;
    b[1].imag = 0;
    c.real= 5;
    
    c.real = c.real * c.real;
    
    MKL_Complex16 y;
    vzPow( n, a, b, &y );
    
    printf("\n %lf + %lf i \n",y.real, y.imag);
    */
    
    
    
    
    int num, *numptr ;
    num =1;
    numptr = &num;
    int sum = *numptr + num;
    printf("\n %d + %d = %d \n", num,*numptr, sum);
    int NumTotSpin = 3;
    int NumTotStates = pow(2,NumTotSpin);
    
    int **Spin;
    Spin = (int **) malloc(NumTotStates * sizeof(int*));
    
    for (int i=0; i < NumTotStates; i++) {
        Spin[i] = (int *) malloc(NumTotSpin * sizeof(int));
    }
    
    int iTemporal; //variable to allocate the last value position, for then use it as a random generator
    
    for (int i=0; i < NumTotStates; i++) {
        for (int j=0; j < NumTotSpin; j++) {
            Spin[i][j]=0;
            //printf("Spin[%d][%d] = %d \n", i, j, Spin[i][j]);
        }
        
    }
    printf("\n Binari states: \n");
    

    for (int i=0; i < NumTotStates; i++) {
        iTemporal = i;
        printf("\n");
        for (int j=0; j < NumTotSpin; j++) {
            Spin[i][j] = iTemporal % 2;
            
            printf("\tSpin[%d][%d]= %d \n", i, j, Spin[i][j]);
            iTemporal /= 2;
        }
    }
    
    
    double complex z1 = 1;
    double complex z2 = 1*I;
    
    printf("\n\n\n complex prueba %lf %lf \n\n\n",creal(z1), cimag(z2) );
    
    
    
    
    
    printf("Parameters \n\t# Spins %d \n\tg value %lf \n\th value %lf \n\tBoundary condition %d \n\tDiagonalization method %d \n",NumTotSpin,gfield,hfield,BoundCond,DiagType);
    printf("%d\n",NumTotSpin*NumTotSpin);
    
    MKL_Complex16 **Ham;
    Ham = (MKL_Complex16 **) malloc(NumTotStates * sizeof(MKL_Complex16*));
    
    for (int i=0; i < NumTotStates; i++) {
        Ham[i] = (MKL_Complex16 *) malloc(NumTotSpin * sizeof(MKL_Complex16));
    }
    
    return 0;
}
