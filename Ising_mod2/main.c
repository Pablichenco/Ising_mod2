//
//  main.c
//  Ising_mod2
//
//  Created by Pablo BalRes on 05/11/2019.
//  Copyright © 2019 Pablo BalRes. All rights reserved.
//



#include <stdio.h>
#include <stdlib.h>  // malloc function
#include <stdbool.h> //library to include boolean algebra
#include <complex.h> //include complex numbers
#include "gnuplot_i.h" //library for using gnuplot in C by N. Devillard
#include "Math_MKL_Complex16.h" //contain the MKL_Complex16 sum multiplication and LAPACKE_zheev rutine


/**************************************************
 *                                                *
 *                   Global Variables             *
 *                                                *
 **************************************************/

/*
 Boundary condition TRUE->PBC FALSE->OBC
 PBC: Periodic Boundary Condition
 OBC: Open Boundary Condition
 */
bool BoundCond;

bool DiagType = false; //In case we want to implement other packs of diagonalization

//static variables useful to fix the parameters of the problem
static int TotSpin; //total number of spin of the system
static int NumTotStates; // total number of states 2^TotSpin, ==> the dimention of the hilber space ( =dim(H) )
int **Spin;
static double hfield; //value of the longitudinal magnetic field h
static double gfield; //value of the trasversal magnetic field g
//static double hfield_copy; //copy of h field, hfieldcopy used only in DAVIDSON algort (see Diag_Hamiltonian() ).

/////////////////////////////////////////////////////////

/**************************************************
 *                                                *
 *       Implemented Functions (PROTOTYPE)        *
 *                                                *
 **************************************************/

//Reading the parameters of the system in a file
void Initialize_Data(void);

//generate binary values for Spin variable
void BinarySpin(void);

//Construc the Hamiltonian of the problem
void Build_Ham(double hh, double gg, MKL_Complex16 **Ham);

//Desition for chose the algorithm Davidson or LAPACK calling Build_Ham and Diag_Matrix (Note: Davidson doesn't implemented for the moment)
void Diag_Hamiltonian(bool DiagType,double *E, MKL_Complex16 **Ham,MKL_Complex16 *Ground_wave);

//Calculate the magnetization of the system
void Magnet(MKL_Complex16 *Ground_wave, double *Mag_x, double *Mag_z);

/////////////////////////////////////////////////////////


int main(int argc, const char * argv[]) {
    
    //local variables
    
    double *E = NULL,*Mag_z,*Mag_x; //Energy, values of Sigma in z and x
    MKL_Complex16 *Ground_wave; //Ground state wave function
    MKL_Complex16 **Ham; //Hamiltonian of the problem
    FILE *SystemOutput;
    gnuplot_ctrl    * GNU_Frame;
    
/////////////////////////////////////////////////////////

    //inicialice the data of the problem, using Initialize_Data()
    Initialize_Data();
    char *BC,*DT;
    if (BoundCond) {
        BC="PBC";
        if (DiagType) {
            DT="Davidson";
        }else DT="Lapack";
    }
    else{
        BC="OBC";
        if (DiagType) {
            DT="Davidson";
        }else DT="Lapack";
    }
    printf("Parameters \n\t# Spins %d \n\tg value %lf \n\th value %lf \n\tBoundary condition %s \n\tDiagonalization method %s \n",TotSpin,gfield,hfield,BC,DT);
    
/////////////////////////////////////////////////////////

    // Hilbert space dimension dim( H ) = 2^#total of spins
    NumTotStates = pow(2,TotSpin);
    
/////////////////////////////////////////////////////////

    //reserve space for the Spin variables
    Spin = (int **) malloc(NumTotStates * sizeof(int*));
    
    for (int i=0; i < NumTotStates; i++) {
        Spin[i] = (int *) malloc(TotSpin * sizeof(int));
    }
 
    //generate binary values for spin using BinarySpin()
    BinarySpin();

/////////////////////////////////////////////////////////
    //Using calloc instead of malloc, because inizialize the values to cero
    Mag_x = (double *) calloc(TotSpin, sizeof(double));
    Mag_z = (double *) calloc(TotSpin, sizeof(double));
    Ground_wave = (MKL_Complex16 *) calloc(NumTotStates, sizeof(MKL_Complex16));
    E = (double *) calloc(NumTotStates, sizeof(double));
    
    Ham = (MKL_Complex16 **) calloc(NumTotStates, sizeof(MKL_Complex16*));

    for(int i=0; i < NumTotStates; i++){
        Ham[i] = (MKL_Complex16 *) calloc(NumTotStates, sizeof(MKL_Complex16));
    }
/////////////////////////////////////////////////////////

    Diag_Hamiltonian( DiagType, E, Ham, Ground_wave);
    
/////////////////////////////////////////////////////////
    
    SystemOutput = fopen("System_out.dat", "w");
    
    if (SystemOutput == NULL)
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    
/////////////////////////////////////////////////////////

    Magnet(Ground_wave,Mag_x,Mag_z);
    //Implemented GNUPLOT to plot data
    GNU_Frame = gnuplot_init();
    
    //FALTAAAAAAAAA salvar los datos obtenidos

    printf("\n ground wave after Diag_Matrix\n");
    for (int i=0; i < NumTotStates; i++) {
        Ground_wave[i] = Ham[i][0];
        printf("%lf , %lf\n",Ground_wave[i].real,Ground_wave[i].imag);
    }
    
    
    gnuplot_cmd(GNU_Frame, "plot cos(x)");
    
    
    
    fclose(SystemOutput);
    gnuplot_close(GNU_Frame);
    
    /*
     ! Prints out the lowest three energy levels
     PRINT *,'First three energy levels:   ', Energy(1:3)
     
     ! Prints out the magnetizations along X and Z
     IF (PBC) THEN
     PRINT *,'GS magnetization along X:   ', MagnetX(1)
     PRINT *,'GS magnetization along Z:   ', MagnetZ(1)
     ELSE
     WRITE (6,111)  MagnetX
     WRITE (6,112)  MagnetZ
     111  FORMAT ('GS magnetization profile along X:   ', 100(ES15.8,3x))
     112  FORMAT ('GS magnetization profile along Z:   ', 100(ES15.8,3x))
     ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     //CALL cpu_time(cputime)
     //PRINT *,'CPU time for the full process:',CPUtime
     
     STOP
     END PROGRAM ISING
     */
    return 0;
}

void Initialize_Data()
{
    FILE *Condition_file;
    int temp1;
    
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
        fscanf(Condition_file, "%d[^\n]", &TotSpin);
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

void BinarySpin(){
    int iTemporal; //variable to allocate the last value position, for then use it as a random generator
    for (int i=0; i < NumTotStates; i++) {
        iTemporal = i;
        for (int j=0; j < TotSpin; j++) {
            Spin[i][j] = iTemporal % 2;
            iTemporal /= 2;
        }
    }
}

void Diag_Hamiltonian(bool DiagType,double *E, MKL_Complex16 **Ham,MKL_Complex16 *Ground_wave){// Build the hamiltonian and also the diagonalization
    //int kMax = 0; //convenient to manage the use of Davidson algorithms
    
        if (DiagType == true) { //Using Davidson
            
            printf("Not developet jet");
            
        }else { //Using LAPACK
            
            Build_Ham(hfield,gfield,Ham);
           
            Diag_Matrix(Ham,NumTotStates,E);

        }
}

void Build_Ham(double hh, double gg, MKL_Complex16 **Ham){
    int Exchange = 0;
    
/********************************************
*                                           *
*   sigma_x = |0  1|     sigma_z = |1  0|   *
*             |1  0|               |0 -1|   *
*                                           *
*********************************************/
    
    //Sigma_z Sigma_z [coupling]
    for (int j = 0; j < TotSpin - 1; j++) {
        for (int i = 0; i < NumTotStates; i++) {
            //printf("\nspin[%d][%d]",i,j);
            if (Spin[i][j] == Spin[i][j+1]) Ham[i][i].real += -1.0;
            if (Spin[i][j] != Spin[i][j+1]) Ham[i][i].real += 1.0;
        }
    }
    // check periodic boundary conditions
    if (BoundCond == true) {//PBC
        for (int i = 0; i < NumTotStates; i++) {
            if (Spin[i][TotSpin-1] == Spin[i][0]) Ham[i][i].real += -1.0;
            if (Spin[i][TotSpin-1] != Spin[i][0]) Ham[i][i].real += -1.0;
        }
    }
    //Sigma_z  [longitudinal field]
    for (int j = 0; j < TotSpin; j++) {
        for (int i = 0; i < NumTotStates; i++) {
            if (Spin[i][j] == 1) Ham[i][i].real += - hh;
            if (Spin[i][j] == 0) Ham[i][i].real += + hh;
        }
    }
    //Sigma_x  [transverse field]
    for (int j = 0; j < TotSpin; j++) {
        for (int i = 0; i < NumTotStates; i++) {
            if (Spin[i][j] == 1){ Exchange = i - pow(2, j-1);}
            if (Spin[i][j] == 0){ Exchange = i + pow(2, j-1);}
            Ham[Exchange][i].real += + gg;
        }
    }
}

void Magnet(MKL_Complex16 *Ground_wave, double *Mag_x, double *Mag_z){
    MKL_Complex16 Mx_Sum;
    
    Mx_Sum = multipl_conj_MKL(Conju_MKL(Ground_wave[0]), Ground_wave[0]);
    double MagZ = 0.0;
    int Exchange = 0;

    for (int j = 0; j < TotSpin; j++) {
        Mx_Sum.real = 0.0;
        Mx_Sum.imag = 0.0;
        for (int i = 0; i < NumTotStates; i++) {
            if (Spin[i][j] == 1) Mag_z[j] += multipl_conj_MKL(Ground_wave[i], Conju_MKL(Ground_wave[i] ) ).real; /*?????sqrt(multipl_conj_MKL(Ground_wave[i], Conju_MKL(Ground_wave[i] ) ).real) ;  //fabs( pow( Ground_wave[i], 2.0 ) );*/
            
            if (Spin[i][j] == 0) Mag_z[j] += - multipl_conj_MKL(Ground_wave[i], Conju_MKL(Ground_wave[i] ) ).real;/*?????sqrt(multipl_conj_MKL(Ground_wave[i], Conju_MKL(Ground_wave[i] ) ).real);*/
            
            if (Spin[i][j] == 1) Exchange = i - pow( 2.0, j - 1 );
            
            if (Spin[i][j] == 0) Exchange = i + pow( 2.0, j - 1 );
            
            Mx_Sum = Sum_MKL(Mx_Sum, multipl_conj_MKL(Conju_MKL(Ground_wave[0]), Ground_wave[0]));// Mx_Sum +=....
        }
        if ( fabs(Mx_Sum.imag) > 1.e-10 ){
            printf("Unreal magnetization");
            break;
        }
        Mag_x[j] = Mx_Sum.real;
    }
    Exchange = 0;
    for (int i = 0; i < NumTotStates; i++) {
        Exchange = 0;
        for (int j = 0; j < TotSpin; j++) {
            if (Spin[i][j] == 1) Exchange += 1;
            if (Spin[i][j] == 0) Exchange += - 1;
        }
        MagZ += fabs( 1.0 * Exchange ) * multipl_conj_MKL( Conju_MKL( Ground_wave[i] ), Ground_wave[i] ).real;//pow( cabs( Ground_wave[i] ), 2.0 );
    }
    printf("\nSymmetry-broken magnetization along Z: %lf \n", MagZ / TotSpin);
}
