#include <stdlib.h>
#include<stdio.h>
#include<math.h>
#include "scaraDynamics.h"

Matrice* getInertiaMatrix(const maniParams* params,double* q , Matrice* mat)
{
    if(mat->row != NUM_DOF || mat->col != NUM_DOF) return NULL ; // inertia matrice size must be 4x4
    if(!(mat->arr)) return NULL;
    if(!clearMat(mat)) return NULL;

    const double a1 = params->a1;
    const double a2 = params->a2;
    const double d4 = params->d4;
    const double m1 = params->m1;
    const double m2 = params->m2;
    const double m3 = params->m3;
    const double m4 = params->m4;
    const double Izz = params->Izz;    
    const double sin_t1 = sin(q[0]);
    const double sin_t2 = sin(q[1]);
    const double sin_d3 = sin(q[2]);
    const double sin_t4 = sin(q[3]);
    const double cos_t1 = cos(q[0]);
    const double cos_t2 = cos(q[1]);
    const double cos_d3 = cos(q[2]);
    const double cos_t4 = cos(q[3]);


    MAT(mat,0,0) = Izz/2 + (a1*a1*m1)/2 + (a1*a1*m2)/2 + (a1*a1*m3)/2 + (a2*a2*m2)/2 + (a1*a1*m4)/2 + (a2*a2*m3)/2 + (a2*a2*m4)/2 + a1*a2*m2*cos(q[1]) + a1*a2*m3*cos(q[1]) + a1*a2*m4*cos(q[1]);
    MAT(mat,1,1) = Izz/2 + (a2*a2*m2)/2 + (a2*a2*m3)/2 + (a2*a2*m4)/2 ; 
    MAT(mat,2,2) = m3 / 2 + m4 / 2;
    MAT(mat,3,3) = Izz / 2 ;  
    MAT(mat,1,0) = MAT(mat,0,1) = Izz/2 + (a2*a2*m2)/2 + (a2*a2*m3)/2 + (a2*a2*m4)/2 + (a1*a2*m2*cos(q[1]))/2 + (a1*a2*m3*cos(q[1]))/2 + (a1*a2*m4*cos(q[1]))/2;
    MAT(mat,3,0) = MAT(mat,0,3) = -Izz / 2;
    MAT(mat,3,1) = MAT(mat,1,3) = -Izz / 2; 
    return mat ; 


    
}

Matrice* getGravityMatrix(const maniParams* params, Matrice* mat )
{
    if(mat->row != NUM_DOF || mat->col != 1) return NULL;
    if(!(mat->arr)) return NULL;    
    if(!clearMat(mat)) return NULL;
    MAT(mat,2,0) = -(0.5 * G_ACC * params->m3 + G_ACC * params->m4);
    return mat  ;

}

Matrice* getChristoffelSyms(const maniParams* params, double* q,double* qd,Matrice* mat)
{
    if(mat->row != NUM_DOF || mat->col != NUM_DOF) return NULL;
    if(!(mat->arr)) return NULL;
    if(!clearMat(mat)) return NULL;

    double total_mass = (params->m2 + params->m3 + params->m4);
    double s_t2 = sin(q[1]);
    double a1 = params->a1 ;
    double a2 = params->a2;

    MAT(mat,0,0) = -a1 * a2 * qd[1] * s_t2 * total_mass ;
    MAT(mat,1,0) = a1 * a2 * qd[0] * s_t2 * total_mass ;
    MAT(mat,0,1) = -a1 * a2 * s_t2 * (qd[0] + qd[1]) * total_mass ; 

    return mat ; 

}



int main()
{

    maniParams params = {12, 32, 10, 2, 3, 2, 4, 4};
    double q[] = {0,1,2,3};
    double qd[] = {0,1,2,3};
    double arr[4][4] = {{1,23,4,1},{4,5,2,6},{3,6,1,2},{12,2,1,1}};
    Matrice D  = {.arr = (double *)arr, .col = NUM_DOF ,.row = NUM_DOF};
    getInertiaMatrix(&params,q,&D);

    double arr_g[4][1] = {0};
    Matrice G = {.arr = (double*) arr_g , .col = 1 , .row = NUM_DOF};
    getGravityMatrix(&params, &G); 

    double arr_c[4][4] = {0};
    Matrice C = {.arr = (double *) arr_c , .row= NUM_DOF ,.col = NUM_DOF};
    getChristoffelSyms(&params , q , qd , &C);

    printMat(&D);
    printf("\n");
    printMat(&G); 
    printf("\n");
    printMat(&C);

    return 0;
}