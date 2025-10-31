#pragma once

#include "MatriceOps/MatriceOps.h"

#define     NUM_DOF          4 
#define     G_ACC           9.81  

typedef struct manipulatorParams_t{
// link lenghts in mm (D-H parameters)
double a1,a2,d4;

// link masses in kg
double m1,m2,m3,m4,Izz;

}maniParams;



Matrice* getInertiaMatrix(const maniParams* params,double* q,Matrice* mat); 
Matrice* getGravityMatrix(const maniParams* params, Matrice* mat); 
Matrice* getChristoffelSyms(const maniParams* params, double* q,double* qd,Matrice* mat);


