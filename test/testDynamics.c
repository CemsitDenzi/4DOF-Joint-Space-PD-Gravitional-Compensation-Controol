#include "scaraDyn/scaraDynamics.h"
#include<stdlib.h>
#include<stdio.h>

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