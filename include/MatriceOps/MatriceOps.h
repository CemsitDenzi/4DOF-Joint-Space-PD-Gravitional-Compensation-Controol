#pragma once

#include<stdlib.h>

#define     MAT(mat,i,j)        ((mat)->arr[(i)* (mat)->col + (j)])


typedef struct  Matrice_t 
{
    size_t row;
    size_t col;
    double* arr;

}Matrice;


Matrice* clearMat(Matrice* mat);

Matrice* matrixSum(Matrice* out_mat,const Matrice* mat1, const Matrice* mat2 );

void printMat(const Matrice* mat);

Matrice* matrixMultiply(Matrice* out_mat,const Matrice* mat1, const Matrice* mat2 );

Matrice* transpose(Matrice* out_mat, const Matrice* mat);

Matrice* elwiseMultiply(Matrice* out_mat,const Matrice* mat, int scalar);

Matrice* inverseMat(Matrice* out_mat,const Matrice* mat);