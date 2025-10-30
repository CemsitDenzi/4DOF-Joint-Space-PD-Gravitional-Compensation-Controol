#include "MatriceOps.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


Matrice* clearMat(Matrice* mat)
{
    for (size_t i = 0; i < mat->row; i++)
    {
        for (size_t j = 0; j < mat->col; j++)
        {
            mat->arr[i * mat->col + j]= 0.0 ; 
        }
        
    }
    return mat; 
}



Matrice* matrixSum(Matrice* out_mat,const Matrice* mat1, const Matrice* mat2 )
{
    if(mat1 ||mat2||out_mat) return NULL;

    for (size_t i = 0; i < out_mat->row; i++)
    {
        for (size_t j = 0; j < out_mat->col; j++)
        {
            size_t idx = i * out_mat->col + j ;
            out_mat->arr[idx]=out_mat->arr[idx] + out_mat->arr[idx] ; 
        }
        
    }
    return out_mat; 
}

Matrice* matrixMultiply(Matrice* out_mat,const Matrice* mat1, const Matrice* mat2 )
{
    if(!mat1 || !mat2|| !out_mat) return NULL;
    if(mat1->col != mat2->row ) return NULL;
    if(out_mat->row != mat1->row || out_mat->col != mat2->col) return NULL; // 


    for (size_t i = 0; i < out_mat->row; i++)
    {
        for (size_t j = 0; j < out_mat->col; j++)
        {
          double sum = 0.0 ; 
          for (size_t k = 0; k < mat1->col; k++)
          {
            sum += mat1->arr[i * mat1->col + k] * mat2->arr[k * mat2->col + j];
          }
          out_mat->arr[i * out_mat->col + j] = sum;
        }
    }
    return out_mat; 
}

Matrice* transpose(Matrice* out_mat, const Matrice* mat)
{
    for (size_t i = 0; i < mat->row; i++)
    {
        for (size_t j = 0; j < mat->col; j++)
        {
            out_mat->arr[j * out_mat->col + i] = mat->arr[i * mat->col + j]; 
        } 
    }
    
    return(out_mat);
}

void printMat(const Matrice* mat)
{
    for (size_t i = 0; i < mat->row; i++)
    {
        for (size_t j = 0; j < mat->col; j++)
        {
            printf("%.3f  ", mat->arr[i * mat->col + j]);
        }
        printf("\n");
        
    }
    
}


Matrice* elwiseMultiply(Matrice* out_mat,const Matrice* mat, int scalar)
{
    if(out_mat->col != mat->col || out_mat->row != mat->row)
        return NULL; 

    for (size_t i = 0; i < mat->row; i++)
    {
        for (size_t j = 0; j < mat->col; j++)
        {
            MAT(out_mat,i,j) = scalar * MAT(mat,i,j);
        }
        
    }
    return out_mat; 
}


//Note input parameters row1 and row2 indexing starts with 0;
static Matrice* swapRows(Matrice* mat,int row1 , int row2)
{
    if(row1 >= mat->row || row2 >= mat->row){
        fprintf(stderr, "\nIndex exceeds max row size.");
        exit(1);
    }

    for (size_t i = 0; i < mat->col; i++) {
        double tmp = MAT(mat,row1,i);
        MAT(mat,row1,i) = MAT(mat,row2,i);
        MAT(mat,row2,i) = tmp;
    }
    return mat; 
}

//Returns index of pivot row(Duzenlencek)
static int getPivotRow(const Matrice* mat,int col){
    size_t n = mat->row < mat->col ? mat->row : mat->col ; 

    double max_val = fabs(MAT(mat,col,col));
    int pivot = col;

    for (size_t j = col+1; j < n; j++)
    {
        double val =fabs(MAT(mat,j,col)) ;

        if(val > max_val)
            {
            max_val = val;
            pivot = j ;
        }

    if(max_val < 1e-12)
        return -1; 
    
}
return pivot;
}

// Need to free memory after usage
Matrice* inverseMat(Matrice* out_mat,const Matrice* mat)
{
    if(mat->col != mat->row) return NULL;
    size_t n = mat->col; 
    Matrice aug_mat = {n,2 * n,(double*) malloc(sizeof(double) * n * n * 2) } ;

    for (size_t i = 0; i < n; i++)
    {

        for (size_t j = 0; j < n; j++)
        {
            MAT(&aug_mat,i,j) = MAT(mat,i,j);
        }

        for (size_t j = n; j < 2*n; j++)
        {
            MAT(&aug_mat,i,j) = i == (j -n ) ? 1.0 : 0.0 ;
        }   
    }

    for (size_t col = 0; col < n; col++)
    {
        int pivot = getPivotRow(mat,col);

        if (fabs(MAT(&aug_mat, pivot, col)) < 1e-12) {
            printf("Matrice singular.\n");
            free(aug_mat.arr);
            exit(1);  //degiscek
        }

        if(pivot != col)
            swapRows(&aug_mat,pivot,col);

        
        double pivot_val = MAT(&aug_mat, col, col);

        //Normalization
        for (int j = 0; j < 2 * n; j++)
            MAT(&aug_mat, col, j) /= pivot_val;

        for (int i = 0; i < n; i++) {
            if (i == col) continue;
            double factor = MAT(&aug_mat, i, col);
            for (int j = 0; j < 2 * n; j++)
                MAT(&aug_mat, i, j) -= factor * MAT(&aug_mat, col, j);
        }
    }
      
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            MAT(out_mat, i, j) = MAT(&aug_mat, i, j + n);
        

    out_mat->col = mat->col;
    out_mat->row = mat->row;

    free(aug_mat.arr);
    aug_mat.arr = NULL;
    return out_mat;    
}

int main()
{
    double arr1[2][3] = {{1,2,3},{4,5,6}};
    double arr2[4][4] = {{4,2,1,2},{7,8,6,4},{5,3,6,2},{1,4,2,3}};
    double out[4][4] = {0};
    Matrice mat1 = {2,3,(double*)arr1};
    Matrice mat2 = {4,4,(double*)arr2};
    Matrice out_mat = {4,4,(double*) out};
    printMat(elwiseMultiply(&out_mat,&mat2,2));
    return 0;
}