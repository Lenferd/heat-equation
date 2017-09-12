//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <stdlib.h>
#include "Task.h"


typedef struct MatrixValue {
    double x1;
    double x2Comp;

    double y1;
    double z1;
} MatrixVal;

typedef struct SparseMatrix {
    int _size;
    int _rows;
    double *values;
    int *columns;   // какой столбец
    int *pointerB;  // указатель на начало строки

    int threads;
} SparseMatrix;

int boundaries_conditions(double *vect, Task *task);

void spMatrixInit(SparseMatrix *sp, int size, int rows, int omp_threads);

// without boundaries
int fillMatrix_3d_MPI_woBound(SparseMatrix *sp, MatrixVal *m_value, ProcSettings *proc);

void multiplicateMatrVector(SparseMatrix *sp, double *input, double *result, int size);

//void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2);
//
//void fillMatrix3d6Expr(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);
//void fillMatrix3d6Expr_wo_boundaries(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);
//
//void boundaries_matrix_fix(double *&vect, int sizeX, int sizeY, int sizeZ);
//
//void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size);
//void spMatrixInit(SparseMatrix &sp, int size, int rows, int threads);
//void printVectors(SparseMatrix &sp);


#endif //SPARSEMATRIX_SPARSEMATRIX_H
