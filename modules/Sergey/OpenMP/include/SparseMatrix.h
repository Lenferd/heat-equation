//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <stdlib.h>
#include "Task.h"
#include <immintrin.h>


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

void spMatrixInit(SparseMatrix *sp, int size, int rows, int threads);
int boundariesFix_forAdditionalXYZ(double *vect, Task *task);
void fillMatrix3d6Expr_wo_boundaries_for_xyz(SparseMatrix *sp, MatrixVal* matrVal, Task *task);
void multiplicateVectorAVXColumn5(SparseMatrix *sp, double *vect, double *result, int size);

#endif //SPARSEMATRIX_SPARSEMATRIX_H
