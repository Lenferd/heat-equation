//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <omp.h>
#include <cstdio>
#include <immintrin.h>
#include "Task.h"
#include "StructDeclamer.h"

const int ENABLE_PARALLEL = 1;

struct MatrixValue;

struct SparseMatrix {
    int _size;
    int _rows;
    double *values;
    int *columns;   // какой столбец
    int *pointerB;  // указатель на начало строки

    int threads;
};


void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2);

void fillMatrix3d6Expr(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);
void fillMatrix3d6Expr_wo_boundaries(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);
void fillMatrix3d6Expr_wo_boundaries_for_xyz(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ);

void boundaries_matrix_fix(double *&vect, int sizeX, int sizeY, int sizeZ);
void boundaries_matrix_fix_for_xyz(double *&vect, int sizeX, int sizeY, int sizeZ);

void multiplicateVectorAVXLine(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVectorAVXColumn(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVectorAVXColumn2(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVectorAVXColumn3(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVectorAVXColumn4(SparseMatrix &sp, double *&vect, double *&result, int size, int sizeX, int sizeY, int sizeZ);
void multiplicateVectorAVXColumn5(SparseMatrix &sp, double *&vect, double *&result, int size);

void multiplicateVectorAVXBlocks(SparseMatrix &sp, double *&vect, double *&result, int size);

void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVector_wo_boundaries(SparseMatrix &sp, double *&vect, double *&result, int size);
void multiplicateVectorRunge(SparseMatrix &sp, double *&vect, double *&additional_vect, double *&result, int size);

void spMatrixInit(SparseMatrix &sp, int size, int rows, int threads);
void printVectors(SparseMatrix &sp);


#endif //SPARSEMATRIX_SPARSEMATRIX_H
