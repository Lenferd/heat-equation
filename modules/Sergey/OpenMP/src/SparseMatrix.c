//
// Created by Lenferd on 17/10/12.
//
#include "SparseMatrix.h"

void spMatrixInit(SparseMatrix *sp, int size, int rows, int threads) {
    sp->_size = size;
    sp->_rows = rows;
    sp->values = (double *) malloc(size * sizeof(double));
    sp->columns = (int *) malloc(size * sizeof(int));
    sp->pointerB = (int *) malloc((rows + 1) * sizeof(int));
    sp->threads = threads;
}

void fillMatrix3d6Expr_wo_boundaries_for_xyz(SparseMatrix *sp, MatrixVal *matrVal, Task *task) {
    int fixedSizeX = task->nX + 2;
    int fixedSizeY = task->nY + 2;
    int fixedSizeZ = task->nZ + 2;

// size without additional boundaries xz xy)
    int realSizeY = task->nX + 2;
    int realSizeZ = realSizeY * fixedSizeY;

// offset for section start
    int offsetSizeZ = realSizeY * fixedSizeY;

    int index = 0;
    int pIndex = 0;


    int sectionStart = 0;
    for (int z = 0; z < fixedSizeZ; ++z) {
        for (int y = 0; y < fixedSizeY; ++y) {
            sectionStart = z * offsetSizeZ + y * realSizeY;

            int fixBounds = 0;

            for (int x = 0; x < fixedSizeX; ++x) {
                if (x == 0) {
                    fixBounds = 1;
                } else if ((x + 1) == fixedSizeX) {
                    fixBounds = -1;
                }
                // Z first
                sp->values[index] = matrVal->z1;
                sp->columns[index] = z == 0 ?
                                    fixBounds + x + sectionStart + realSizeZ * (fixedSizeZ - 1) :
                                    fixBounds + x + sectionStart - realSizeZ;
                sp->pointerB[pIndex++] = index;
                ++index;


                // Y first
                sp->values[index] = matrVal->y1;
                sp->columns[index] = y == 0 ?
                                    fixBounds + x + sectionStart + realSizeY * (fixedSizeY - 1) :
                                    fixBounds + x + sectionStart - realSizeY;
                ++index;

                // X Group center
                sp->values[index] = matrVal->x1;
                sp->columns[index] = sectionStart + fixBounds + x - 1;
                ++index;

                sp->values[index] = matrVal->x2Comp;
                sp->columns[index] = sectionStart + fixBounds + x;
                ++index;

                sp->values[index] = matrVal->x1;
                sp->columns[index] = sectionStart + fixBounds + x + 1;
                ++index;

                // Y second
                sp->values[index] = matrVal->y1;
                sp->columns[index] = y == fixedSizeY - 1 ?
                                    fixBounds + x + sectionStart - realSizeY * (fixedSizeY - 1) :
                                    fixBounds + x + sectionStart + realSizeY;
                ++index;

                // Z second
                sp->values[index] = matrVal->z1;
                sp->columns[index] = z == fixedSizeZ - 1 ?
                                    fixBounds + x + sectionStart - realSizeZ * (fixedSizeZ - 1) :
                                    fixBounds + x + sectionStart + realSizeZ;
                ++index;

                // afterloop bound fix value clearing
                fixBounds = 0;
            }
        }
    }

    sp->pointerB[pIndex] = index + 1;   //end
}

int boundariesFix_forAdditionalXYZ(double *vect, Task *task) {

    int realSizeX = task->nX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * (task->nY + 2);

    int fixedSizeX = task->nX + 2;
    int fixedSizeY = task->nY + 2;
    int fixedSizeZ = task->nZ + 2;

    int sectionStart;
    int sectionEnd;

    // z == 0 && z == sizeZ - 1
    for (int y = 0; y < fixedSizeY; ++y) {
        sectionStart = y * realSizeY;
        sectionEnd = (fixedSizeZ - 1) * realSizeZ + y * realSizeY;
        for (int x = 0; x < fixedSizeX; ++x) {
            vect[sectionStart + x] = vect[sectionEnd - realSizeZ + x];
            vect[sectionEnd + x] = vect[sectionStart + realSizeZ + x];
        }
    }

    // y == 0 && y == sizeY - 1
    for (int z = 0; z < fixedSizeZ; ++z) {
        sectionStart = z * realSizeZ;
        sectionEnd = (fixedSizeY - 1) * realSizeY + z * realSizeZ;
        for (int x = 0; x < realSizeX; ++x) {
            vect[sectionStart + x] = vect[sectionEnd - realSizeY + x];
            vect[sectionEnd + x] = vect[sectionStart + realSizeY + x];
        }
    }

    // x == 0 && x == sizeX -1
    for (int p = 0; p < realSizeZ * fixedSizeZ; p += realSizeY) {
        vect[p] = vect[p + 1];
        vect[p + realSizeY - 1] = vect[p + realSizeY - 2];
    }

    return 0;
}

void multiplicateVectorAVXColumn5(SparseMatrix *sp, double *vect, double *result, int size) {

    int *point_b;
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

    for (int i = 0; i < size; i += 4) {
        point_b = (sp->pointerB + i);
        point_b_first = *point_b;

// first, get vect values
        __m256d z1f = _mm256_loadu_pd(vect + sp->columns[point_b_first]);
        __m256d y1f = _mm256_loadu_pd(vect + sp->columns[point_b_first + 1]);

        __m256d x1f = _mm256_loadu_pd(vect + sp->columns[point_b_first + 2]);
        __m256d x2Comp = _mm256_loadu_pd(vect + sp->columns[point_b_first + 3]);
        __m256d x1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 4]);

        __m256d y1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 5]);
        __m256d z1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 6]);


// 8 __m256d

// second, get sp->values
        first_values = (sp->values + point_b_first);
        second_values = (sp->values + point_b_first + 7);
        third_values = (sp->values + point_b_first + 14);
        fourth_values = (sp->values + point_b_first + 21);

// and fit it to __m256d
        __m256d z1fval = _mm256_set_pd(first_values[0], second_values[0], third_values[0], fourth_values[0]);
        __m256d y1fval = _mm256_set_pd(first_values[1], second_values[1], third_values[1], fourth_values[1]);

        __m256d x1fval = _mm256_set_pd(first_values[2], second_values[2], third_values[2], fourth_values[2]);
        __m256d x2Compval = _mm256_set_pd(first_values[3], second_values[3], third_values[3], fourth_values[3]);
        __m256d x1sval = _mm256_set_pd(first_values[4], second_values[4], third_values[4], fourth_values[4]);

        __m256d y1sval = _mm256_set_pd(first_values[5], second_values[5], third_values[5], fourth_values[5]);
        __m256d z1sval = _mm256_set_pd(first_values[6], second_values[6], third_values[6], fourth_values[6]);

// 7 __m256d

        __m256d result1 = _mm256_mul_pd(z1f, z1fval);
        __m256d result2 = _mm256_mul_pd(y1f, y1fval);
        __m256d res_sum1 = _mm256_add_pd(result1, result2);

        result1 = _mm256_mul_pd(x1f, x1fval);
        __m256d res_sum2 = _mm256_add_pd(res_sum1, result1);

        result1 = _mm256_mul_pd(x2Comp, x2Compval);
        res_sum1 = _mm256_add_pd(res_sum2, result1);

        result1 = _mm256_mul_pd(x1s, x1sval);
        res_sum2 = _mm256_add_pd(res_sum1, result1);

        result1 = _mm256_mul_pd(y1s, y1sval);
        res_sum1 = _mm256_add_pd(res_sum2, result1);

        result1 = _mm256_mul_pd(z1s, z1sval);
        res_sum2 = _mm256_add_pd(res_sum1, result1);

        _mm256_storeu_pd(result + i, res_sum2);
    }
}
