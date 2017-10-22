//
// Created by lenferd on 27.10.16.
//

#include "SparseMatrix.h"
#include "../../MPI/include/SparseMatrix.h"

void spMatrixInit(SparseMatrix &sp, int size, int rows, int threads) {
    sp._size = size;
    sp._rows = rows;
    sp.values = new double[size];
    sp.columns = new int[size];
    sp.pointerB = new int[rows + 1];
    sp.threads = threads;
}

void multiplicateVector(SparseMatrix &sp, double *&vect, double *&result, int size) {

    omp_set_num_threads(sp.threads);

#pragma omp parallel for
    for (int i = 0; i < size; i++) {  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp.pointerB[i]; j < sp.pointerB[i + 1]; j++) {
            local_result += sp.values[j] * vect[sp.columns[j]];
        }
        result[i] = local_result;
    }
}

//void multiplicateVector_wo_boundaries(SparseMatrix &sp, double *&vect, double *&result, int sizeX, int sizeY, int sizeZ) {
//
//    int index = 0;
//    int sectionStart = 0;
//    for (int z = 1; z < sizeZ + 1; ++z) {
//        for (int y = 1; y < sizeY + 1; ++y) {
//            sectionStart = z * (sizeX + 2) *(sizeY+2) + y * (sizeX+2);
//            for (int x = 1 ; x < sizeX + 1; ++x) {
//            }
//
//    for (int i = 0; i < size; i++) {  // iteration FOR RESULT VECTOR!!!
//        double local_result = 0;
//        for (int j = sp.pointerB[i]; j < sp.pointerB[i + 1]; j++) {
//            local_result += sp.values[j] * vect[sp.columns[j]];
//        }
//        result[i] = local_result;
//    }
//}

void multiplicateVectorAVXLine(SparseMatrix &sp, double *&vect, double *&result, int size) {
    int first_point_b;
    double *temp_vect = new double[8];
    double *result_vect;

    for (int i = 0; i < size; i++) {  // iteration FOR RESULT VECTOR!!!
        first_point_b = sp.pointerB[i];

        __m256d first_values = _mm256_loadu_pd(sp.values + first_point_b);
        __m256d second_values = _mm256_loadu_pd(sp.values + first_point_b + 4);

        temp_vect[0] = vect[sp.columns[first_point_b]];
        temp_vect[1] = vect[sp.columns[first_point_b + 1]];
        temp_vect[2] = vect[sp.columns[first_point_b + 2]];
        temp_vect[3] = vect[sp.columns[first_point_b + 3]];

        temp_vect[4] = vect[sp.columns[first_point_b + 4]];
        temp_vect[5] = vect[sp.columns[first_point_b + 5]];
        temp_vect[6] = vect[sp.columns[first_point_b + 6]];
        temp_vect[7] = 0;

        __m256d first_vect = _mm256_loadu_pd(temp_vect);
        __m256d second_vect = _mm256_loadu_pd(temp_vect + 4);

        __m256d first_result = _mm256_mul_pd(first_values, first_vect);
        __m256d second_result = _mm256_mul_pd(second_values, second_vect);
        __m256d third_result = _mm256_add_pd(first_result, second_result);

        result_vect = (double *) &third_result;
        result[i] = result_vect[0] + result_vect[1] + result_vect[2] + result_vect[3];
    }
}

void multiplicateVectorAVXColumn(SparseMatrix &sp, double *&vect, double *&result, int size) {
    int *point_b = new int[4];
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

//    __m256d zero = _mm256_set1_pd(0);

    for (int i = 0; i < size; i += 4) {
        point_b = (sp.pointerB + i);
        point_b_first = *point_b;

        // first, get vect values
        __m256d z1f = _mm256_loadu_pd(vect + sp.columns[point_b_first]);
        __m256d y1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 1]);

        __m256d x1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 2]);
        __m256d x2Comp = _mm256_loadu_pd(vect + sp.columns[point_b_first + 3]);
        __m256d x1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 4]);

        __m256d y1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 5]);
        __m256d z1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 6]);


        // 8 __m256d

        // second, get sp.values
        first_values = (sp.values + point_b_first);
        second_values = (sp.values + point_b_first + 7);
        third_values = (sp.values + point_b_first + 14);
        fourth_values = (sp.values + point_b_first + 21);

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

        __m256d result3 = _mm256_mul_pd(x1f, x1fval);
        __m256d res_sum2 = _mm256_add_pd(res_sum1, result3);

        __m256d result4 = _mm256_mul_pd(x2Comp, x2Compval);
        __m256d res_sum3 = _mm256_add_pd(res_sum2, result4);

        __m256d result5 = _mm256_mul_pd(x1s, x1sval);
        __m256d res_sum4 = _mm256_add_pd(res_sum3, result5);

        __m256d result6 = _mm256_mul_pd(y1s, y1sval);
        __m256d res_sum5 = _mm256_add_pd(res_sum4, result6);

        __m256d result7 = _mm256_mul_pd(z1s, z1sval);
        __m256d res_sum6 = _mm256_add_pd(res_sum5, result7);

        _mm256_storeu_pd(result + i, res_sum6);
        // 7.9 vs 5.3, but abs 0.56 and rel 1535.
    }
}

void multiplicateVectorAVXColumn2(SparseMatrix &sp, double *&vect, double *&result, int size) {
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

//    __m256d zero = _mm256_set1_pd(0);

    for (int i = 0; i < size; i += 4) {
        point_b_first = *(sp.pointerB + i);
        first_values = (sp.values + point_b_first);
        second_values = (sp.values + point_b_first + 7);
        third_values = (sp.values + point_b_first + 14);
        fourth_values = (sp.values + point_b_first + 21);

        __m256d z1f = _mm256_loadu_pd(vect + sp.columns[point_b_first]);
        __m256d z1fval = _mm256_set_pd(first_values[0], second_values[0], third_values[0], fourth_values[0]);
        __m256d result_mul1 = _mm256_mul_pd(z1f, z1fval);

        __m256d y1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 1]);
        __m256d y1fval = _mm256_set_pd(first_values[1], second_values[1], third_values[1], fourth_values[1]);
        __m256d result_mul2 = _mm256_mul_pd(y1f, y1fval);
        __m256d res_sum1 = _mm256_add_pd(result_mul1, result_mul2);


        __m256d x1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 2]);
        __m256d x1fval = _mm256_set_pd(first_values[2], second_values[2], third_values[2], fourth_values[2]);
        result_mul1 = _mm256_mul_pd(x1f, x1fval);
        __m256d res_sum2 = _mm256_add_pd(result_mul1, res_sum1);


        __m256d x2Comp = _mm256_loadu_pd(vect + sp.columns[point_b_first + 3]);
        __m256d x2Compval = _mm256_set_pd(first_values[3], second_values[3], third_values[3], fourth_values[3]);
        result_mul1 = _mm256_mul_pd(x2Comp, x2Compval);
        res_sum1 = _mm256_add_pd(result_mul1, res_sum2);

        __m256d x1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 4]);
        __m256d x1sval = _mm256_set_pd(first_values[4], second_values[4], third_values[4], fourth_values[4]);
        result_mul1 = _mm256_mul_pd(x1s, x1sval);
        res_sum2 = _mm256_add_pd(result_mul1, res_sum1);

        __m256d y1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 5]);
        __m256d y1sval = _mm256_set_pd(first_values[5], second_values[5], third_values[5], fourth_values[5]);
        result_mul1 = _mm256_mul_pd(y1s, y1sval);
        res_sum1 = _mm256_add_pd(result_mul1, res_sum2);

        __m256d z1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 6]);
        __m256d z1sval = _mm256_set_pd(first_values[6], second_values[6], third_values[6], fourth_values[6]);
        result_mul1 = _mm256_mul_pd(z1s, z1sval);
        res_sum2 = _mm256_add_pd(result_mul1, res_sum1);

        _mm256_storeu_pd(result + i, res_sum2);
    }
}

void multiplicateVectorAVXColumn3(SparseMatrix &sp, double *&vect, double *&result, int size) {
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

//    __m256d zero = _mm256_set1_pd(0);

    for (int i = 0; i < size; i += 4) {
        point_b_first = *(sp.pointerB + i);
        first_values = (sp.values + point_b_first);
        second_values = (sp.values + point_b_first + 7);
        third_values = (sp.values + point_b_first + 14);
        fourth_values = (sp.values + point_b_first + 21);

        __m256d v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first]);
        __m256d v_sp = _mm256_set_pd(first_values[0], second_values[0], third_values[0], fourth_values[0]);
        __m256d result_mul1 = _mm256_mul_pd(v_vect, v_sp);

        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 1]);
        v_sp = _mm256_set_pd(first_values[1], second_values[1], third_values[1], fourth_values[1]);
        __m256d result_mul2 = _mm256_mul_pd(v_vect, v_sp);
        __m256d res_sum1 = _mm256_add_pd(result_mul1, result_mul2);


        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 2]);
        v_sp = _mm256_set_pd(first_values[2], second_values[2], third_values[2], fourth_values[2]);
        result_mul1 = _mm256_mul_pd(v_vect, v_sp);
        __m256d res_sum2 = _mm256_add_pd(result_mul1, res_sum1);


        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 3]);
        v_sp = _mm256_set_pd(first_values[3], second_values[3], third_values[3], fourth_values[3]);
        result_mul1 = _mm256_mul_pd(v_vect, v_sp);
        res_sum1 = _mm256_add_pd(result_mul1, res_sum2);

        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 4]);
        v_sp = _mm256_set_pd(first_values[4], second_values[4], third_values[4], fourth_values[4]);
        result_mul1 = _mm256_mul_pd(v_vect, v_sp);
        res_sum2 = _mm256_add_pd(result_mul1, res_sum1);

        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 5]);
        v_sp = _mm256_set_pd(first_values[5], second_values[5], third_values[5], fourth_values[5]);
        result_mul1 = _mm256_mul_pd(v_vect, v_sp);
        res_sum1 = _mm256_add_pd(result_mul1, res_sum2);

        v_vect = _mm256_loadu_pd(vect + sp.columns[point_b_first + 6]);
        v_sp = _mm256_set_pd(first_values[6], second_values[6], third_values[6], fourth_values[6]);
        result_mul1 = _mm256_mul_pd(v_vect, v_sp);
        res_sum2 = _mm256_add_pd(result_mul1, res_sum1);

        _mm256_storeu_pd(result + i, res_sum2);
    }
}

void multiplicateVectorPart(SparseMatrix &sp, double *&vect, double *&result, int start, int size) {

    omp_set_num_threads(sp.threads);

    #pragma omp parallel for
    for (int i = start; i < start + size; i++) {  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp.pointerB[i]; j < sp.pointerB[i + 1]; j++) {
            local_result += sp.values[j] * vect[sp.columns[j]];
        }
        result[i] = local_result;
    }
}


void multiplicateVectorAVXColumn4(SparseMatrix &sp, double *&vect, double *&result, int size, int sizeX, int sizeY,
                                  int sizeZ) {
    int realSizeX = sizeX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * sizeY;

    int *point_b = new int[4];
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

    multiplicateVectorPart(sp, vect, result, 0, realSizeZ); // z == 0
    for (int z = 1; z < sizeY - 1; ++z) {
        multiplicateVectorPart(sp, vect, result, z * realSizeZ, realSizeY); // y == 0
        for (int y = 1; y < sizeY - 1; ++y) {
            multiplicateVectorPart(sp, vect, result, z * realSizeZ + y * realSizeY, realSizeX); // x == 0
            for (int x = 1; x < realSizeX - 1; x += 4) {
                int sectionStart = z * realSizeZ + y * realSizeY;
                if (x + 4 >= realSizeX - 1) {
                    multiplicateVectorPart(sp, vect, result, sectionStart + x, realSizeX - x); // x == 0
                }
                point_b = (sp.pointerB + x);
                point_b_first = *point_b;

                // first, get vect values
                __m256d z1f = _mm256_loadu_pd(vect + sp.columns[point_b_first]);
                __m256d y1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 1]);

                __m256d x1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 2]);
                __m256d x2Comp = _mm256_loadu_pd(vect + sp.columns[point_b_first + 3]);
                __m256d x1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 4]);

                __m256d y1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 5]);
                __m256d z1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 6]);


                // 8 __m256d

                // second, get sp.values
                first_values = (sp.values + point_b_first);
                second_values = (sp.values + point_b_first + 7);
                third_values = (sp.values + point_b_first + 14);
                fourth_values = (sp.values + point_b_first + 21);

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

                _mm256_storeu_pd(result + x, res_sum2);
                // 7.9 vs 5.3, but abs 0.56 and rel 1535.
            }
            multiplicateVectorPart(sp, vect, result, (z * realSizeZ) + ((y+1) * realSizeY) - realSizeX, realSizeX); // x == sizeX - 1
        }
        multiplicateVectorPart(sp, vect, result, ((z + 1) * realSizeZ) - realSizeY, realSizeY); // y == sizeY - 1
    }
    multiplicateVectorPart(sp, vect, result, size - realSizeZ, realSizeZ); // z == sizeZ -1
}


void multiplicateVectorAVXColumn5(SparseMatrix &sp, double *&vect, double *&result, int size) {

    int *point_b;
    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

    for (int i = 0; i < size; i += 4) {
        point_b = (sp.pointerB + i);
        point_b_first = *point_b;

        // first, get vect values
        __m256d z1f = _mm256_loadu_pd(vect + sp.columns[point_b_first]);
        __m256d y1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 1]);

        __m256d x1f = _mm256_loadu_pd(vect + sp.columns[point_b_first + 2]);
        __m256d x2Comp = _mm256_loadu_pd(vect + sp.columns[point_b_first + 3]);
        __m256d x1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 4]);

        __m256d y1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 5]);
        __m256d z1s = _mm256_loadu_pd(vect + sp.columns[point_b_first + 6]);


        // 8 __m256d

        // second, get sp.values
        first_values = (sp.values + point_b_first);
        second_values = (sp.values + point_b_first + 7);
        third_values = (sp.values + point_b_first + 14);
        fourth_values = (sp.values + point_b_first + 21);

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

void multiplicateVectorRunge(SparseMatrix &sp, double *&vect, double *&additional_vect, double *&result, int size) {

    omp_set_num_threads(sp.threads);

#pragma omp parallel for
    for (int i = 0; i < size; i++) {  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp.pointerB[i]; j < sp.pointerB[i + 1]; j++) {
            local_result += sp.values[j] * vect[sp.columns[j]];
        }
        result[i] = local_result + additional_vect[i];
    }
}


void fillMatrix2Expr(SparseMatrix &sp, int size, double expr1, double expr2) {
    int index = 0;
    int pIndex = 0;

    sp.values[index] = 1;
    sp.columns[index] = 0;
    sp.pointerB[pIndex++] = 0;
    ++index;

    for (int i = 1; i < size - 1; ++i) {

        //printf("index %d \n", index);
        sp.values[index] = expr1;
        sp.columns[index] = i - 1;
        sp.pointerB[pIndex++] = index;
        ++index;

        sp.values[index] = expr2;
        sp.columns[index] = i;
        ++index;

        sp.values[index] = expr1;
        sp.columns[index] = i + 1;
        ++index;
    }

    sp.values[index] = 1;
    sp.columns[index] = size - 1;
    sp.pointerB[pIndex++] = index;

    sp.pointerB[pIndex] = index + 1;   //end
}

void fillMatrix3d6Expr(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ) {
    int realSizeX = sizeX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * sizeY;
    int index = 0;
    int pIndex = 0;


    int sectionStart = 0;
    for (int z = 0; z < sizeZ; ++z) {
        for (int y = 0; y < sizeY; ++y) {
            sectionStart = z * realSizeZ + y * realSizeY;

            /** Boundaries rule
             *  If we on the edge, we should use same expression (line with parametrs), as the line after.
             *  If it's first line the pattern for her is line two. (+1)
             *  If it's last line, pattern - previous line.         (-1)
             *  Realization - fixes value, whose start to work, if we on the boundaries, joins @var x
             *  @var fixBounds
             */

            int fixBounds = 0;

            for (int x = 0; x < realSizeX; ++x) {
                if (x == 0) {
                    fixBounds = 1;
                } else if ((x + 1) == realSizeX) {
                    fixBounds = -1;
                }
                // Z first
                sp.values[index] = taskexpr.z1;
                sp.columns[index] = z == 0 ?
                                    fixBounds + x + sectionStart + realSizeZ * (sizeZ - 1) :
                                    fixBounds + x + sectionStart - realSizeZ;
                sp.pointerB[pIndex++] = index;
                ++index;


                // Y first
                sp.values[index] = taskexpr.y1;
                sp.columns[index] = y == 0 ?
                                    fixBounds + x + sectionStart + realSizeY * (sizeY - 1) :
                                    fixBounds + x + sectionStart - realSizeY;
                ++index;

                // X Group center
                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + fixBounds + x - 1;
                ++index;

                sp.values[index] = taskexpr.x2Comp;
                sp.columns[index] = sectionStart + fixBounds + x;
                ++index;

                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + fixBounds + x + 1;
                ++index;

                // Y second
                sp.values[index] = taskexpr.y1;
                sp.columns[index] = y == sizeY - 1 ?
                                    fixBounds + x + sectionStart - realSizeY * (sizeY - 1) :
                                    fixBounds + x + sectionStart + realSizeY;
                ++index;

                // Z second
                sp.values[index] = taskexpr.z1;
                sp.columns[index] = z == sizeZ - 1 ?
                                    fixBounds + x + sectionStart - realSizeZ * (sizeZ - 1) :
                                    fixBounds + x + sectionStart + realSizeZ;
                ++index;

                // afterloop bound fix value clearing
                fixBounds = 0;
            }
        }
    }

    sp.pointerB[pIndex] = index + 1;   //end
}

void fillMatrix3d6Expr_wo_boundaries(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ) {
    int realSizeX = sizeX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * sizeY;
    int index = 0;
    int pIndex = 0;


    int sectionStart = 0;
    for (int z = 0; z < sizeZ; ++z) {
        for (int y = 0; y < sizeY; ++y) {
            sectionStart = z * realSizeZ + y * realSizeY;

            for (int x = 0; x < realSizeX; ++x) {
                // Z first
                sp.values[index] = taskexpr.z1;
                sp.columns[index] = z == 0 ?
                                    x + sectionStart + realSizeZ * (sizeZ - 1) :
                                    x + sectionStart - realSizeZ;
                sp.pointerB[pIndex++] = index;
                ++index;


                // Y first
                sp.values[index] = taskexpr.y1;
                sp.columns[index] = y == 0 ?
                                    x + sectionStart + realSizeY * (sizeY - 1) :
                                    x + sectionStart - realSizeY;
                ++index;

                // X Group center
                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + x - 1;
                ++index;

                sp.values[index] = taskexpr.x2Comp;
                sp.columns[index] = sectionStart + x;
                ++index;

                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + x + 1;
                ++index;

                // Y second
                sp.values[index] = taskexpr.y1;
                sp.columns[index] = y == sizeY - 1 ?
                                    x + sectionStart - realSizeY * (sizeY - 1) :
                                    x + sectionStart + realSizeY;
                ++index;

                // Z second
                sp.values[index] = taskexpr.z1;
                sp.columns[index] = z == sizeZ - 1 ?
                                    x + sectionStart - realSizeZ * (sizeZ - 1) :
                                    x + sectionStart + realSizeZ;
                ++index;

            }
        }
    }

    sp.pointerB[pIndex] = index + 1;   //end
}

void fillMatrix3d6Expr_wo_boundaries_for_xyz(SparseMatrix &sp, MatrixValue &taskexpr, int sizeX, int sizeY, int sizeZ) {
    int fixedSizeX = sizeX + 2;
    int fixedSizeY = sizeY + 2;
    int fixedSizeZ = sizeZ + 2;

    // size without additional boundaries xz xy)
    int realSizeY = sizeX+2;
    int realSizeZ = realSizeY * fixedSizeY;

    // offset for section start
    int offsetSizeZ = realSizeY * fixedSizeY;

    int index = 0;
    int pIndex = 0;


    int sectionStart = 0;
    for (int z = 0; z < fixedSizeZ; ++z) {
        for (int y = 0; y < fixedSizeY; ++y) {
            sectionStart = z * offsetSizeZ + y * realSizeY;

            for (int x = 0; x < fixedSizeX; ++x) {
                sp.values[index] = taskexpr.z1;
                sp.columns[index] = x + sectionStart - realSizeZ;
                sp.pointerB[pIndex++] = index;
                ++index;


                // Y first
                sp.values[index] = taskexpr.y1;
                sp.columns[index] =
                                    x + sectionStart - realSizeY;
                ++index;

                // X Group center
                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + x - 1;
                ++index;

                sp.values[index] = taskexpr.x2Comp;
                sp.columns[index] = sectionStart + x;
                ++index;

                sp.values[index] = taskexpr.x1;
                sp.columns[index] = sectionStart + x + 1;
                ++index;

                // Y second
                sp.values[index] = taskexpr.y1;
                sp.columns[index] =
                                    x + sectionStart + realSizeY;
                ++index;

                // Z second
                sp.values[index] = taskexpr.z1;
                sp.columns[index] =
                                    x + sectionStart + realSizeZ;
                ++index;

            }
        }
    }

    sp.pointerB[pIndex] = index + 1;   //end
}

void printVectors(SparseMatrix &sp) {
    printf("values\n");
    for (int i = 0; i < sp._size; ++i) {
        printf("%lf ", sp.values[i]);
    }
    printf("\n");

}

void boundaries_matrix_fix(double *&vect, int sizeX, int sizeY, int sizeZ) {
    int realSizeX = sizeX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * sizeY;


    int sectionStart = 0;
#pragma omp parallel for
    for (int z = 0; z < sizeZ; ++z) {
        for (int y = 0; y < sizeY; ++y) {
            sectionStart = z * realSizeZ + y * realSizeY;

            for (int x = 0; x < realSizeX; ++x) {
                if (x == 0) {
                    vect[sectionStart + x] = vect[sectionStart + x + 1];
                } else if ((x + 1) == realSizeX) {
                    vect[sectionStart + x] = vect[sectionStart + x - 1];
                } else { ;
                }
            }
        }
    }
}

void boundaries_matrix_fix_for_xyz(double *&vect, int sizeX, int sizeY, int sizeZ) {
    int realSizeX = sizeX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * (sizeY + 2);

    int fixedSizeX = sizeX + 2;
    int fixedSizeY = sizeY + 2;
    int fixedSizeZ = sizeZ + 2;

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
    for (int p = 0; p < realSizeZ * fixedSizeZ; p+=realSizeY) {
        vect[p] = vect[p+1];
        vect[p + realSizeY - 1] = vect[p + realSizeY - 2];
    }

}

#define _MM_TRANSPOSE4_PD(row0,row1,row2,row3)                                 \
                {                                                                \
                    __m256d tmp3, tmp2, tmp1, tmp0;                              \
                                                                                 \
                    tmp0 = _mm256_shuffle_pd((row0),(row1), 0x0);                    \
                    tmp2 = _mm256_shuffle_pd((row0),(row1), 0xF);                \
                    tmp1 = _mm256_shuffle_pd((row2),(row3), 0x0);                    \
                    tmp3 = _mm256_shuffle_pd((row2),(row3), 0xF);                \
                                                                                 \
                    (row0) = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);   \
                    (row1) = _mm256_permute2f128_pd(tmp2, tmp3, 0x20);   \
                    (row2) = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);   \
                    (row3) = _mm256_permute2f128_pd(tmp2, tmp3, 0x31);   \
                }

void multiplicateVectorAVXColumn5_shuffle(SparseMatrix *sp, double *vect, double *result, int size) {

    int point_b_first;

    double *first_values;
    double *second_values;
    double *third_values;
    double *fourth_values;

    for (int i = 0; i < size; i += 4) {
        point_b_first = *(sp->pointerB + i);

        // first, get vect values
        __m256d z1f = _mm256_loadu_pd(vect + sp->columns[point_b_first]);
        __m256d y1f = _mm256_loadu_pd(vect + sp->columns[point_b_first + 1]);

        __m256d x1f = _mm256_loadu_pd(vect + sp->columns[point_b_first + 2]);
        __m256d x2Comp = _mm256_loadu_pd(vect + sp->columns[point_b_first + 3]);
        __m256d x1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 4]);

        __m256d y1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 5]);
        __m256d z1s = _mm256_loadu_pd(vect + sp->columns[point_b_first + 6]);


        __m256d z1fval = _mm256_loadu_pd(sp->values + point_b_first);

        __m256d y1fval = _mm256_loadu_pd(sp->values + point_b_first + 7);

        __m256d x1fval = _mm256_loadu_pd(sp->values + point_b_first + 14);

        __m256d x2Compval = _mm256_loadu_pd(sp->values + point_b_first + 21);

        _MM_TRANSPOSE4_PD(z1fval, y1fval, x1fval, x2Compval);
        __m256d x1sval = _mm256_loadu_pd(sp->values + point_b_first + 4);

        __m256d y1sval = _mm256_loadu_pd(sp->values + point_b_first + 11);

        __m256d z1sval = _mm256_loadu_pd(sp->values + point_b_first + 18);

        __m256d null = _mm256_loadu_pd(sp->values + point_b_first + 25);

        _MM_TRANSPOSE4_PD(x1sval, y1sval, z1sval, null);

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

void multiplicateVector_values_AVX(MatrixValue *value, double *vect, double *result, int size, Task *task) {

    int fixedSizeX = task->nX + 2;
    int fixedSizeY = task->nY + 2;
    int fixedSizeZ = task->nZ + 2;

    // size without additional boundaries xz xy)
    int realSizeY = fixedSizeX;
    int realSizeZ = realSizeY * fixedSizeY;

    // offset for section start

    __m256d z1val = _mm256_set1_pd(value->z1);

    __m256d y1val = _mm256_set1_pd(value->y1);

    __m256d x1val = _mm256_set1_pd(value->x1);

    __m256d x2Compval = _mm256_set1_pd(value->x2Comp);

    for (int i = realSizeZ; i < size - realSizeZ; i += 4) {

        // first, get vect values
        __m256d z1f = _mm256_loadu_pd(vect + i - realSizeZ);
        __m256d y1f = _mm256_loadu_pd(vect + i - realSizeY);

        __m256d x1f = _mm256_loadu_pd(vect + i - 1);
        __m256d x2Comp = _mm256_loadu_pd(vect + i);
        __m256d x1s = _mm256_loadu_pd(vect + i + 1);

        __m256d y1s = _mm256_loadu_pd(vect + i + realSizeY);
        __m256d z1s = _mm256_loadu_pd(vect + i + realSizeZ);


        __m256d result1 = _mm256_mul_pd(z1f, z1val);
        __m256d result2 = _mm256_mul_pd(y1f, y1val);
        __m256d res_sum1 = _mm256_add_pd(result1, result2);

        result1 = _mm256_mul_pd(x1f, x1val);
        __m256d res_sum2 = _mm256_add_pd(res_sum1, result1);

        result1 = _mm256_mul_pd(x2Comp, x2Compval);
        res_sum1 = _mm256_add_pd(res_sum2, result1);

        result1 = _mm256_mul_pd(x1s, x1val);
        res_sum2 = _mm256_add_pd(res_sum1, result1);

        result1 = _mm256_mul_pd(y1s, y1val);
        res_sum1 = _mm256_add_pd(res_sum2, result1);

        result1 = _mm256_mul_pd(z1s, z1val);
        res_sum2 = _mm256_add_pd(res_sum1, result1);

        _mm256_storeu_pd(result + i, res_sum2);
    }
}
