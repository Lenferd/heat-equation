//
// Created by lenferd on 03.04.17.
//

#include <iostream>
#include <omp.h>
#include <cmath>
#include "Task_old.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);

double normVect(double *&vect1, double *&vect2, int size) {
    double norm = 0;
    norm = fabs(vect1[0] - vect2[0]);
    for (int h = 0; h < size; h++) {
        if (fabs(vect1[h] - vect2[h]) > norm)
            norm = fabs(vect1[h] - vect2[h]);
    }
    return norm;
}

int main(int argc, char** argv) {

    double eps = 1e-6;

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;

    int threads = 0;

    if (argc != 5) {
        printf("Input data error!\nFormat: setting.txt function.txt out.txt <threads>\n");
        exit(0);
    }

    string settingFile = argv[1];
    string functionFile = argv[2];
    string outfilename = argv[3];
    threads = atoi(argv[4]);

    // Read task settings
    Task task;
    initTaskUsingFile(task, settingFile);
    setTimestep(task);

    // Init memory & read function file
    double** vect;
    initMemoryReadData(vect, functionFile, task);

    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    boundaries_matrix_fix(vect[prevTime], task.nX, task.nY, task.nZ);

    // value for the matrix
    MatrixValue matrixValue;

    double addit_up_value_x = (task.sigma * task.dt) / (task.stepX * task.stepX);
    double addit_up_value_y = (task.sigma * task.dt) / (task.stepY * task.stepY);
    double addit_up_value_z = (task.sigma * task.dt) / (task.stepZ * task.stepZ);
    printf("x %lf\n", addit_up_value_x);
    printf("y %lf\n", addit_up_value_y);
    printf("z %lf\n", addit_up_value_z);

    double addit_dw_value = (1 + 2 * addit_up_value_x + 2 * addit_up_value_y + 2 * addit_up_value_z); // a_ii in jacobi method
    printf("dw %lf\n", addit_dw_value);

    double addit_add_value = 1 / addit_dw_value;

    matrixValue.x1 = addit_up_value_x / addit_dw_value; // this value must be with minus, but we already apply it to
                                                        // addit_up_values
    matrixValue.y1 = addit_up_value_y / addit_dw_value;
    matrixValue.z1 = addit_up_value_z / addit_dw_value;
    matrixValue.x2Comp = 0;

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMat, sparseMatrixSize, task.fullVectSize, threads);
    fillMatrix3d6Expr(spMat, matrixValue, task.nX, task.nY, task.nZ);


    // prevTime is b vect now
    double* const_vect = new double[task.fullVectSize];

    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {

//        #pragma omp parallel for
        for (int i = 0; i < task.fullVectSize; i++) {
            const_vect[i] = vect[prevTime][i] * addit_add_value ;
        }
        int iteration_counter = 0;

        do {
            multiplicateVector(spMat, vect[prevTime], vect[currTime], task.fullVectSize);

//            #pragma omp parallel for
            for (int i = 0; i < task.fullVectSize; i++) {
                vect[currTime][i] += const_vect[i];
            }

            prevTime = (prevTime + 1) % 2;
            currTime = (currTime + 1) % 2;

            iteration_counter++;
        } while (normVect(vect[prevTime], vect[currTime], task.fullVectSize) > eps);  // vect[prevTime] = k + 1 now

    }

    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);


    // Output
    FILE *outfile = fopen(outfilename.c_str(), "w");

//    double outData;
    for (int i = 0; i < task.fullVectSize; ++i) {
        if (i % (task.nX + 2) != 0 && i % (task.nX + 2) != task.nX + 1)
            fprintf(outfile, "%2.15le\n", vect[0][i]);
    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}

