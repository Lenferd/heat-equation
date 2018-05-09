//
// Created by lenferd on 28.03.17.
//

#include <iostream>
#include <omp.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);
int main(int argc, char** argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;

    int threads;

    if (argc != 5) {
        printf("input data error!\n Format: setting.txt function.txt out.txt <threads>\n");
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
    // initMemoryReadData(vect, functionFile, task);
    initMemoryReadData_for_additional_xyz(vect, functionFile, task);

    double* vectK1 = new double[task.fullVectSize];
    double* vectK2 = new double[task.fullVectSize];
    double* vectK3 = new double[task.fullVectSize];
    double* vectK4 = new double[task.fullVectSize];


    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    // boundaries_matrix_fix(vect[prevTime], task.nX, task.nY, task.nZ);
    boundaries_matrix_fix_for_xyz(vect[0], task.nX, task.nY, task.nZ);
    /***
     * K1 Vector
     */
    // value for the matrix
    MatrixValue matrixValueK1;
    matrixValueK1.x1 = (task.sigma) / (task.stepX * task.stepX);
    matrixValueK1.y1 = (task.sigma) / (task.stepY * task.stepY);
    matrixValueK1.z1 = (task.sigma) / (task.stepZ * task.stepZ);
    matrixValueK1.x2Comp = (- 2 * matrixValueK1.x1 - 2 * matrixValueK1.y1 - 2 * matrixValueK1.z1);

//    printf("\nk1\n");
//    printf("x1 %lf\t", matrixValueK1.x1);
//    printf("y1 %lf\t", matrixValueK1.y1);
//    printf("z1 %lf\t", matrixValueK1.z1);
//    printf("x2C %lf\t", matrixValueK1.x2Comp);

    // init and fill sparseMatrix
    SparseMatrix spMatK1;
    // int sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;
    int sparseMatrixSize = 7 * (task.nX + 2) * (task.nY + 2) * (task.nZ + 2);

    spMatrixInit(spMatK1, sparseMatrixSize, task.fullVectSize, threads);
    // fillMatrix3d6Expr(spMatK1, matrixValueK1, task.nX, task.nY, task.nZ);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(spMatK1, matrixValueK1, task.nX, task.nY, task.nZ);
    /***
    * K2
    */
    // value for the matrix
    MatrixValue matrixValueK2;
    matrixValueK2.x1 = (task.sigma * task.dt * 0.5) / (task.stepX * task.stepX);
    matrixValueK2.y1 = (task.sigma * task.dt * 0.5) / (task.stepY * task.stepY);
    matrixValueK2.z1 = (task.sigma * task.dt * 0.5) / (task.stepZ * task.stepZ);
    matrixValueK2.x2Comp = (1 - 2 * matrixValueK2.x1 - 2 * matrixValueK2.y1 - 2 * matrixValueK2.z1);

//    printf("\nk2\n");
//    printf("x1 %lf\t", matrixValueK2.x1);
//    printf("y1 %lf\t", matrixValueK2.y1);
//    printf("z1 %lf\t", matrixValueK2.z1);
//    printf("x2C %lf\t", matrixValueK2.x2Comp);
    // init and fill sparseMatrix
    SparseMatrix spMatK2;
    // sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK2, sparseMatrixSize, task.fullVectSize, threads);
    // fillMatrix3d6Expr(spMatK2, matrixValueK2, task.nX, task.nY, task.nZ);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(spMatK2, matrixValueK2, task.nX, task.nY, task.nZ);

    /***
    * K3
    */
    // value for the matrix
    MatrixValue matrixValueK3;
    matrixValueK3.x1 = (task.sigma * task.dt * 0.5) / (task.stepX * task.stepX);
    matrixValueK3.y1 = (task.sigma * task.dt * 0.5) / (task.stepY * task.stepY);
    matrixValueK3.z1 = (task.sigma * task.dt * 0.5) / (task.stepZ * task.stepZ);
    matrixValueK3.x2Comp = (-2 * matrixValueK3.x1 - 2 * matrixValueK3.y1 - 2 * matrixValueK3.z1);

//    printf("\nk3\n");
//    printf("x1 %lf\t", matrixValueK3.x1);
//    printf("y1 %lf\t", matrixValueK3.y1);
//    printf("z1 %lf\t", matrixValueK3.z1);
//    printf("x2C %lf\t", matrixValueK3.x2Comp);
    // init and fill sparseMatrix
    SparseMatrix spMatK3;
    // sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK3, sparseMatrixSize, task.fullVectSize, threads);
    // fillMatrix3d6Expr(spMatK3, matrixValueK3, task.nX, task.nY, task.nZ);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(spMatK3, matrixValueK3, task.nX, task.nY, task.nZ);
    /***
    * K4 Vector
    */
    // value for the matrix
    MatrixValue matrixValueK4;
    matrixValueK4.x1 = (task.sigma * task.dt) / (task.stepX * task.stepX);
    matrixValueK4.y1 = (task.sigma * task.dt) / (task.stepY * task.stepY);
    matrixValueK4.z1 = (task.sigma * task.dt) / (task.stepZ * task.stepZ);
    matrixValueK4.x2Comp = (- 2 * matrixValueK4.x1 - 2 * matrixValueK4.y1 - 2 * matrixValueK4.z1);

//    printf("\nk4\n");
//    printf("x1 %lf\t", matrixValueK4.x1);
//    printf("y1 %lf\t", matrixValueK4.y1);
//    printf("z1 %lf\t", matrixValueK4.z1);
//    printf("x2C %lf\t", matrixValueK4.x2Comp);

    // init and fill sparseMatrix
    SparseMatrix spMatK4;
    // sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMatK4, sparseMatrixSize, task.fullVectSize, threads);
    // fillMatrix3d6Expr(spMatK4, matrixValueK4, task.nX, task.nY, task.nZ);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(spMatK4, matrixValueK4, task.nX, task.nY, task.nZ);
    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
        // multiplicateVector(spMatK1, vect[prevTime], vectK1, task.fullVectSize);
        // multiplicateVector(spMatK2, vectK1, vectK2, task.fullVectSize);
        // multiplicateVectorRunge(spMatK3, vectK2, vectK1, vectK3, task.fullVectSize);
        // multiplicateVectorRunge(spMatK4, vectK3, vectK1, vectK4, task.fullVectSize);

        multiplicateVectorAVXColumn5_shuffle(&spMatK1, vect[prevTime], vectK1, task.fullVectSize);
        boundaries_matrix_fix_for_xyz(vectK1, task.nX, task.nY, task.nZ);
        multiplicateVectorAVXColumn5_shuffle(&spMatK2, vectK1, vectK2, task.fullVectSize);
        boundaries_matrix_fix_for_xyz(vectK2, task.nX, task.nY, task.nZ);
        // multiplicateVectorAVXColumn5_shuffle(&spMatK3, vectK2, vectK3, task.fullVectSize);
        multiplicateVectorRunge(spMatK3, vectK2, vectK1, vectK3, task.fullVectSize);
        boundaries_matrix_fix_for_xyz(vectK3, task.nX, task.nY, task.nZ);
        // multiplicateVectorAVXColumn5_shuffle(&spMatK4, vectK3, vectK4, task.fullVectSize);
        multiplicateVectorRunge(spMatK4, vectK3, vectK1, vectK4, task.fullVectSize);
        boundaries_matrix_fix_for_xyz(vectK4, task.nX, task.nY, task.nZ);
        #pragma omp parallel for
        for (int i = 0; i < task.fullVectSize; ++i) {
//            if (i % (task.nX + 2) != 0 && i % (task.nX + 2) != task.nX + 1) {
                vect[currTime][i] =
                        vect[prevTime][i] + task.dt / 6 * (vectK1[i] + 2 * vectK2[i] + 2 * vectK3[i] + vectK4[i]);
//            } else if (i % (task.nX)){
//                vect[currTime][i] = vect[prevTime][i];
//            }
        }
        boundaries_matrix_fix_for_xyz(vect[currTime], task.nX, task.nY, task.nZ);
//        boundaries_matrix_fix(vect[currTime], vect[prevTime], task.nX, task.nY, task.nZ);

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E-time_S);
    printf("On %d threads\n", threads);

    // Output
    FILE *outfile = fopen(outfilename.c_str(), "w");

    int realSizeX = task.nX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * (task.nY + 2);

    int offset;
    for (int z = 1; z < task.nZ + 1; ++z) {
        for (int y = 1; y < task.nY +1; ++y) {
            offset = z * realSizeZ + y * realSizeY;
            for (int x = 1; x < task.nX + 1; ++x) {
                fprintf(outfile, "%2.15le\n", vect[prevTime][offset+x]);
            }
        }
    }
}

double getVectorValue(double *vect, int x, int y, int z, Task task) {
    return vect[x + (task.nX + 2) * y + (task.nX+2)*task.nY*z];
}
