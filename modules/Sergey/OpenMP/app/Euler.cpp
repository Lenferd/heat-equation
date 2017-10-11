#include <iostream>
#include <cstdlib>
#include "Task.h"
#include "omp.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);

int main(int argc, char **argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;
    int threads = 0;

    if (argc != 5) {
        printf("input data error!\n Format: setting.txt function.txt out.txt <threads>\n");
        return 0;
    }


    string settingFile = argv[1];
    string functionFile = argv[2];
    string outfilename = argv[3];
    threads = atoi(argv[4]);

    // Read task settings
    Task task;
    initTaskUsingFile(task, settingFile);
    setTimestep(task);

    string consoleInput = "";

    // Init memory & read function file
    double **vect;
    initMemoryReadData_for_additional_xyz(vect, functionFile, task);

    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    boundaries_matrix_fix_for_xyz(vect[0], task.nX, task.nY, task.nZ);

    // value for the matrix
    MatrixValue matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.stepX * task.stepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.stepY * task.stepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.stepZ * task.stepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 7 * (task.nX + 2) * (task.nY + 2) * (task.nZ + 2);

    spMatrixInit(spMat, sparseMatrixSize, task.fullVectSize, threads);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(spMat, matrixValue, task.nX, task.nY, task.nZ);

    // Calculating
    time_S = omp_get_wtime();
//tFinish
    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(spMat, vect[prevTime], vect[currTime], task.fullVectSize);
        boundaries_matrix_fix_for_xyz(vect[currTime], task.nX, task.nY, task.nZ);
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E - time_S);
    printf("On %d threads\n", threads);

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