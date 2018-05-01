#include <iostream>
#include <omp.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;

double getVectorValue(double *vect, int x, int y, int z, Task task);

int main(int argc, char **argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;
    int threads = 0;

    if (argc != 5) {
        printf("input data error!\n Format: setting.txt function.txt out.txt");
        exit(0);
    }


    string settingFile = argv[1];
    string functionFile = argv[2];
    string outfilename = argv[3];
    threads = atoi(argv[4]);


    // File variables
//    string functionFile = "../../../../../initial/function2.txt";
    //string functionFile = "function2.txt";
//    string functionFile = "../../../../../initial_for_tests/function5.txt";
//    string functionFile = "../../initial_test/function.txt";
//    string settingFile = "../../../../../initial/setting2.ini";
    //string settingFile = "setting2.ini";
//    string settingFile = "../../../../../initial_for_tests/setting5.ini";
//    string settingFile = "../../initial_test/setting.ini";

    // Read task settings
    Task task;
    initTaskUsingFile(task, settingFile);
    setTimestep(task);

    string consoleInput = "";

    // Init memory & read function file
    double **vect;
    initMemoryReadData(vect, functionFile, task);

    // vector time-index for loop
    prevTime = 0;
    currTime = 1;

    boundaries_matrix_fix(vect[prevTime], task.nX, task.nY, task.nZ);

    // value for the matrix
    MatrixValue matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.stepX * task.stepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.stepY * task.stepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.stepZ * task.stepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 9 * task.nX * task.nY * task.nZ;

    spMatrixInit(spMat, sparseMatrixSize, task.fullVectSize, threads);
    fillMatrix3d6Expr(spMat, matrixValue, task.nX, task.nY, task.nZ);



    // Calculating
    time_S = omp_get_wtime();

    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVector(spMat, vect[prevTime], vect[currTime], task.fullVectSize);
        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E - time_S);
    printf("On %d threads\n", threads);

    // Output
//    string outfilename = "../../../../../result/Sergey-N/Openmp_Euler_1.txt";

    //string outfilename = "accurancy_test/stepX_" + std::to_string(task.stepX) + "_stepT_" +
    //        std::to_string(task.dt)+ ".txt";
    //string outfilename = "accurancy_test/res.txt";

    FILE *outfile = fopen(outfilename.c_str(), "w");

    for (int i = 0; i < task.fullVectSize; ++i) {
        if (i % (task.nX + 2) != 0 && i % (task.nX + 2) != task.nX + 1)
            fprintf(outfile, "%2.15le\n", vect[prevTime][i]);
    }

}