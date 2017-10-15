#include "omp.h"
#include "SparseMatrix.h"
#include "Task.h"


int main(int argc, char **argv) {

    // Timing variables
    double time_S, time_E;
    int prevTime, currTime;
    int threads = 0;

    if (argc != 5) {
        printf("input data error!\n Format: setting.txt function.txt out.txt <threads>\n");
        return 0;
    }


    char* settingFile = argv[1];
    char*  functionFile = argv[2];
    char* outfilename = argv[3];
    threads = atoi(argv[4]);

    Task task;
    // Read task settings
    if (read_settings(settingFile, &task) != 0) {
        printf("Settings reading error\n");
        return -1;
    }
    // Init memory & read function file
    task.fullVectSize = (task.nX + 2) * (task.nY + 2) * (task.nZ + 2);
    double* vect = (double *)calloc((size_t) task.fullVectSize, sizeof(double));
    double* next_vect = (double *)calloc((size_t) task.fullVectSize, sizeof(double));
    double* tmp_vect;
    int errors = 0;
    errors = initFunctionData_forAdditionalXYZ(functionFile, vect, &task);
    if (errors != 0) {
        printf("Function file reading error\n");
        return -2;
    }

    // vector time-index for loop
    prevTime = 0;
    currTime = 1;


    // Boundaries fix
    errors = boundariesFix_forAdditionalXYZ(vect, &task);
    if (errors != 0) {
        printf("Boundaries error\n");
        return -3;
    }

    // value for the matrix
    MatrixVal matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 7 * (task.nX + 2) * (task.nY + 2) * (task.nZ + 2);

    spMatrixInit(&spMat, sparseMatrixSize, task.fullVectSize, threads);
    fillMatrix3d6Expr_wo_boundaries_for_xyz(&spMat, &matrixValue, &task);

    // Calculating
    time_S = omp_get_wtime();
//tFinish
    for (double j = 0; j < task.tFinish; j += task.dt) {
        multiplicateVectorAVXColumn5_shuffle(&spMat, vect, next_vect, task.fullVectSize);
        boundariesFix_forAdditionalXYZ(next_vect, &task);
        tmp_vect = vect;
        vect = next_vect;
        next_vect = tmp_vect;
    }
    time_E = omp_get_wtime();
    printf("Run time %.15lf\n", time_E - time_S);
    printf("On %d threads\n", threads);

    FILE *outfile = fopen(outfilename, "w");

    int realSizeX = task.nX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * (task.nY + 2);

    int offset;
    for (int z = 1; z < task.nZ + 1; ++z) {
        for (int y = 1; y < task.nY +1; ++y) {
            offset = z * realSizeZ + y * realSizeY;
            for (int x = 1; x < task.nX + 1; ++x) {
                fprintf(outfile, "%2.15le\n", vect[offset+x]);
            }
        }
    }

}