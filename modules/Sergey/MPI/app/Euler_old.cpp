#include <iostream>
#include <mpi.h>
#include "Task.h"
#include "SparseMatrix.h"

using std::string;
using std::cout;
using std::endl;
using std::flush;

int main(int argc, char **argv) {

    // Timing variables
    double startTime, endTime;
    const int ROOT = 0;

    int threads = 0;
    int sizeP, rankP;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    Task task;

    double *vect;

    string outfilename;

    if (rankP == ROOT) {

        printf("Start\n");

        if (argc != 5) {
            printf("input data error!\n Format: setting.txt function.txt out.txt <threads>");
            exit(0);
        }


        string settingFile = argv[1];
        string functionFile = argv[2];
        outfilename = argv[3];
        threads = atoi(argv[4]);

        // Read task settings
        initTaskUsingFile(task, settingFile);
        setTimestep(task);

        // Init memory & read function file
        initMemoryReadDataMPI(vect, functionFile, task);
        boundaries_matrix_fix(vect, task.nX, task.nY, task.nZ);
    }

    MPI_Bcast(&threads, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.xStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.xEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.yStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.yEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.zStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.zEnd, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.sigma, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.bc, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepX, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepY, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.timeStepZ, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nX, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nY, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.nZ, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.fullVectSize, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.tStart, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.tFinish, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&task.dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);


    if (sizeP % 2 != 0) {
        cout << "Error! Wrong proc count" << endl;
        exit(0);
    }

    if (rankP == ROOT) {
        printf("Define variables\n");
    }
    int full_sizeX = task.nX + 2;    // define addition boundaries as default condition
    int sizeY = full_sizeX;          // size of one "y" cell
    int sizeZ = sizeY * task.nY;    // size of full plast (layer) of values

    int lineP = 2;                      // proc in line
    int rowP = sizeP / 2;
    if (sizeP % 2 != 0) {
        printf("%d don't divided by 2", sizeP);
        exit(0);
    }
    int proc_nX = full_sizeX;               // already added as default boundaries cell
    int proc_nY = task.nY / lineP + 2;  // there and after +2 as boundaries conditions
    if (task.nY % lineP != 0) {
        printf("nY %d don't divided by %d", task.nY, sizeP);
        exit(0);
    }
    int proc_nZ = task.nZ / rowP + 2;
    if (task.nZ % rowP != 0) {
        printf("nZ %d don't divided by rowP %d", task.nZ, rowP);
        exit(0);
    }
    int proc_vect_size = proc_nX * proc_nY * proc_nZ;




    int full_proc_sizeX = proc_nX;
    int proc_sizeY = full_proc_sizeX;
    int proc_sizeZ = proc_sizeY * proc_nY;  // this is size of one line of data

    int block_size_without_boundaries = proc_nX * (task.nY / lineP);
    int block_size = proc_sizeZ;

    // Generate data for send
    double *sendvect;

    if (rankP == ROOT) {
        sendvect = new double[sizeP * proc_vect_size];

        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
            sendvect[i] = -9000000000;
        }
//
        cout << "sizeP: " << sizeP << endl;
        cout << "proc_nX: " << proc_nX << endl;
        cout << "proc_nY: " << proc_nY << endl;
        cout << "proc_nZ: " << proc_nZ << endl;
        cout << "block_size: " << block_size << endl;
        cout << "proc_vect_size: " << proc_vect_size << endl;

//        printf("SEND!!!!!!\n");
//        for (int i = 0; i < task.fullVectSize; ++i) {
//            printf( "%2.15le\n", vect[i]);
//        }


        int position = 0;
        for (int p = 0; p < sizeP; ++p) {

            int y_offset = (p % lineP) * block_size_without_boundaries;
            int z_offset = (p / lineP) * lineP * (proc_nZ - 2) * block_size_without_boundaries;

            for (int z = 0; z < proc_nZ; ++z) {
                for (int i = 0, j = 0; i < block_size; ++i) {
                    position = p * proc_vect_size + z * block_size + i;
                    if (z == 0 || z == proc_nZ - 1) {   // for each process (top and bottom)
                        sendvect[position] = 0;
                    } else if (i / full_proc_sizeX == 0 || (i + full_proc_sizeX) / block_size == 1) { // first and last
                        sendvect[position] = 0;
                    } else {
                        sendvect[position] = vect[z_offset + y_offset + (z - 1) * sizeZ + j];
                        ++j;
                    }
                }
            }
        }

        std::cout << "Flag 2" << std::endl;
        fflush(stdout);
//        printf("=========\n");
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            printf( "%2.15le\n", sendvect[i]);
//        }

    }

    int *displs = new int[sizeP];
    int *sendcounts = new int[sizeP];

    double **proc_vect = new double *[2];
    proc_vect[0] = new double[proc_vect_size];
    proc_vect[1] = new double[proc_vect_size];

    for (int l = 0; l < sizeP; ++l) {
        displs[l] = proc_vect_size * l;
        sendcounts[l] = proc_vect_size;     // we send half part
    }


    MPI_Scatterv(sendvect, sendcounts, displs, MPI_DOUBLE,
                 proc_vect[0], proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);




    MPI_Barrier(MPI_COMM_WORLD);
    // value for the matrix
    MatrixValue matrixValue;
    matrixValue.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrixValue.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrixValue.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrixValue.x2Comp = (1 - 2 * matrixValue.x1 - 2 * matrixValue.y1 - 2 * matrixValue.z1);

    // init and fill sparseMatrix
    SparseMatrix spMat;
    int sparseMatrixSize = 9 * proc_nX * proc_nY * proc_nZ;

    spMatrixInit(spMat, sparseMatrixSize, proc_vect_size, threads);

    fillMatrix3d6Expr_wo_boundaries(spMat, matrixValue, proc_nX, proc_nY, proc_nZ);


    int prevTime = 1;
    int currTime = 0;

    startTime = MPI_Wtime();

    int send_vect_size = proc_nZ * proc_nX; // Because nY don't involved (left right)
    double *left_vect = new double[send_vect_size];
    double *right_vect = new double[send_vect_size];

    double *left_vect_get = new double[send_vect_size];
    double *right_vect_get = new double[send_vect_size];

    for (double j = 0; j < task.tFinish; j += task.dt) {

        // Before each iteration - swap data
        // 1 and 2 iter - IT'S TIME TO REPACK
//
//        for (int z = 0; z < proc_nZ; ++z) {
//            for (int i = 0; i < proc_nX; ++i) {
//                //z * proc_sizeZ + i
//                left_vect[z * proc_nX + i] = proc_vect[currTime][z * proc_sizeZ + proc_sizeY + i];
//                //z * proc_sizeZ + (proc_nY - 1) * proc_sizeY + i
//                right_vect[z * proc_nX + i] = proc_vect[currTime][z * proc_sizeZ + (proc_nY - 2) * proc_sizeY + i];
//            }
//        }
//
//        int left_proc = (rankP / lineSizeP) * lineSizeP + (rankP - 1 + lineSizeP) % lineSizeP;
//        int right_proc = (rankP / lineSizeP) * lineSizeP + (rankP + 1 + lineSizeP) % lineSizeP;
//
//        MPI_Sendrecv(left_vect, send_vect_size, MPI_DOUBLE, left_proc, 1,
//            left_vect_get, send_vect_size, MPI_DOUBLE, right_proc, 1, MPI_COMM_WORLD, &status);
////
//        MPI_Sendrecv(right_vect, send_vect_size, MPI_DOUBLE, right_proc, 1,
//            right_vect_get, send_vect_size, MPI_DOUBLE, left_proc, 1, MPI_COMM_WORLD, &status);
////
//        for (int z = 0; z < proc_nZ; ++z) {
//            for (int i = 0; i < proc_nX; ++i) {
//                proc_vect[currTime][z * proc_sizeZ + i] = left_vect_get[z * proc_nX + i];
//                proc_vect[currTime][z * proc_sizeZ + (proc_nY - 1) * proc_sizeY + i] =
//                    right_vect_get[z * proc_nX + i];
//            }
//        }
////
//        int top_proc = (rankP - lineSizeP + sizeP) % sizeP;
//        int bottom_proc = (rankP + lineSizeP) % sizeP;
//
////        cout << "top_proc: " << top_proc << endl;
////        cout << "bottom_proc: " << bottom_proc << endl;
//
//        // 1 and 2 iter - in one field in the memory
//        send_vect_size = proc_nY * proc_nX;
//        // top
//        MPI_Sendrecv(proc_vect[currTime] + proc_sizeZ, send_vect_size, MPI_DOUBLE, top_proc, 2,
//                     proc_vect[currTime], send_vect_size, MPI_DOUBLE, bottom_proc, 2, MPI_COMM_WORLD, &status);
//
//        // bottom
//        int offset = proc_vect_size - send_vect_size;
//        MPI_Sendrecv(proc_vect[currTime] + offset - proc_sizeZ, send_vect_size, MPI_DOUBLE, bottom_proc, 2,
//                     proc_vect[currTime] + offset, send_vect_size, MPI_DOUBLE, top_proc, 2, MPI_COMM_WORLD, &status);
//
////
//        int boundaries_offset = proc_sizeY;
//
//        // BOUNDARIES!
//        for (int k = 0; k < proc_nY * proc_nZ; ++k) {
//            proc_vect[currTime][k * boundaries_offset] = proc_vect[currTime][k * boundaries_offset + 1];
//            proc_vect[currTime][k * boundaries_offset + proc_nX - 1] =
//                    proc_vect[currTime][k * boundaries_offset + proc_nX - 2];
//        }

//        if (rankP != 2) {
////            multiplicateVector(spMat, proc_vect[currTime], proc_vect[prevTime], proc_vect_size);
//            ;
//        } else {
//            for (int i = 0; i < proc_vect_size; ++i) {
//                proc_vect[prevTime][i] == -999999999;
//            }
//        }
//        if (rankP == 5) {
//            for (int i = 0; i < proc_vect_size; ++i) {
//                proc_vect[prevTime][i] = -999999999;
//                proc_vect[currTime][i] = -999999999;
//            }
//        }

        prevTime = (prevTime + 1) % 2;
        currTime = (currTime + 1) % 2;
    }
    std::cout << "Flag" << std::endl;
    fflush(stdout);

    // And now we should scatter it in one vector...

    sendvect = new double[sizeP * proc_vect_size];

    for (int i = 0; i < sizeP * proc_vect_size; ++i) {
        sendvect[i] = -9000000000;
    }

    MPI_Gather(proc_vect[prevTime], proc_vect_size, MPI_DOUBLE,
               sendvect, proc_vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

//    if (rankP == ROOT) {
//        cout << "print sendvect vect" << endl;
//        cout.precision(15);
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            cout << i << "  " << sendvect[i] << endl;
//        }
//        cout << "end" << endl;
//    }in_Vect_offset



    if (rankP == ROOT) {
        for (int j = 0; j < task.fullVectSize; ++j) {
            vect[j] = -330000;
        }

//        printf("=========\n");
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            printf( "%2.15le\n", sendvect[i]);
//        }

        cout << "full vect size: " << task.fullVectSize << endl;
        cout << "nX: " << task.nX << endl;
        cout << "nY: " << task.nY << endl;
        cout << "nZ: " << task.nZ << endl;
        cout << endl;
        cout << "sizeP: " << sizeP << endl;
        cout << "proc_nX: " << proc_nX << endl;
        cout << "proc_nY: " << proc_nY << endl;
        cout << "proc_nZ: " << proc_nZ << endl;
        cout << "block_size: " << block_size << endl;
        cout << "proc_vect_size: " << proc_vect_size << endl;
        cout << "sizeZ: " << sizeZ << endl;
        cout << "lineSizeP: " << lineP << endl;


        int temp;
        int position = 0;
        for (int p = 0; p < sizeP; ++p) {

            int y_offset = (p % lineP) * block_size_without_boundaries;
            int z_offset = (p / lineP) * lineP * (proc_nZ - 2) * block_size_without_boundaries;

            for (int z = 0; z < proc_nZ; ++z) {
                for (int i = 0, j = 0; i < block_size; ++i) {
                    position = p * proc_vect_size + z * block_size + i;
                    if (z == 0 || z == proc_nZ - 1) {   // for each process (top and bottom)
                        ;
                    } else if (i / full_proc_sizeX == 0 || (i + full_proc_sizeX) / block_size == 1) { // first and last
                        ;
                    } else {
//                        if (temp > z_offset + y_offset + (z - 1) * sizeZ + j) {
//                            printf("****\nblock_size %d\nbswithout_boundaries %d\n", block_size, block_size_without_boundaries);
//                            printf("===\np %d\ny %d\nzo %d\nz %d\n", position, y_offset, z_offset, z);
//                        }
//                        printf("%d\n", temp = z_offset + y_offset + (z - 1) * sizeZ + j);
                        vect[z_offset + y_offset + (z - 1) * sizeZ + j] = sendvect[position];
                        ++j;
                    }
                }
            }
        }


        printf("Finish\n");
    }


    if (rankP == ROOT) {
        endTime = MPI_Wtime();
        printf("Run time %.15lf\n", endTime - startTime);


//        printf("GET!!!!\n");
//        for (int i = 0; i < task.fullVectSize; ++i) {
//            printf( "%2.15le\n", vect[i]);
//        }

        // Output
        FILE *outfile = fopen(outfilename.c_str(), "w");
//        FILE *outfile = fopen("../../result/Sergey/Sergey_Euler_MPI_full_test2.txt", "w");

        for (int i = 0; i < task.fullVectSize; ++i) {
            if (i % (task.nX + 2) != 0 && i % (task.nY + 2) != task.nZ + 1)
                fprintf(outfile, "%2.15le\n", vect[i]);
        }
//        for (int i = 1; i <= task.nX; ++i) {
//            fprintf(outfile, "%2.15le\n", getVectorValue(vect, i, 0, 0, task));
//
//        }
    }
    MPI_Finalize();
}
