//
// Created by lenferd on 08.09.17.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "Task.h"
#include "SparseMatrix.h"

#define ROOT 0

int main(int argc, char **argv) {

    // time variables
    double startTime, endTime;

    // task settings struct
    Task task;

    // MPI Stuff
    int sizeP, rankP, omp_threads;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &sizeP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankP);

    double *vect = NULL;
    char* out_filename;

    if (rankP == ROOT) {
        int errors = 0;
        printf("Start\n");

        if (argc != 5) {
            printf("input data error!\n Format: setting.txt function.txt out.txt <omp_threads>");
            exit(0);
        }

        char* setting_file = argv[1];
        char* function_file = argv[2];
        out_filename = argv[3];
        omp_threads = atoi(argv[4]);

        // Read task settings
        if (read_settings(setting_file, &task) != 0) {
            printf("Settings reading error\n");
            return -1;
        }

        // Init memory & read fuction file
        task.fullVectSize = (size_t)(task.nX + 2) * (task.nY) * (task.nZ);
        vect = (double *)calloc(task.fullVectSize, sizeof(double));

        errors = initFunctionMPI(function_file, vect, &task);
        if (errors != 0) {
            printf("Function file reading error\n");
            return -2;
        }

        // Boundaries fix
        errors = boundaries_conditions(vect, &task);
        if (errors != 0) {
            printf("Boundaries error\n");
            return -3;
        }
    }

    MPI_Bcast(&omp_threads, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

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
        printf("Processors count cannot be divided by 2\n");
    }

    if (rankP == ROOT) {
        printf("Define variables\n");
    }

    // Default matrix (not already splitted)
    task.sizeY = task.nX + 2;  // define addition boundaries as default condition
    task.sizeZ = task.sizeY * task.nY;   // size of full plast (layer) of values

    int lineP = 2;                      // proc in line
    int rowP = sizeP / 2;

    if (sizeP % 2 != 0) {
        printf("%d don't divided by 2", sizeP);
        exit(0);
    }

    // Proc settings
    ProcSettings proc;

    proc.nX = task.sizeY;               // already added as default boundaries cell
    proc.nY = task.nY / lineP + 2;  // there and after +2 as boundaries conditions
    if (task.nY % lineP != 0) {
        printf("nY %d don't divided by %d", task.nY, sizeP);
        exit(0);
    }
    proc.nZ = task.nZ / rowP + 2;
    if (task.nZ % rowP != 0) {
        printf("nZ %d don't divided by rowP %d", task.nZ, rowP);
        exit(0);
    }

    proc.vect_size = (size_t)proc.nX * proc.nY * proc.nZ;
    proc.sizeY = proc.nX;
    proc.sizeZ = proc.sizeY * proc.nY;  // this is size of one line of data


//    int block_size_without_boundaries = proc.nX * (task_settings.nY / lineP);
//    int block_size = proc_sizeZ;

    // Generate data for send
    ///// SEND THIS SHIT //////
    double *send_vect;

    if (rankP == ROOT) {
        send_vect = (double*)calloc((size_t)sizeP * proc.vect_size, sizeof(double));
//        for (int i = 0; i < sizeP * proc_vect_size; ++i) {
//            send_vect[i] = -9000000000;
//        }

        for (int p = 0; p < sizeP; ++p) {

            // We searching offset for original block
            // That's means, to the y offset we should get original
            int y_block_size = task.sizeY * (task.nY / lineP); // half xy (one proc layer)
            int y_offset = (p % lineP) * y_block_size;
            int z_block_size = task.sizeZ * (task.nZ / rowP);
            int z_offset = (p / lineP) * z_block_size;

            for (int z = 1; z < proc.nZ - 1; ++z) {
                for (int y = 1; y < proc.nY - 1; ++y) {
                    for (int x = 0; x < proc.nX; ++x) {
                        send_vect[p * proc.vect_size + x + (y * proc.sizeY) + (z * proc.sizeZ)] =
                                vect[y_offset + z_offset +
                                        x + ((y - 1) * task.sizeY) + ((z - 1) * task.sizeZ)];
                    }
                }
            }
        }
    }

    int *displ = (int*)calloc(sizeP, sizeof(int));
    int *sendcounts = (int*)calloc(sizeP, sizeof(int));

    double *proc_vect = (double*)calloc(proc.vect_size, sizeof(double));
    double *proc_vect_next = (double*)calloc(proc.vect_size, sizeof(double));
    double *tmp;

    for (int i = 0; i < sizeP; ++i) {
        displ[i] = (int)proc.vect_size * i;
        sendcounts[i] = (int)proc.vect_size;
    }

    MPI_Scatterv(send_vect, sendcounts, displ, MPI_DOUBLE,
                 proc_vect, (int)proc.vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    ////// CALCULATE THIS SHIT /////

    if (rankP == ROOT) {
        printf("Define matrix\n");
    }

    // Matrix fill
    MatrixVal matrix_val;
    matrix_val.x1 = (task.sigma * task.dt) / (task.timeStepX * task.timeStepX);
    matrix_val.y1 = (task.sigma * task.dt) / (task.timeStepY * task.timeStepY);
    matrix_val.z1 = (task.sigma * task.dt) / (task.timeStepZ * task.timeStepZ);
    matrix_val.x2Comp = (1 - 2 * matrix_val.x1 - 2 * matrix_val.y1 - 2 * matrix_val.z1);

    // Init and fill sparseMatrix
    SparseMatrix spMatrix;
    int sp_matrix_size = 9 * proc.nX * proc.nY * proc.nZ;
    spMatrixInit(&spMatrix, sp_matrix_size, (int)proc.vect_size, omp_threads);

    fillMatrix_3d_MPI_woBound(&spMatrix, &matrix_val, &proc);


    if (rankP == ROOT) {
        printf("Define sendrecv info\n");
    }

    // MPI Sendrecv information
    int send_vect_size_lr = proc.nZ * proc.nX; // left-right
    int send_vect_size_td = proc.nY * proc.nX;
    double *left_vect = (double*)malloc(send_vect_size_lr * sizeof(double));
    double *right_vect = (double*)malloc(send_vect_size_lr * sizeof(double));

    double *left_vect_get = (double*)malloc(send_vect_size_lr * sizeof(double));
    double *right_vect_get = (double*)malloc(send_vect_size_lr * sizeof(double));

    for (double j=0; j < task.tFinish; j+=task.dt) {
        // Before each iteration - swap data

        MPI_Barrier(MPI_COMM_WORLD);
        // 1 and 2 iter - IT'S TIME TO REPACK

        /// LEFT-RIGHT (xz)

        for (int z = 0; z < proc.nZ; ++z) {
            for (int i = 0; i < proc.nX; ++i) {
                //z * proc_sizeZ + i
                left_vect[z * proc.nX + i] = proc_vect[z * proc.sizeZ + proc.sizeY + i];
                //z * proc_sizeZ + (proc_nY - 1) * proc_sizeY + i
                right_vect[z * proc.nX + i] = proc_vect[(z+1) * proc.sizeZ - (2 * proc.sizeY) + i];
            }
        }

//        for (int p = 1, i=0; p < proc.vect_size; p+=proc.sizeY, ++i) {
//            left_vect[i] = proc_vect[p];
//            right_vect[i] = proc_vect[p + proc.sizeY - 3];  // 2 is imaginary boundaries
//        }

        // (rankP / lineSizeP) * lineSizeP find out start line
        int left_proc = (rankP / lineP) * lineP + (rankP - 1 + lineP) % lineP; // + lineP for 0-1
        int right_proc = (rankP / lineP) * lineP + (rankP + 1 + lineP) % lineP;

        MPI_Sendrecv(left_vect, send_vect_size_lr, MPI_DOUBLE, left_proc, 1,
                     right_vect_get, send_vect_size_lr, MPI_DOUBLE, right_proc, 1, MPI_COMM_WORLD, &status);

        MPI_Sendrecv(right_vect, send_vect_size_lr, MPI_DOUBLE, right_proc, 2,
                     left_vect_get, send_vect_size_lr, MPI_DOUBLE, left_proc, 2, MPI_COMM_WORLD, &status);

        for (int z = 0; z < proc.nZ; ++z) {
            for (int i = 0; i < proc.nX; ++i) {
                //z * proc_sizeZ + i
                proc_vect[z * proc.sizeZ + i] = left_vect_get[z * proc.nX + i];
                        //z * proc_sizeZ + (proc_nY - 1) * proc_sizeY + i
                proc_vect[(z+1) * proc.sizeZ - (1 * proc.sizeY) + i] = right_vect_get[z * proc.nX + i];
            }
        }

        for (int p = 0, i=0; p < proc.vect_size; p+=proc.sizeY, ++i) {
            proc_vect[p] = left_vect_get[i];
            proc_vect[p + proc.sizeY - 1] = right_vect_get[i];  // -1 for end in prev vector
        }
        /// TOP-DOWN (xy)

        int top_proc = (rankP - lineP + sizeP) % sizeP;
        int bottom_proc = (rankP + lineP) % sizeP;
        int offset = (int)proc.vect_size - send_vect_size_td;

//         top
        MPI_Sendrecv(proc_vect + proc.sizeZ, send_vect_size_td, MPI_DOUBLE, top_proc, 3,
                     proc_vect + offset, send_vect_size_td, MPI_DOUBLE, bottom_proc, 3, MPI_COMM_WORLD, &status);
//
//         bottom
        MPI_Sendrecv(proc_vect + offset - proc.sizeZ, send_vect_size_td, MPI_DOUBLE, bottom_proc, 4,
                     proc_vect, send_vect_size_td, MPI_DOUBLE, top_proc, 4, MPI_COMM_WORLD, &status);
//
        /// FRONT-BACK (zy)
        for (int p = 0; p < proc.vect_size; p+=proc.nX) {
            proc_vect[p] = proc_vect[p+1];
            proc_vect[p + proc.nX - 1] = proc_vect[p + proc.nX - 2];
        }

        /// AWESOME CALCULATING STAGE

        multiplicateMatrVector(&spMatrix, proc_vect, proc_vect_next, (int)proc.vect_size);
        tmp = proc_vect;
        proc_vect = proc_vect_next;
        proc_vect_next = tmp;
    }

    ////// GET THIS SHIT //////
    double* get_vect = (double *)calloc(sizeP*proc.vect_size, sizeof(double));

    double* vect_out = (double *)calloc(task.fullVectSize, sizeof(double));
    MPI_Gather(proc_vect, (int)proc.vect_size, MPI_DOUBLE,
               get_vect, (int)proc.vect_size, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

    if (rankP == ROOT) {
        for (int p = 0; p < sizeP; ++p) {

            // We searching offset for original block
            // That's means, to the y offset we should get original
            int y_block_size = task.sizeY * (task.nY / lineP); // half xy (one proc layer)
            int y_offset = (p % lineP) * y_block_size;
            int z_block_size = task.sizeZ * (task.nZ / rowP);
            int z_offset = (p / lineP) * z_block_size;

            for (int z = 1; z < proc.nZ - 1; ++z) {
                for (int y = 1; y < proc.nY - 1; ++y) {
                    for (int x = 0; x < proc.nX; ++x) {
                        vect_out[y_offset + z_offset + x + ((y - 1) * task.sizeY) +
                                ((z - 1) * task.sizeZ)] =
                                get_vect[p * proc.vect_size + x + (y * proc.sizeY) + (z * proc.sizeZ)];
                    }
                }
            }
        }
    }

    printf("finish\n");

    if (rankP == ROOT) {
        FILE *outfile = fopen(out_filename, "w");

        for (int i = 0; i < task.fullVectSize; ++i) {
            if (i % (task.nX + 2) != 0 && i % (task.nY + 2) != task.nZ + 1)
                fprintf(outfile, "%2.15le\n", vect_out[i]);
        }
    }

    MPI_Finalize();
    return 0;
}