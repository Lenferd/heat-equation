//
// Created by lenferd on 10.09.17.
//

#include "SparseMatrix.h"

int boundaries_conditions(double *vect, Task *task) {
    int realSizeX = task->nX + 2;
    int realSizeY = realSizeX;
    int realSizeZ = realSizeY * task->nZ;


    int sectionStart = 0;

//    for (int p = 0; p < proc.vect_size; p+=proc.sizeY) {
//        proc_vect[p] = proc_vect[p+1];
//        proc_vect[p + proc.sizeY - 1] = proc_vect[p + proc.sizeY - 2];
//    }

    #pragma omp parallel for
    for (int z = 0; z < task->nZ; ++z) {
        for (int y = 0; y < task->nY ; ++y) {
            // get index first and last x (boundaries)
            sectionStart = z * realSizeZ + y * realSizeY;
            vect[sectionStart] = vect[sectionStart+1];
            vect[sectionStart + realSizeX - 1] = vect[sectionStart + realSizeX -2];
        }
    }

    return 0;
}

void spMatrixInit(SparseMatrix *sp, int size, int rows, int threads) {
    sp->_size= size;
    sp->_rows = rows;
    sp->values = (double*)malloc(size*sizeof(double));
    sp->columns = (int*)malloc(size*sizeof(int));
    sp->pointerB = (int*)malloc((rows+1)*sizeof(int));

    sp->threads = threads;
}

int fillMatrix_3d_MPI_woBound(SparseMatrix *sp, MatrixVal *m_value, ProcSettings *proc){
//    int realSizeX = sizeX;
//    int realSizeY = realSizeX;
//    int realSizeZ = realSizeY * sizeY;
    int index = 0;
    int pIndex = 0;


    int sectionStart = 0;
    for (int z = 0; z < proc->nZ; ++z) {
        for (int y = 0; y < proc->nY ; ++y) {
            sectionStart = z * proc->sizeZ + y * proc->sizeY;
            for (int x = 0; x < proc->nX; ++x) {
                // Z first
                sp->values[index] = m_value->z1;
                sp->columns[index] = z == 0 ?
                                    x + sectionStart + proc->sizeZ * (proc->nZ - 1) :
                                    x + sectionStart - proc->sizeZ;
                sp->pointerB[pIndex++] = index;
                ++index;


                // Y first
                sp->values[index] = m_value->y1;
                sp->columns[index] = y == 0 ?
                                    x + sectionStart + proc->sizeY * (proc->nY - 1) :
                                    x + sectionStart - proc->sizeY;
                ++index;

                // X Group center
                sp->values[index] = m_value->x1;
                sp->columns[index] = sectionStart + x - 1;
                ++index;

                sp->values[index] = m_value->x2Comp;
                sp->columns[index] = sectionStart + x;
                ++index;

                sp->values[index] = m_value->x1;
                sp->columns[index] = sectionStart + x + 1;
                ++index;

                // Y second
                sp->values[index] = m_value->y1;
                sp->columns[index] = y == proc->nY - 1?
                                    x + sectionStart - proc->sizeY * (proc->nY - 1) :
                                    x + sectionStart + proc->sizeY;
                ++index;

                // Z second
                sp->values[index] = m_value->z1;
                sp->columns[index] = z == proc->nZ - 1 ?
                                    x + sectionStart - proc->sizeZ * (proc->nZ - 1) :
                                    x + sectionStart + proc->sizeZ;
                ++index;
            }
        }
    }

    sp->pointerB[pIndex] = index + 1;   //end
}

void multiplicateMatrVector(SparseMatrix *sp, double *input, double *result, int size){
//    omp_set_num_threads(sp->threads);

    #pragma omp parallel for
    for (int i = 0; i < size; i++){  // iteration FOR RESULT VECTOR!!!
        double local_result = 0;
        for (int j = sp->pointerB[i]; j < sp->pointerB[i+1]; j++) {
            local_result += sp->values[j] * input[sp->columns[j]];
        }
        result[i] = local_result;
    }
}