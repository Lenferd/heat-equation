#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sp_mat.h>

#include "sp_mat.h"

void initSpMat(SpMatrix *mat, size_t nz, size_t nRows) {
  mat->nz = nz;
  mat->nRows = nRows;
  mat->value = (TYPE *)malloc(sizeof(TYPE) * nz);
  mat->col = (int *)malloc(sizeof(int) * nz);
  mat->rowIndex = (int *)malloc(sizeof(int) * (nRows) + 1);
  memset(mat->rowIndex, 0, nRows + 1);
}

void freeSpMat(SpMatrix* mat) {
  free(mat->value);
  free(mat->col);
  free(mat->rowIndex);
}

void multMV(TYPE** result, SpMatrix mat, TYPE* vec) {
  TYPE localSum;
  #pragma omp parallel private(localSum) num_threads(2) if (ENABLE_PARALLEL)
  {
    #pragma omp for nowait
    for (int i = 0; i < mat.nRows; i++) {
      localSum = 0.0;
      for (int j = mat.rowIndex[i]; j < mat.rowIndex[i + 1]; j++)
        localSum += mat.value[j] * vec[mat.col[j]];
      (*result)[i] = localSum;
    }
  }
}

void sumV(TYPE **result, TYPE *U, TYPE *k1, TYPE *k2, TYPE *k3, TYPE *k4, size_t N, double h) {
  #pragma omp parallel for num_threads(2) if (ENABLE_PARALLEL)
  for (int i = 0; i < N; i++)
    (*result)[i] = U[i] + h*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}


void printSpMat(SpMatrix mat) {
  for (int i = 0; i < mat.nRows; i++) {
    for (int j = 0; j < mat.nRows; j++)
      printf("%.0lf", procedure(mat, i, j));
    printf("\n");
  }
}

TYPE procedure(SpMatrix mat, int i, int j) {
  TYPE result = 0;
  int N1 = mat.rowIndex[i];
  int N2 = mat.rowIndex[i+1];
  for(int k = N1; k < N2; k++) {
    if (mat.col[k] == j) {
      result = mat.value[k];
      break;
    }
  }
  return result;
}


void denseMult(double **result, double **mat, double *vec, size_t dim) {
  memset(*result, 0, dim*sizeof(double));
  for (int x = 0; x < dim; x++) {
    for (int i = 0;i < dim;i++)
      (*result)[x]+=mat[x][i]*vec[i];
  }
}
void createExplicitSpMat(SpMatrix *mat, TYPE coeffs[5], int dim, int NX, int NXY) {
  int index = 0, j, k;
  mat->rowIndex[0] = 0;
  for (int i = 0; i < dim; i++) {
    if (i % NX != 0 && i % NX != NX - 1) {
      // Смещение на x + 1 и x - 1 с учётом граничных условий
      // ***************************************
      mat->col[index] = i - 1;
      mat->value[index] = coeffs[0];
//      mat->value[index] = 1;
      index++;

      mat->col[index] = i;
      mat->value[index] = coeffs[1];
//      mat->value[index] = 2;
      index++;

      mat->col[index] = i + 1;
      mat->value[index] = coeffs[0];
//      mat->value[index] = 1;
      index++;
      // ***************************************

      // Смещение на y + 1 и y - 1 с учётом цикличности условия
      // ***************************************
      j = i - NX;
      if (j <= 0) {
        mat->col[index] = dim + j;
        mat->value[index] = coeffs[2];
//        mat->value[index] = 3;

        index++;
      } else {
        mat->col[index] = j;
        mat->value[index] = coeffs[2];
//        mat->value[index] = 3;

        index++;
      }
      mat->col[index] = (i + NX) % NXY;
      mat->value[index] = coeffs[2];
//      mat->value[index] = 3;

      index++;
      // ***************************************

      // Смещение на z + 1 и z - 1 с учётом цикличности условия
      // ***************************************
      k = i - NXY;
      if (k <= 0) {
        mat->col[index] = (int) dim + k;
        mat->value[index] = coeffs[3];
        index++;
      } else {
        mat->col[index] = k;
        mat->value[index] = coeffs[3];
        index++;

      }
      mat->col[index] = (i + NXY) % (int) dim;
      mat->value[index] = coeffs[3];
      index++;
      // ***************************************
      
      mat->rowIndex[i + 1] = mat->rowIndex[i] + 7;
    } else if (i % NX == 0) {
      // ГРАНИЧНЫЕ УСЛОВИЯ СЛЕВА
      mat->col[index] = i;
      mat->value[index] = coeffs[0];
//      mat->value[index] = 1;
      index++;

      mat->col[index] = i + 1;
      mat->value[index] = coeffs[1];
//      mat->value[index] = 2;
      index++;

      mat->col[index] = i + 2;
//      mat->value[index] = 1;
      mat->value[index] = coeffs[0];
      index++;

      mat->rowIndex[i + 1] = mat->rowIndex[i] + 3;
    } else if (i % NX == NX - 1) {
      // ГРАНИЧНЫЕ УСЛОВИЯ СПРАВА
      mat->col[index] = i - 2;
      mat->value[index] = coeffs[0];
//      mat->value[index] = 1;
      index++;

      mat->col[index] = i - 1;
      mat->value[index] = coeffs[1];
//      mat->value[index] = 2;
      index++;

      mat->col[index] = i;
//      mat->value[index] = 1;
      mat->value[index] = coeffs[0];
      index++;

      mat->rowIndex[i + 1] = mat->rowIndex[i] + 3;
    }
  }

}
