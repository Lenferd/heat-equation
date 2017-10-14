//
// Created by lenferd on 08.09.17.
//

#ifndef HEAT_EQUATION_TASK_H
#define HEAT_EQUATION_TASK_H

#include <stdio.h>
typedef struct TaskStruct {
    double  xStart, xEnd;   // range

    double  yStart, yEnd;
    double  zStart, zEnd;

    double  sigma;          //
    int     bc;             // not used

    double  timeStepX;       // time time between calculating
    double  timeStepY;
    double  timeStepZ;

    int     nX;             // count of initial elements
    int     nY;
    int     nZ;

    int     sizeY;
    int     sizeZ;
    int  fullVectSize;

    double  tStart, tFinish;
    double  dt;
} Task;


typedef struct ProcSettings {
    int nX;
    int nY;
    int nZ;

    int sizeY;
    int sizeZ;

    size_t vect_size;
} ProcSettings;

int read_settings(char *settings_path, Task *task);
int initFunctionData_forAdditionalXYZ(char *filename, double *vect, Task *task);
#endif //HEAT_EQUATION_TASK_H
