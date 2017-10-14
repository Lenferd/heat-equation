//
// Created by lenferd on 27.03.17.
//

#ifndef HEAT_EQUATION_TASK_H
#define HEAT_EQUATION_TASK_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
using std::string;

struct Task_cpp {
    double  xStart, xEnd;   // range

    double  yStart, yEnd;
    double  zStart, zEnd;

    double  sigma;          //
    int     bc;             // not used

    double  stepX;       // time time between calculating
    double  stepY;
    double  stepZ;

    int     nX;             // count of initial elements
    int     nY;
    int     nZ;
    int     fullVectSize;

    double  tStart, tFinish;
    double  dt;
};

int initTaskUsingFile(Task_cpp &task, string settingFile);
int initMemoryReadData(double **& vect, string file, Task_cpp &task);
int initMemoryReadData_for_additional_xyz(double **& vect, string file, Task_cpp &task);
void setTimestep(Task_cpp &task);
#endif //HEAT_EQUATION_TASK_H
