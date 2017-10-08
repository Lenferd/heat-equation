//
// Created by lenferd on 27.03.17.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "Task.h"

using std::string;

int initTaskUsingFile(Task &task, string settingFile) {
    FILE *inSettingfile = fopen(settingFile.c_str(), "r");
    if (inSettingfile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        exit(-1);
    }

//    XSTART=-1.0
//    XEND=1.0
//    YSTART=-1.0
//    YEND=1.0
//    ZSTART=-1.0
//    ZEND=1.0
//    SIGMA=1.0
//    NX=50
//    NY=50
//    NZ=50
//    TSTART=0.000
//    TFINISH=8e-3
//    dt=8e-8
//    BC=2

    // File reading
    int scaned_values = 0;
    scaned_values += fscanf(inSettingfile, "XSTART=%lf\n", &task.xStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "XEND=%lf\n", &task.xEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "YSTART=%lf\n", &task.yStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "YEND=%lf\n", &task.yEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "ZSTART=%lf\n", &task.zStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "ZEND=%lf\n", &task.zEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "SIGMA=%lf\n", &task.sigma);      // coef of heat conduction

    scaned_values += fscanf(inSettingfile, "NX=%d\n", &task.nX);             // count of initial elements
    scaned_values += fscanf(inSettingfile, "NY=%d\n", &task.nY);             //
    scaned_values += fscanf(inSettingfile, "NZ=%d\n", &task.nZ);             //

    scaned_values += fscanf(inSettingfile, "TSTART=%lf\n", &task.tStart);    // start time
    scaned_values += fscanf(inSettingfile, "TFINISH=%lf\n", &task.tFinish);   // finish time
    scaned_values += fscanf(inSettingfile, "dt=%lf\n", &task.dt);            // delta of time difference
    scaned_values += fscanf(inSettingfile, "BC=%d\n", &task.bc);         // Not using right now

    fclose(inSettingfile);
    if (scaned_values != 14) {
        printf("values scanned %d, must be 14\n", scaned_values);
        printf("File data reading error\n");
        exit(-2);
    }
    return 0;
}

void setTimestep(Task &task){
    task.stepX = (fabs(task.xStart) + fabs(task.xEnd)) / task.nX;

    task.stepY = (fabs(task.yStart) + fabs(task.yEnd)) / task.nY;
    task.stepZ = (fabs(task.zStart) + fabs(task.zEnd)) / task.nZ;
}

int initMemoryReadData(double **& vect, string file, Task &task) {
    FILE *inFunctionfile = fopen(file.c_str(), "r");

    vect = new double*[2];
    task.fullVectSize = (task.nX + 2) * (task.nY) * (task.nZ);
    vect[0] = new double[task.fullVectSize];
    vect[1] = new double[task.fullVectSize];

    int scan_value = 0;
    /// Read file
    for (int k = 0; k < task.nZ; k++) {
        for (int j = 0; j < task.nY; ++j) {
            for (int i = 1; i < task.nX + 1; ++i) {
                scan_value += fscanf(inFunctionfile, "%lf\n",
                                     &vect[0][i + (task.nX + 2) * j + (task.nX + 2) * task.nY * k]);
            }
        }
    }

    fclose(inFunctionfile);
    if (scan_value != task.nX * task.nY * task.nZ) {
        printf("Data reading error\n");
        exit(-3);
    }
    return 0;
}

int initMemoryReadDataMPI(double *& vect, string file, Task &task) {
    FILE *inFunctionfile = fopen(file.c_str(), "r");

    task.fullVectSize = (task.nX + 2) * (task.nY) * (task.nZ);
    vect = new double[task.fullVectSize];

    int error;
    /// Read file
    for (int k = 0; k < task.nZ; k++) {
        for (int j = 0; j < task.nY; ++j) {
            for (int i = 1; i < task.nX + 1; ++i) {
                error += fscanf(inFunctionfile, "%lf\n", &vect[i + (task.nX + 2) * j + (task.nX+2) * task.nY * k]);
            }
        }
    }

    fclose(inFunctionfile);

    return 0;
}

int read_settings(char *settings_path, Task *task) {
    return 0;
}
