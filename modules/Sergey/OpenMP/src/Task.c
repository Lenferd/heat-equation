//
// Created by lenferd on 08.09.17.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Task.h"

int read_settings(char *settings_path, Task *task) {
    FILE *inSettingfile = fopen(settings_path, "r");

    if (inSettingfile == NULL) {
        printf("File reading error. Try to relocate input file\n");
        return -1;
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
    scaned_values += fscanf(inSettingfile, "XSTART=%lf\n", &task->xStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "XEND=%lf\n", &task->xEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "YSTART=%lf\n", &task->yStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "YEND=%lf\n", &task->yEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "ZSTART=%lf\n", &task->zStart);    // start coordinate
    scaned_values += fscanf(inSettingfile, "ZEND=%lf\n", &task->zEnd);        // end coordinate

    scaned_values += fscanf(inSettingfile, "SIGMA=%lf\n", &task->sigma);      // coef of heat conduction

    scaned_values += fscanf(inSettingfile, "NX=%d\n", &task->nX);             // count of initial elements
    scaned_values += fscanf(inSettingfile, "NY=%d\n", &task->nY);             //
    scaned_values += fscanf(inSettingfile, "NZ=%d\n", &task->nZ);             //

    scaned_values += fscanf(inSettingfile, "TSTART=%lf\n", &task->tStart);    // start time
    scaned_values += fscanf(inSettingfile, "TFINISH=%lf\n", &task->tFinish);   // finish time
    scaned_values += fscanf(inSettingfile, "dt=%lf\n", &task->dt);            // delta of time difference
    scaned_values += fscanf(inSettingfile, "BC=%d\n", &task->bc);         // Not using right now

    fclose(inSettingfile);
    if (scaned_values != 14) {
        printf("values scanned %d, must be 14\n", scaned_values);
        printf("File data reading error\n");
        return -1;
    }

    // Set timestep
    task->timeStepX = (fabs(task->xStart) + fabs(task->xEnd)) / task->nX;
    task->timeStepY = (fabs(task->yStart) + fabs(task->yEnd)) / task->nY;
    task->timeStepZ = (fabs(task->zStart) + fabs(task->zEnd)) / task->nZ;
    return 0;
}

int initFunctionData_forAdditionalXYZ(char *filename, double *vect, Task *task) {
    FILE *in_function_file = fopen(filename, "r");

    int scan_value = 0;
    /// Read file
    for (int z = 1; z < task->nZ + 1; z++) {
        for (int y = 1; y < task->nY + 1; ++y) {
            for (int i = 1; i < task->nX + 1; ++i) {
                scan_value += fscanf(in_function_file, "%lf\n",
                                     &vect[i + (task->nX + 2) * y + (task->nX + 2) * (task->nY + 2) * z]);
            }
        }
    }

    fclose(in_function_file);
    if (scan_value != task->nX * task->nY * task->nZ) {
        printf("Data reading error\n");
        exit(-3);
    }
    return 0;
}