/****************************************************************************
 *
 * projector.c
 *
 * Copyright (C) 2024 Lorenzo Colletta
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************************/


/***
This program implements the Siddon algorithm to compute 2-dimensional projections of a three dimensional voxel grid.

COMPILE:
This program can be compiled so that: it returns the output in a text file in pgm format, in this case the calculated
values ​​are converted to a gray scale from 0 to 255; returns the output in a binary file, in this case the values ​​are
stored as double values ​​for maximum accuracy, the structure of the file is described below.

gcc command:

    gcc -Wall -Wpedantic -std=c99 -fopenmp projector.c ./source/voxel.c -I./source/ -o projector

Compilation to use a binary file as output:

    gcc -Wall -Wpedantic -std=c99 -fopenmp -DBINARY projector.c ./source/voxel.c -I./source/ -o projector

RUN:

    projector input.dat output.dat/output.pgm

- First parameter is the name of the input file;
- Second parameter is the name of a text or a binary file to store the output in.

BINARY FILE STRUCTURE:

The computed output is a set of projection images each with a resolution of 'n'*'n'.
The output file is structured as follows:
- the first value is the number of images produced (integer type)
- the second value is 'n' (integer type), the resolution of the image side
- the third value is the maximum value computed (type double)
- the third value is the minimum value computed (type double)
- each sequence of values ​​representing an image is preceded by a value (type double) which indicates the angle from
  which the image was computed
- an image is a sequence of n*n double values; the image is stored as a one-dimensional array, sorted first by the x
  coordinate and then by the z coordinate (considering a three-dimensional Cartesian space with the x axis from left
  to right, the y axis oriented upwards and z perpendicular to them)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "projection.h"
#include "common.h"

#ifdef _OPENMP
#include <omp.h>
#else
double omp_get_wtime( void ) { return 0; }
#endif

#define OBJ_BUFFER 100                          // limits the number of voxel alogn the y axis computed per time
#define TABLES_DIM 1024

/**
 * The following global variables are defined as according to common.h header file.
 */
int gl_pixelDim;
int gl_angularTrajectory;
int gl_positionsAngularDistance;
int gl_objectSideLenght;
int gl_detectorSideLength;
int gl_distanceObjectDetector;
int gl_distanceObjectSource;
int gl_voxelXDim;
int gl_voxelYDim;
int gl_voxelZDim;
int gl_nVoxel[3];
int gl_nPlanes[3];

/**
 * The following global variables are defined as according to projection.h header file.
 */
double *sineTable;
double *cosineTable;

/**
 * Reeds the environment values used to compute the voxel grid from the specified binary file.
 * 'filePointer' the file pointer to read the values from.
 * @returns 0 in case of writing failure, 1 otherwise.
 */
int readSetUP(FILE* filePointer)
{
    int buffer[16];
    //read object parameters set up
    if(!fread(buffer, sizeof(int), 16, filePointer)){
        return 0;
    }

    gl_pixelDim = buffer[0];
    gl_angularTrajectory = buffer[1];
    gl_positionsAngularDistance = buffer[2];
    gl_objectSideLenght = buffer[3];
    gl_detectorSideLength = buffer[4];
    gl_distanceObjectDetector = buffer[5];
    gl_distanceObjectSource = buffer[6];
    gl_voxelXDim = buffer[7];
    gl_voxelYDim = buffer[8];
    gl_voxelZDim = buffer[9];
    gl_nVoxel[0] = buffer[10];
    gl_nVoxel[1] = buffer[11];
    gl_nVoxel[2] = buffer[12];
    gl_nPlanes[0] = buffer[13];
    gl_nPlanes[1] = buffer[14];
    gl_nPlanes[2] = buffer[15];

    return 1;
}

int main(int argc, char *argv[])
{
    FILE* inputFilePointer;
    FILE* outputFilePointer;

    if(argc != 3){
#ifdef BINARY
        fprintf(stderr,"Usage: %s input.dat output.dat\nFirst parameter is the name of the input file;\nSecond parameter is the name of a binary file to store the output at.\n",argv[0]);
#else
        fprintf(stderr,"Usage: %s input.dat output.pgm\nFirst parameter is the name of the input file;\nSecond parameter is the name of a text file to store the output at.\n",argv[0]);
#endif
        return EXIT_FAILURE;
    }
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];


    inputFilePointer = fopen(inputFileName,"rb");
    if (!inputFilePointer){
        printf("Unable to open file %s!\n",inputFileName);
        exit(1);
    }

    if(!readSetUP(inputFilePointer)){
        printf("Unable to read from file!\n");
        exit(1);
    }

    int nSidePixels = gl_detectorSideLength / gl_pixelDim;

    //number of angular positions
    const int nTheta = (int)(gl_angularTrajectory / gl_positionsAngularDistance);
    //array containing the coefficents of each voxel
    double *grid = (double*)malloc(sizeof(double) * gl_nVoxel[X] * gl_nVoxel[Z] * OBJ_BUFFER);
    //array containing the computed absorption detected in each pixel of the detector
    double *absorbment = (double*)calloc(nSidePixels * nSidePixels * (nTheta + 1), sizeof(double));
    //each thread has its own variable to store its minimum and maximum absorption computed
    double absMaxValue, absMinValue;

    sineTable = (double*)malloc(sizeof(double)*TABLES_DIM);
    cosineTable = (double*)malloc(sizeof(double)*TABLES_DIM);

    initTables(sineTable, cosineTable, TABLES_DIM);

    double totalTime = omp_get_wtime();

#ifdef BINARY
    outputFilePointer = fopen(outputFileName,"wb");
    if (!outputFileName){
        printf("Unable to open file!\n");
        exit(1);
    }
#else
    outputFilePointer = fopen(outputFileName,"w");
    if (!outputFileName){
        printf("Unable to open file %s!\n",outputFileName);
        exit(1);
    }
#endif

    //iterates over object subsection
    for(int slice = 0; slice < gl_nVoxel[Y]; slice += OBJ_BUFFER){

        //read voxels coefficients
        if(slice < gl_nVoxel[Y]){
            if(!fread(grid, sizeof(double), gl_nVoxel[X] * gl_nVoxel[Z] * OBJ_BUFFER, inputFilePointer)){
                printf("Unable to read on file!\n");
                exit(1);
            }
        } else {
            if(!fread(grid, sizeof(double), gl_nVoxel[X] * gl_nVoxel[Z] * (slice - gl_nVoxel[Y]), inputFilePointer)){
                printf("Unable to write on file!\n");
                exit(1);
            }
        }

        //computes subsection projection
        computeProjections(slice, grid, absorbment, &absMaxValue, &absMinValue);

    }
    fprintf(stderr,"Execution time: %lf\n", omp_get_wtime() - totalTime);
    fflush(stderr);

//write on file
#ifdef BINARY
    int matrixDetails[] = {nTheta + 1, nSidePixels};
    if(!fwrite(matrixDetails, sizeof(int), 2, outputFilePointer)){
        printf("Unable to write on file!\n");
        exit(1);
    }
    double minMaxValues[] = {absMaxValue, absMinValue};
    if(!fwrite(minMaxValues, sizeof(double), 2, outputFilePointer)){
        printf("Unable to write on file!\n");
        exit(1);
    }
    for(int i = 0; i <= nTheta; i++){
        double angle = -gl_angularTrajectory / 2 + i * gl_positionsAngularDistance;
        if(!fwrite(&angle, sizeof(double), 1, outputFilePointer)){
            printf("Unable to write on file!\n");
            exit(1);
        }
        if(!fwrite(absorbment + i * nSidePixels * nSidePixels, sizeof(double), nSidePixels * nSidePixels, outputFilePointer)){
            printf("Unable to write on file!\n");
            exit(1);
        }
    }
#else

    //iterates over each absorption value computed, prints a value between [0-255]
    fprintf(outputFilePointer,"P2\n%d %d\n255", nSidePixels, nSidePixels * (nTheta + 1));
    for(double positionIndex = 0; positionIndex <= nTheta; positionIndex ++){
        double angle = - gl_angularTrajectory / 2 + positionIndex * gl_positionsAngularDistance;
        fprintf(outputFilePointer,"\n#%lf",angle);
        for(int i = 0; i < nSidePixels; i++ ){
            fprintf(outputFilePointer,"\n");
            for(int j = 0; j < nSidePixels; j++ ){
                int pixelIndex = positionIndex * nSidePixels * nSidePixels + i *nSidePixels + j;
                int color =  (absorbment[pixelIndex] - absMinValue) * 255 / (absMaxValue - absMinValue);
                fprintf(outputFilePointer,"%d ", color);
            }
        }
    }

#endif


    free(grid);
    free(absorbment);

    fclose(inputFilePointer);
    fclose(outputFilePointer);
}
