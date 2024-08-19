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

#define OBJ_BUFFER 1000
#define TABLES_DIM 1024

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

double *sineTable;
double *cosineTable;

int readSetUP(FILE* filePointer){
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
    int header[] = {nTheta + 1, nSidePixels};
    if(!fwrite(header, sizeof(int), 2, outputFilePointer)){
        printf("Unable to write on file!\n");
        exit(1);
    }
    for(int i = 0; i < nTheta; i++){
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
