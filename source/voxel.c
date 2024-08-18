#include "voxel.h"
#include "common.h"
#include <math.h>
#include <assert.h>

void environmentParametersInit(int pixelDim,
                                int angularTrajectory,
                                int positionsAngularDistance,
                                int objectSideLenght,
                                int detectorSideLength,
                                int distanceDetectorObject,
                                int distanceObjectSource,
                                int voxelXDim,
                                int voxelYDim,
                                int voxelZDim,
                                int nVoxel[3],
                                int nPlanes[3]){
    
    gl_pixelDim = pixelDim;
    gl_angularTrajectory = angularTrajectory;
    gl_positionsAngularDistance = positionsAngularDistance;
    gl_objectSideLenght = objectSideLenght;
    gl_detectorSideLength = detectorSideLength;
    gl_distanceObjectDetector = distanceDetectorObject;
    gl_distanceObjectSource = distanceObjectSource;
    gl_voxelXDim = voxelXDim;
    gl_voxelYDim = voxelYDim;
    gl_voxelZDim = voxelZDim;
    gl_nVoxel[0] = nVoxel[0];
    gl_nVoxel[1] = nVoxel[1];
    gl_nVoxel[2] = nVoxel[2];
    gl_nPlanes[0] = nPlanes[0];
    gl_nPlanes[1] = nPlanes[1];
    gl_nPlanes[2] = nPlanes[2];
}


void initTables( double *sinTable, double *cosTable, int length){
    
    const int nTheta = (int)(gl_angularTrajectory / gl_positionsAngularDistance);                      //number of angular position
    assert(nTheta < length/sizeof(sinTable[0]));

    //iterates over each source  Ntheta
    for(int positionIndex = 0; positionIndex <= nTheta; positionIndex++){
        sinTable[positionIndex] = sin((gl_angularTrajectory / 2 - positionIndex * gl_positionsAngularDistance) * M_PI / 180);
        cosTable[positionIndex] = cos((gl_angularTrajectory / 2 - positionIndex * gl_positionsAngularDistance) * M_PI / 180);
    }
}



/**
 * Generates a sub-section of a solid cubic object given its side length.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'sideLength' length of the side of the cubic object
*/
void generateCubeSlice(double *f, int nOfSlices, int offset, int sideLength){
    const int innerToOuterDiff = gl_nVoxel[X] / 2 - sideLength / 2;
    const int rightSide = innerToOuterDiff + sideLength;

#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, gl_nVoxel, offset, sideLength, innerToOuterDiff, rightSide)
    for(int n = 0 ; n < nOfSlices; n++){
        for(int i = 0; i < gl_nVoxel[Z]; i++){
            for(int j = 0; j < gl_nVoxel[X]; j++){
                f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 0;
                if((i >= innerToOuterDiff) && 
                    (i <= rightSide) && 
                    (j >= innerToOuterDiff) && 
                    (j <= rightSide) && 
                    (n + offset >= innerToOuterDiff) && 
                    (n + offset <=  gl_nVoxel[Y] - innerToOuterDiff)
                ){
                    f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 1.0;
                } else {
                    f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 0.0;
                }
            }
        }
    }
}

/**
 * Generates a sub-section of a solid spherical object given its diameter.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'diameter' is the diameter of the sphere.
*/
void generateSphereSlice(double *f, int nOfSlices, int offset, int diameter)
{
#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, gl_nVoxel, offset, diameter, gl_objectSideLenght, gl_voxelYDim, gl_voxelXDim, gl_voxelZDim)
    for (int n = 0; n < nOfSlices; n++) {
        for (int r = 0; r < gl_nVoxel[Z]; r++) {
            for (int c = 0; c < gl_nVoxel[X]; c++) {
                Point temp;
                temp.y = -(gl_objectSideLenght / 2) + (gl_voxelYDim / 2) + (n + offset) * gl_voxelYDim;
                temp.x = -(gl_objectSideLenght / 2) + (gl_voxelXDim / 2) + (c) * gl_voxelXDim;
                temp.z = -(gl_objectSideLenght / 2) + (gl_voxelZDim / 2) + (r) * gl_voxelZDim;
                const double distance = sqrt(pow(temp.x, 2) + pow(temp.y, 2) + pow(temp.z, 2));
                if(distance <= diameter && c < gl_nVoxel[Z] / 2){
                    f[(gl_nVoxel[Z]) * r + c + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 1;
                } else {
                    f[(gl_nVoxel[Z]) * r + c + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 0.0;
                }
            }
        }
    }
}

/**
 * Generates a sub-section of a solid cubic object with an internal spherical cavity.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'sideLength' lenght of the side of the object.
*/
void generateCubeWithSphereSlice(double *f, int nOfSlices, int offset, const int sideLength){
    const int innerToOuterDiff = gl_nVoxel[X] / 2 - sideLength / 2;
    const int rightSide = innerToOuterDiff + sideLength;
    const Point sphereCenter = {-sideLength * gl_voxelXDim / 4, -sideLength * gl_voxelYDim / 4, -sideLength * gl_voxelZDim / 4};

    // iterates over each voxel of the grid 
#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, offset, sideLength, gl_nVoxel, innerToOuterDiff, rightSide, sphereCenter, gl_objectSideLenght, gl_voxelYDim, gl_voxelXDim, gl_voxelZDim)
    for(int n = 0 ; n < nOfSlices; n++){
        for(int i = 0; i < gl_nVoxel[Z]; i++){
            for(int j = 0; j < gl_nVoxel[X]; j++){

                f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 0;
                if ( (i >= innerToOuterDiff) &&
                     (i <= rightSide) &&
                     (j >= innerToOuterDiff) &&
                     (j <= rightSide) &&
                     (n + offset >= innerToOuterDiff) &&
                     (n + offset <=  gl_nVoxel[Y] - innerToOuterDiff) ) {
                    
                    // voxel position is inside the cubic object
                    Point temp;
                    temp.y = -(gl_objectSideLenght / 2) + (gl_voxelYDim / 2) + (n + offset) * gl_voxelYDim;
                    temp.x = -(gl_objectSideLenght / 2) + (gl_voxelXDim / 2) + (j) * gl_voxelXDim;
                    temp.z = -(gl_objectSideLenght / 2) + (gl_voxelZDim / 2) + (i) * gl_voxelZDim;
                    const double distance = sqrt(pow(temp.x - sphereCenter.x, 2) + pow(temp.y - sphereCenter.y, 2) + pow(temp.z - sphereCenter.z, 2));
                    
                    if(distance > sideLength * gl_voxelXDim / 6){
                        // voxel position is inside the cubic object and inside the sperical cavity
                        f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 1.0;
                    }
                } else {
                    // voxel position is outside the cubic object
                    f[(gl_nVoxel[Z]) * i + j + n * gl_nVoxel[X] * gl_nVoxel[Z]] = 0.0;
                }
            }
        }
    }
}