/**
 * Usage:
 *  compile:  gcc -std=c99 -Wall -Wpedantic -fopenmp  projector.c -lm -o projector
 *  run:      ./projector 0 1 > CubeWithSphere.pgm
 *  convert:  convert CubeWithSphere.pgm CubeWithSphere.jpeg
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#ifdef _OPENMP
#include <omp.h>
#else
double omp_get_wtime( void ) { return 0; }
#endif

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define VOXEL_X 100             // voxel side along x-axis
#define VOXEL_Y 100             // voxel side along y-axis
#define VOXEL_Z 100             // voxel side along z-axis

#define PIXEL 85                // detector's pixel side lenght

#define VOXEL_MAT 100000        // entire voxel matrix side lenght

#define DETECTOR 200000         // detector lenght

#define DOD 150000              //distance from object center to detector
#define DOS 600000              //distance from source to object center

#define AP 90                   //source path angle
#define STEP_ANGLE 15           //angular distance between each source step

#define OBJ_BUFFER 100          //voxel coefficients buffer size

//cartesian axis
enum axis{
    X,
    Y,
    Z
};

//models a point of coordinates (x,y,z) in the cartesian coordinate system
struct point
{
    double x;
    double y;
    double z;
};

//models a structure containing the range of indices of the planes to compute the intersection with
struct ranges{
    int xMinIndx;
    int xMaxIndx;
    int yMinIndx;
    int yMaxIndx;
    int zMinIndx;
    int zMaxIndx;
};

double sin_table[1024], cos_table[1024];

// number of voxel along X axis [0], Y axis [1], Z Axis [2]
const int nVoxel[3] = {VOXEL_MAT / VOXEL_X, VOXEL_MAT / VOXEL_Y, VOXEL_MAT / VOXEL_Z};

//number of parallel planes alogn X axis [0], Y axis [1], Z axis [2]
const int nPlanes[3] = {(VOXEL_MAT / VOXEL_X) + 1, (VOXEL_MAT / VOXEL_Y) + 1, (VOXEL_MAT / VOXEL_Z) + 1};

//position of the top-left pixel of the detector relative to the center of the detector
const double elementOffset = DETECTOR / 2 - PIXEL / 2;

//number of pixel along a detector's side
const int nSidePixels = DETECTOR / PIXEL;

//flag indicating whether the detector rotates so that it's always orthogonal to the
//line passing between the source and the center of the detector
int stationaryDetector = 0;

void init_tables( void )
{
    const int nTheta = (int)(AP / STEP_ANGLE);                      //number of angular position
    assert(nTheta < sizeof(sin_table)/sizeof(sin_table[0]));
    //iterates over each source  Ntheta
    for(int angle = 0; angle <= nTheta; angle++){
        sin_table[angle] = sin((AP / 2 - angle * STEP_ANGLE) * M_PI / 180);
        cos_table[angle] = cos((AP / 2 - angle * STEP_ANGLE) * M_PI / 180);
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
    const int innerToOuterDiff = nVoxel[X] / 2 - sideLength / 2;
    const int rightSide = innerToOuterDiff + sideLength;

#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, nVoxel, offset, sideLength, innerToOuterDiff, rightSide)
    for(int n = 0 ; n < nOfSlices; n++){
        for(int i = 0; i < nVoxel[Z]; i++){
            for(int j = 0; j < nVoxel[X]; j++){
                f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 0;
                if( (i >= innerToOuterDiff) && (i <= rightSide) && (j >= innerToOuterDiff) && (j <= rightSide) && (n + offset >= innerToOuterDiff) && (n + offset <=  nVoxel[Y] - innerToOuterDiff)){
                    f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 1.0;
                } else {
                    f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 0.0;
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
#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, nVoxel, offset, diameter)
    for (int n = 0; n < nOfSlices; n++) {
        for (int r = 0; r < nVoxel[Z]; r++) {
            for (int c = 0; c < nVoxel[X]; c++) {
                struct point temp;
                temp.y = -(VOXEL_MAT / 2) + (VOXEL_Y / 2) + (n + offset) * VOXEL_Y;
                temp.x = -(VOXEL_MAT / 2) + (VOXEL_X / 2) + (c) * VOXEL_X;
                temp.z = -(VOXEL_MAT / 2) + (VOXEL_Z / 2) + (r) * VOXEL_Z;
                const double distance = sqrt(pow(temp.x, 2) + pow(temp.y, 2) + pow(temp.z, 2));
                if(distance <= diameter && c < nVoxel[Z] / 2){
                    f[(nVoxel[Z]) * r + c + n * nVoxel[X] * nVoxel[Z]] = 1;
                } else {
                    f[(nVoxel[Z]) * r + c + n * nVoxel[X] * nVoxel[Z]] = 0.0;
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
void generateCubeWithSphereSlice(double *f, int nOfSlices, int offset, int sideLength){
    const int innerToOuterDiff = nVoxel[X] / 2 - sideLength / 2;
    const int rightSide = innerToOuterDiff + sideLength;
    const struct point sphereCenter = {-15000, -15000, 1500};

#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, offset, sideLength, nVoxel, innerToOuterDiff, rightSide, sphereCenter)
    for(int n = 0 ; n < nOfSlices; n++){
        for(int i = 0; i < nVoxel[Z]; i++){
            for(int j = 0; j < nVoxel[X]; j++){
                f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 0;
                if ( (i >= innerToOuterDiff) &&
                     (i <= rightSide) &&
                     (j >= innerToOuterDiff) &&
                     (j <= rightSide) &&
                     (n + offset >= innerToOuterDiff) &&
                     (n + offset <=  nVoxel[Y] - innerToOuterDiff) ) {
                    struct point temp;
                    temp.y = -(VOXEL_MAT / 2) + (VOXEL_Y / 2) + (n + offset) * VOXEL_Y;
                    temp.x = -(VOXEL_MAT / 2) + (VOXEL_X / 2) + (j) * VOXEL_X;
                    temp.z = -(VOXEL_MAT / 2) + (VOXEL_Z / 2) + (i) * VOXEL_Z;
                    const double distance = sqrt(pow(temp.x - sphereCenter.x, 2) + pow(temp.y - sphereCenter.y, 2) + pow(temp.z - sphereCenter.z, 2));
                    if(distance > 10000)
                        f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 1.0;
                } else {
                    f[(nVoxel[Z]) * i + j + n * nVoxel[X] * nVoxel[Z]] = 0.0;
                }
            }
        }
    }
}

/**
 * returns the coordinate of a plane parallel to the YZ plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getXPlane(int index){
    return -((VOXEL_MAT) / 2) + index * VOXEL_X;
}

/**
 * returns the coordinate of a plane parallel to the XZ plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getYPlane(int index){
    return -((VOXEL_MAT) / 2) + index * VOXEL_Y;
}

/**
 * returns the coordinate of a plane parallel to the XY plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getZPlane(int index){
    return -((VOXEL_MAT) / 2) + index * VOXEL_Z;
}

/**
 * Computes the parametrical value of the intersection between a ray and a plane given the coordinate of the plane.
 * Returns the axis to which the ray is orthogonal, -1 otherwise.
 * 'a' and 'b' are the points that identify the ray.
 * 'value' is the coordinate of the plane.
 * 'ax' is the axis orthogonal to the plane.
 * 'inters' is the parametrical value that identify the intersection point along the ray.
*/
int getIntersection(const struct point a, const struct point b, double value, enum axis ax, double *inters){
    switch (ax){
        case X:
            if(a.x - b.x != 0){
                *inters = (value - a.x) / (b.x - a.x);
                return -1;
            } else {
                return X;
            }
            break;
        case Y:
            if(a.y - b.y != 0){
                *inters = (value - a.y) / (b.y - a.y);
                return -1;
            } else {
                return Y;
            }
            break;
        case Z:
            if(a.z - b.z != 0){
                *inters = (value - a.z) / (b.z - a.z);
                return -1;
            } else {
                return Z;
            }
        default:
                return -1;
            break;
    }
    return -1;
}

/**
 * Computes the parametrical value of the intersection between a ray and a plane given the index of the plane.
 * Returns the axis to which the ray is orthogonal, -1 otherwise.
 * 'a' and 'b' are the points that identify the ray.
 * 'planeIndex' is the index of the plane.
 * 'ax' is the axis orthogonal to the plane.
 * 'inters' is the parametrical value that identify the intersection point along the ray.
*/
int getIntersectionFromIndex(const struct point a, const struct point b, int planeIndex, enum axis ax, double *inters){
    double value;
    if (ax == X){
        value = getXPlane(planeIndex);
    } else if(ax == Y){
        value = getYPlane(planeIndex);
    } else if(ax == Z){
        value = getZPlane(planeIndex);
    } else
        assert(0);
    return getIntersection(a, b, value, ax, inters);
}

/**
 * Computes the parametrical values of the intersections between a ray and the object sides along a given axis.
 * Returns the axis to which the ray is orthogonal, -1 otherwise.
 * 'a' and 'b' are the points that identify the ray.
 * 'sideA' is the intersection with the side with the smallest-valued coordinate.
 * 'sideB' is the intersection with the side with the greatest-valued coordinate.
 * 'ax' is the axis orthogonal to the plane.
 * 'slice' index of the sub-section of the object currently available.
*/
int getSidesIntersection(const struct point a, const struct point b, double *sideA, double *sideB, enum axis ax, int slice){
    double firstSide = 0;
    double lastSide = nPlanes[ax] - 1;
    if(ax == Y){
        firstSide = slice;
        lastSide = slice + OBJ_BUFFER;
    }
    getIntersectionFromIndex(a, b, firstSide, ax, sideA);
    return getIntersectionFromIndex(a, b, lastSide, ax, sideB);
}

/**
 * Returns the maximum parametric value a, representing the last intersection between ray and object
 * 'a' is the array containing the parametrical value of the intersection between the ray and the object's side along each axis.
 * 'isParallel' as a value corrisponding to the axis to which the array is orthogonal, -1 otherwise.
*/
double getAMax(double a[3][2], int isParallel){
    double tempMax[3];
    double aMax = 1;
    for(int i = 0; i < 3; i++){
        if(i != isParallel)
            tempMax[i] = a[i][0] > a[i][1] ? a[i][0] : a[i][1];
    }
    for(int i = 0; i < 3; i++){
        if(i != isParallel) /* TODO: aMax sembra il _minimo_ di una serie di valor, è ok? */
            aMax = aMax < tempMax[i] ? aMax : tempMax[i];
    }
    return aMax;
}

/**
 * Returns the minimum parametric value a, representing the fist intersection between ray and object
 * 'a' is the array containing the parametrical value of the intersection between the ray and the object's side along each axis.
 * 'isParallel' has a value corrisponding to the axis to which the array is orthogonal, -1 otherwise.
*/
double getAMin(double a[3][2], int isParallel ){
    double tempMin[3];
    double aMin = 0;
    for(int i = 0; i < 3; i++){
        if(i != isParallel)
            tempMin[i] = a[i][0] < a[i][1] ? a[i][0] : a[i][1];
    }
    for(int i = 0; i < 3; i++){
        if(i != isParallel) /* aMin sembra il _massimo_ di una serie di valori, è ok? */
            aMin = aMin > tempMin[i] ? aMin : tempMin[i];
    }
    return aMin;
}

/**
 * Returns the range of indices of the planes.
 * 'source' and 'pixel' are the points that identify the ray.
 * 'isParallel'  has a value corrisponding to the axis to which the array is orthogonal, -1 otherwise.
 * 'aMin' is the minimum parametrical value of the intersection between the ray and the object.
 * 'aMax' is the maximum parametrical value of the intersection between the ray and the object.
*/
struct ranges getRangeOfIndex(const struct point source, const struct point pixel, int isParallel, double aMin, double aMax){
    struct ranges idxs;

    //gets range of indeces of XZ parallel planes
    if(isParallel != Y){
        if(pixel.y - source.y >= 0){
            idxs.yMinIndx = nPlanes[Y] - ceil((getYPlane(nPlanes[Y] - 1) - aMin * (pixel.y - source.y) - source.y) / VOXEL_Y);
            idxs.yMaxIndx = 1 + floor((aMax * (pixel.y - source.y) + source.y - getYPlane(0)) / VOXEL_Y);
        } else {
            idxs.yMinIndx = nPlanes[Y] - ceil((getYPlane(nPlanes[Y] - 1) - aMax * (pixel.y - source.y) - source.y) / VOXEL_Y);
            idxs.yMaxIndx = floor((aMin * (pixel.y - source.y) + source.y - getYPlane(0)) / VOXEL_Y);
        }
    } else {
        idxs.yMinIndx = 0;
        idxs.yMaxIndx = 0;
    }

    //gets range of indeces of YZ parallel planes
    if(isParallel != X){
        if(pixel.x - source.x >= 0){
            idxs.xMinIndx = nPlanes[X] - ceil((getXPlane(nPlanes[X] - 1) - aMin * (pixel.x - source.x) - source.x) / VOXEL_X);
            idxs.xMaxIndx = 1 + floor((aMax * (pixel.x - source.x) + source.x - getXPlane(0)) / VOXEL_X);
        } else {
            idxs.xMinIndx = nPlanes[X] - ceil((getXPlane(nPlanes[X] - 1) - aMax * (pixel.x - source.x) - source.x) / VOXEL_X);
            idxs.xMaxIndx = floor((aMin * (pixel.x - source.x) + source.x - getXPlane(0)) / VOXEL_X);
        }
    } else {
        idxs.xMinIndx = 0;
        idxs.xMaxIndx = 0;
    }

    //gets range of indeces of XY parallel planes
    if(isParallel != Z){
        if(pixel.z - source.z >= 0){
            idxs.zMinIndx = nPlanes[Z] - ceil((getZPlane(nPlanes[Z] - 1) - aMin * (pixel.z - source.z) - source.z) / VOXEL_Z);
            idxs.zMaxIndx = 1 + floor((aMax * (pixel.z - source.z) + source.z - getZPlane(0)) / VOXEL_Z);
        } else {
            idxs.zMinIndx = nPlanes[Z] - ceil((getZPlane(nPlanes[Z] - 1) - aMax * (pixel.z - source.z) - source.z) / VOXEL_Z);
            idxs.zMaxIndx = floor((aMin * (pixel.z - source.z) + source.z - getZPlane(0)) / VOXEL_Z);
        }

    } else {
        idxs.zMinIndx = 0;
        idxs.zMaxIndx = 0;
    }

    return idxs;
}

/**
 * Computes each parametric value of the intersection between the ray and the planes whose index is in the range planeIndexRange.
 * 'source' and 'pixel' are the points that identify the ray.
 * 'planeIndexRange' is a structure containing the ranges of indeces of planes for each axis.
 * 'a' is a pointer to the array on which to store the parametrical values.
 * 'ax' is the axis orthogonal to the set of planes to which compute the intersection.
*/
void getAllIntersections(const struct point source, const struct point pixel, const struct ranges planeIndexRange, double *a, enum axis ax){
    int start = 0, end = 0;
    int direction = 1;
    double plane;
    double d;
    if(ax == X){
        start = planeIndexRange.xMinIndx;
        end = planeIndexRange.xMaxIndx;
        plane = getXPlane(start);
        d = VOXEL_X;
        if(pixel.x - source.x < 0){
            plane = getXPlane(end);
            direction = -1;
        }
    } else if(ax == Y){
        start = planeIndexRange.yMinIndx;
        end = planeIndexRange.yMaxIndx;
        plane = getYPlane(start);
        d = VOXEL_Y;
        if(pixel.y - source.y < 0){
            plane = getYPlane(end);
            direction = -1;
        }
    } else if(ax == Z) {
        start = planeIndexRange.zMinIndx;
        end = planeIndexRange.zMaxIndx;
        plane = getZPlane(start);
        d = VOXEL_Z;
        if(pixel.z - source.z < 0){
            plane = getZPlane(end);
            direction = -1;
        }
    } else
        assert(0);

    for (int i = start; i < end; i++){
        getIntersection(source, pixel, plane, ax, &a[i - start]);
        plane += d * direction;
    }

}

/**
 * Merges two sorted arrays into one single sorted array.
 * Returns the length of the merged array.
 * 'a' is a pointer to a sorted array.
 * 'b' is a pointer to a sorted array.
 * 'lenA' is the length of the array pointed by 'a'
 * 'lenB' is the length of the array pointed by 'b'
 * 'c' is a pointer to the array to store the results.
*/
int merge(double *a, double *b, int lenA, int lenB, double *c){
    int i = 0;
    int j = 0;
    int k = 0;
    while (j < lenA && k < lenB) {
        if (a[j] < b[k]) {
            c[i] = a[j];
            j++;
        } else {
            c[i] = b[k];
            k++;
        }
        i++;
    }
    while (j < lenA) {
        c[i] = a[j];
        i++;
        j++;
    }
    while (k < lenB) {
        c[i] = b[k];
        i++;
        k++;
    }
    return lenA + lenB;
}

/**
 * Merges three sorted arrays into one single sorted array.
 * Returns the length of the merged array.
 * 'a' is a pointer to a sorted array.
 * 'b' is a pointer to a sorted array.
 * 'c' is a pointer to a sorted array.
 * 'lenA' is the length of the array pointed by 'a'
 * 'lenB' is the length of the array pointed by 'b'
 * 'lenC' is the length of the array pointed by 'b'
 * 'merged' is a pointer to the array to store the results.
*/
int mergeABC(double *a, double *b, double *c, int lenA, int lenB, int lenC, double *merged){
    double ab[lenA + lenB];
    merge(a, b, lenA, lenB, ab);
    return merge(ab, c, lenA + lenB, lenC, merged);
}

/**
 * Returns the cartesian coordinates of a pixel located in row 'r' and column 'c' of the detector.
 * The detector is located at its angular position of angle 'angle'.
 * 'r' is the row of the pixel on the detector matrix.
 * 'c' is the column of the pixel on the detector matrix.
 * 'angle' is the index of the angle the detector is rotated of. Index '0' corresponds to the least
 * value the detector can be rotated of.
*/
struct point getPixel(int r, int c, int angle){
    struct point pixel;
    const double sinAngle = sin_table[angle];
    const double cosAngle = cos_table[angle];

    pixel.x = (DOD * sinAngle) + cosAngle * (-elementOffset + PIXEL * c);
    pixel.y = (-DOD) * cosAngle + sinAngle * (-elementOffset + PIXEL * c);
    pixel.z = -elementOffset + PIXEL * r;

    return pixel;
}

/**
 * Returns the minimum value between 'a' and 'b'.
 */
int min(int a, int b){
    return a < b ? a : b;
}

/**
 * Computes the projection of a sub-section of the object onto the detector for each source position.
 * 'slice' is the index of the sub-section of the object.
 * 'f' stores the coefficients of the voxels cointained in the sub-section.
 * 'absorbment' is the resulting array, contains the value of absorbtion for each pixel.
 * 'absMax' is the maximum absorbtion computed.
 * 'absMax' is the minimum absorbtion computed.
*/
void computeProjections(int slice, double *f, double *absorbment, double *absMax, double *absMin){
    const int nTheta = (int)(AP / STEP_ANGLE);                      //number of angular position
    double amax = -INFINITY;
    double amin = INFINITY;
    double temp[3][2];
    double aMerged[nPlanes[X] + nPlanes[X] + nPlanes[X]];
    double segments;
    double aX[nPlanes[X]];
    double aY[nPlanes[Y]];
    double aZ[nPlanes[Z]];

    //iterates over each source
    for(int angle = 0; angle <= nTheta; angle++){
        struct point source;
        source.z = 0;
        source.x = -sin((AP / 2 - angle * STEP_ANGLE) * M_PI / 180) * DOS;
        source.y = cos((AP / 2 - angle * STEP_ANGLE) * M_PI / 180) * DOS;

        //iterates over each pixel of the detector 
#pragma omp parallel for collapse(2) schedule(dynamic) default(none) shared(nSidePixels, angle, source, slice, f, absorbment, stationaryDetector, nTheta, nVoxel) private(temp, aX, aY, aZ, aMerged, segments) reduction(min:amin) reduction(max:amax)
        for(int r = 0; r < nSidePixels; r++){
            for(int c = 0; c < nSidePixels; c++){
                struct point pixel;

                //gets the pixel position based on whether the detector rotates or not
                if(stationaryDetector){
                    pixel = getPixel(r,c,angle);
                } else {
                    pixel = getPixel(r,c, nTheta / 2);
                }

                //computes Min-Max parametric value 
                double aMin, aMax;
                int isXParallel = 0, isYParallel = 0, isZParallel = 0;
                isXParallel = getSidesIntersection(source, pixel, &temp[X][0], &temp[X][1], X, slice);
                isYParallel = getSidesIntersection(source, pixel, &temp[Y][0], &temp[Y][1], Y, slice);
                isZParallel = getSidesIntersection(source, pixel, &temp[Z][0], &temp[Z][1], Z, slice);
                int isParallel = -1;
                isParallel = isXParallel == -1 ? isParallel : isXParallel; 
                isParallel = isYParallel == -1 ? isParallel : isYParallel; 
                isParallel = isZParallel == -1 ? isParallel : isZParallel; 

                aMin = getAMin(temp, isParallel);
                aMax = getAMax(temp, isParallel);

                if(aMin < aMax){
                    //computes Min-Max plane indexes 
                    const struct ranges indeces = getRangeOfIndex(source, pixel, isParallel, aMin, aMax);

                    //computes lenghts of the arrays containing parametric value of the intersection with each set of parallel planes
                    int lenX = indeces.xMaxIndx - indeces.xMinIndx;
                    int lenY = indeces.yMaxIndx - indeces.yMinIndx;
                    int lenZ = indeces.zMaxIndx - indeces.zMinIndx;
                    lenX = lenX < 0 ? 0 : lenX;
                    lenY = lenY < 0 ? 0 : lenY;
                    lenZ = lenZ < 0 ? 0 : lenZ;
                    const int lenA = lenX + lenY + lenZ;

                    //computes ray-planes intersection Nx + Ny + Nz
                    getAllIntersections(source, pixel, indeces, aX, X);
                    getAllIntersections(source, pixel, indeces, aY, Y);
                    getAllIntersections(source, pixel, indeces, aZ, Z);

                    //computes segments Nx + Ny + Nz
                    mergeABC(aX, aY, aZ, lenX, lenY, lenZ, aMerged);

                    //associates each segment to the respective voxel Nx + Ny + Nz
                    const double d12 = sqrt(pow(pixel.x - source.x, 2) + pow(pixel.y - source.y, 2) + pow(pixel.z - source.z, 2));
                    const int pixelIndex = angle * nSidePixels * nSidePixels + r *nSidePixels + c;
                    for(int i = 0; i < lenA - 1; i ++){
                        segments = d12 * (aMerged[i + 1] - aMerged[i]);
                        const double aMid = (aMerged[i + 1] + aMerged[i]) / 2;
                        const int xRow = min((int)((source.x + aMid * (pixel.x - source.x) - getXPlane(0)) / VOXEL_X), nVoxel[X] - 1);
                        const int yRow = min((int)((source.y + aMid * (pixel.y - source.y) - getYPlane(0)) / VOXEL_Y), nVoxel[Y] - 1);
                        const int zRow = min((int)((source.z + aMid * (pixel.z - source.z) - getZPlane(0)) / VOXEL_Z), nVoxel[Z] - 1);

                        absorbment[pixelIndex] += f[(yRow - slice) * nVoxel[X] * nVoxel[Z] + zRow * nVoxel[Z] + xRow] * segments;
                    }
                    amax = fmax(amax, absorbment[pixelIndex]);
                    amin = fmin(amin, absorbment[pixelIndex]);

                }
            }
        }
        // fprintf(stderr,"Time[%d]: %lf\n",angle, omp_get_wtime() - time);
    }
    *absMax = amax;
    *absMin = amin;
}

int main(int argc, char *argv[])
{

    //number of angular positions
    const int nTheta = (int)(AP / STEP_ANGLE);
    //array containing the coefficents of each voxel
    double *f = (double*)malloc(sizeof(double) * nVoxel[X] * nVoxel[Z] * OBJ_BUFFER);
    //array containing the computed absorption detected in each pixel of the detector
    double *absorbment = (double*)calloc(nSidePixels * nSidePixels * (nTheta + 1), sizeof(double));
    //each thread has its own variable to store its minimum and maximum absorption computed
    double absMaxValue, absMinValue;

    init_tables();
    
    double totalTime = omp_get_wtime();

    int objectType = 0;
    if(argc > 3){
        fprintf(stderr,"Usage: %s [0-1] [object Type]\n First parameter must be 0 or 1, it indicates whether to rotate the detector (0) or keep it stationary (1); object type must be 1,2 or 3.",argv[0]);
        return EXIT_FAILURE;
    }
    if(argc > 1){
        stationaryDetector = atoi(argv[1]);
    }
    if(argc > 2){
        objectType = atoi(argv[2]);
    }

    //iterates over object subsection
    for(int slice = 0; slice < nVoxel[Y]; slice += OBJ_BUFFER){
        //generate object subsection
        switch (objectType){
            case 1:
                generateCubeWithSphereSlice(f, OBJ_BUFFER, slice, nVoxel[X]);
                break;
            case 2:
                generateSphereSlice(f, OBJ_BUFFER, slice, VOXEL_MAT / 2);
                break;
            default:
                generateCubeSlice(f, OBJ_BUFFER, slice, nVoxel[X]);
                break;
        }

        //computes subsection projection
        computeProjections(slice, f, absorbment, &absMaxValue, &absMinValue);
    }
    fprintf(stderr,"Execution time: %lf\n", omp_get_wtime() - totalTime);
    fflush(stderr);

    //iterates over each absorption value computed, prints a value between [0-255]
    printf("P2\n%d %d\n255", nSidePixels, nSidePixels * (nTheta + 1));
    for(double angle = 0; angle <= nTheta; angle ++){
        for(int i = 0; i < nSidePixels; i++ ){
            printf("\n");
            for(int j = 0; j < nSidePixels; j++ ){
                int pixelIndex = angle * nSidePixels * nSidePixels + i *nSidePixels + j;
                int color =  (absorbment[pixelIndex] - absMinValue) * 255 / (absMaxValue - absMinValue);
                printf("%d ", color);
            }
        }
    }

    free(f);
    free(absorbment);

}
