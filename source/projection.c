#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "common.h"
#include "projection.h"

extern double *sineTable;
extern double *cosineTable;

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
        sinTable[positionIndex] = sin((- gl_angularTrajectory / 2 + positionIndex * gl_positionsAngularDistance) * M_PI / 180);
        cosTable[positionIndex] = cos((- gl_angularTrajectory / 2 + positionIndex * gl_positionsAngularDistance) * M_PI / 180);
    }
}


/**
 * Returns the minimum value between 'a' and 'b'.
 */
int min(int a, int b){
    return a < b ? a : b;
}

/**
 * Returns the minimum value between 'a', 'b' and 'c'
 */
int min3(int a, int b, int c){
    return min(a,min(b,c));
}

/**
 * returns the coordinate of a plane parallel to the YZ plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getXPlane(int index){
    return -((gl_objectSideLenght) / 2) + index * gl_voxelXDim;
}

/**
 * returns the coordinate of a plane parallel to the XZ plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getYPlane(int index){
    return -((gl_objectSideLenght) / 2) + index * gl_voxelYDim;
}

/**
 * returns the coordinate of a plane parallel to the XY plane
 * 'index' is the index of the plane to be returned where '0' is the index of the smallest-valued coordinate plane
*/
double getZPlane(int index){
    return -((gl_objectSideLenght) / 2) + index * gl_voxelZDim;
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
 * Computes the parametrical value of the intersection between a ray and a plane given the coordinate of the plane.
 * Returns the axis to which the ray is orthogonal, -1 otherwise.
 * 'a' and 'b' are the components along the axis ortogonal to the plane of the two points defining the ray.
 * 'value' is the coordinate of the plane.
 * 'ax' is the axis orthogonal to the plane.
 * 'inters' is the parametrical value that identify the intersection point along the ray.
*/
int getIntersection(double a, double b, double *plane, int gl_nPlanes, double *inters){
    if(a - b != 0){
        for(int i = 0; i < gl_nPlanes; i++){
            inters[i] = (plane[i] - a) / (b - a);
        }
        return 1;
    }
    return 0;
}

/**
 * Returns the range of indices of the planes.
 * 'source' and 'pixel' are the components along the axis ortogonal to the plane of the two points defining the ray .
 * 'isParallel'  has a value corrisponding to the axis to which the array is orthogonal, -1 otherwise.
 * 'aMin' is the minimum parametrical value of the intersection between the ray and the object.
 * 'aMax' is the maximum parametrical value of the intersection between the ray and the object.
 * 'ax' is the axis orthogonal to the plane.
*/
Ranges getRangeOfIndex(const double source, const double pixel, int isParallel, double aMin, double aMax, Axis ax){
    Ranges idxs;
    double firstPlane, lastPlane;
    int voxelDim;

    if(ax == X){
        voxelDim = gl_voxelXDim;
        firstPlane = getXPlane(0);
        lastPlane = getXPlane(gl_nPlanes[X] - 1);
    } else if( ax == Y){
        voxelDim = gl_voxelYDim;
        firstPlane = getYPlane(0);
        lastPlane = getYPlane(gl_nPlanes[Y] - 1);
    } else {
        voxelDim = gl_voxelZDim;
        firstPlane = getZPlane(0);
        lastPlane = getZPlane(gl_nPlanes[Z] - 1);
    }

    //gets range of indeces of XZ parallel planes
    if(isParallel != Y){
        if(pixel - source >= 0){
            idxs.minIndx = gl_nPlanes[ax] - ceil((lastPlane - aMin * (pixel - source) - source) / voxelDim);
            idxs.maxIndx = 1 + floor((aMax * (pixel - source) + source - firstPlane) / voxelDim);
        } else {
            idxs.minIndx = gl_nPlanes[ax] - ceil((lastPlane - aMax * (pixel - source) - source) / voxelDim);
            idxs.maxIndx = floor((aMin * (pixel - source) + source - firstPlane) / voxelDim);
        }
    } else {
        idxs.minIndx = 0;
        idxs.maxIndx = 0;
    }
    return idxs;
}



/**
 * Computes each parametric value of the intersection between the ray and the planes whose index is in the range planeIndexRange.
 * 'source' and 'pixel' are the components along the axis ortogonal to the plane of the two points defining the ray .
 * 'planeIndexRange' is a structure containing the ranges of indeces of planes.
 * 'a' is a pointer to the array on which to store the parametrical values.
 * 'ax' is the axis orthogonal to the set of planes to which compute the intersection.
*/
void getAllIntersections(const double source, const double pixel, const Ranges planeIndexRange, double *a, Axis ax){
    int start = 0, end = 0;
    double d;

    start = planeIndexRange.minIndx;
    end = planeIndexRange.maxIndx;
    double plane[end - start];
    if(ax == X){
        plane[0] = getXPlane(start);
        d = gl_voxelXDim;
        if(pixel - source < 0){
            plane[0] = getXPlane(end);
            d = -gl_voxelXDim;
        }
    } else if(ax == Y){
        plane[0] = getYPlane(start);
        d = gl_voxelYDim;
        if(pixel - source < 0){
            plane[0] = getYPlane(end);
            d = -gl_voxelYDim;
        }
    } else if(ax == Z){
        plane[0] = getZPlane(start);
        d = gl_voxelZDim;
        if(pixel - source < 0){
            plane[0] = getZPlane(end);
            d = -gl_voxelZDim;
        }
    } else assert(0);

    for (int i = 1; i < end - start; i++){
        plane[i] = plane[i-1] + d;
    }
    getIntersection(source, pixel, plane, end - start, a);
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
int merge3(double *a, double *b, double *c, int lenA, int lenB, int lenC, double *merged){
    double ab[lenA + lenB];
    merge(a, b, lenA, lenB, ab);
    return merge(ab, c, lenA + lenB, lenC, merged);
}

/**
 * Returns the cartesian coordinates of the source.
 * 'index' is the index of the position starting from the position with the least angular distance from the y-axis.
 */
Point getSource(int index){
    Point source;
    
    source.z = 0;
    source.x = sineTable[index] * gl_distanceObjectSource;
    source.y = cosineTable[index] * gl_distanceObjectSource;
    
    return source;
}

/**
 * Returns the cartesian coordinates of a pixel located in row 'r' and column 'c' of the detector.
 * The detector is located at index-th angular position.
 * 'r' is the row of the pixel on the detector matrix.
 * 'c' is the column of the pixel on the detector matrix.
 * 'index' is the index of the position of the detector starting from the position with the 
 * least angular distance from the y-axis.
*/
Point getPixel(int r, int c, int index){
    Point pixel;
    const double sinAngle = sineTable[index];
    const double cosAngle = cosineTable[index];
    const double elementOffset =  gl_detectorSideLength / 2 - gl_pixelDim / 2;

    pixel.x = (-gl_distanceObjectDetector * sinAngle) + cosAngle * (-elementOffset + gl_pixelDim * c);
    pixel.y = (-gl_distanceObjectDetector) * cosAngle - sinAngle * (-elementOffset + gl_pixelDim * c);
    pixel.z = -elementOffset + gl_pixelDim * r;

    return pixel;
}

/**
 * Returns the planes of the object's sides orthogonal to the 'x' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesXPlanes(double *planes){
    planes[0] = getXPlane(0);
    planes[1] = getXPlane(gl_nPlanes[X] - 1);
}

/**
 * Returns the planes of the object's sides orthogonal to the 'y' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesYPlanes(double *planes, int slice){
    planes[0] = getYPlane(slice);
    planes[1] = getYPlane(min(gl_nPlanes[Y] - 1, OBJ_BUFFER + slice));
}

/**
 * Returns the planes of the object's sides orthogonal to the 'z' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesZPlanes(double *planes){
    planes[0] = getZPlane(0);
    planes[1] = getZPlane(gl_nPlanes[Z] - 1);
}

/**
 * Computes the absorption the radiological path of a ray given two points.
 * 'source' and 'pixel' are the points defining the ray.
 * 'angle' is an index that numbers the current position.
 * 'a' is an array containing the paramerter values of the intersections between the ray and voxels.
 * 'lenA' is the length of the array 'a'.
 * 'slice' is the index of the sub-section of the object.
 * 'f' is an array that stores the coefficients of the voxels cointained in the sub-section.
 */
double computeAbsorption(Point source, Point pixel, int angle, double *a, int lenA, int slice, double *f){
    const double d12 = sqrt(pow(pixel.x - source.x, 2) + pow(pixel.y - source.y, 2) + pow(pixel.z - source.z, 2));
    
    double absorbment = 0.0;

    for(int i = 0; i < lenA - 1; i ++){
        const double segments = d12 * (a[i + 1] - a[i]);
        const double aMid = (a[i + 1] + a[i]) / 2;
        const int xRow = min((int)((source.x + aMid * (pixel.x - source.x) - getXPlane(0)) / gl_voxelXDim), gl_nVoxel[X] - 1);
        const int yRow = min3((int)((source.y + aMid * (pixel.y - source.y) - getYPlane(slice)) / gl_voxelYDim), gl_nVoxel[Y] - 1, OBJ_BUFFER - 1);
        const int zRow = min((int)((source.z + aMid * (pixel.z - source.z) - getZPlane(0)) / gl_voxelZDim), gl_nVoxel[Z] - 1);

        absorbment += f[(yRow) * gl_nVoxel[X] * gl_nVoxel[Z] + zRow * gl_nVoxel[Z] + xRow] * segments;
    }
    return absorbment;
}


void computeProjections(int slice, double *f, double *absorbment, double *absMax, double *absMin){
    const int nTheta = (int)(gl_angularTrajectory / gl_positionsAngularDistance);                      //number of angular position
    const int nSidePixels = gl_detectorSideLength / gl_pixelDim;

    double amax = -INFINITY;
    double amin = INFINITY;
    double temp[3][2];
    double aMerged[gl_nPlanes[X] + gl_nPlanes[X] + gl_nPlanes[X]];
    double aX[gl_nPlanes[X]];
    double aY[gl_nPlanes[Y]];
    double aZ[gl_nPlanes[Z]];

    //iterates over each source
    for(int positionIndex = 0; positionIndex <= nTheta; positionIndex++){
        const Point source = getSource(positionIndex);

        //iterates over each pixel of the detector 
#pragma omp parallel for collapse(2) schedule(dynamic) default(none) shared(nSidePixels, positionIndex, source, slice, f, absorbment, nTheta, gl_nVoxel) private(temp, aX, aY, aZ, aMerged) reduction(min:amin) reduction(max:amax)
        for(int r = 0; r < nSidePixels; r++){
            for(int c = 0; c < nSidePixels; c++){
                Point pixel;

                //gets the pixel's center cartesian coordinates
                pixel = getPixel(r,c,positionIndex);

                //computes Min-Max parametric values
                double aMin, aMax;
                double sidesPlanes[2];
                int isParallel = -1;
                getSidesXPlanes(sidesPlanes);
                if(!getIntersection(source.x, pixel.x, sidesPlanes, 2, &temp[X][0])){
                    isParallel = X;
                }
                getSidesYPlanes(sidesPlanes, slice);
                if(!getIntersection(source.y, pixel.y, sidesPlanes, 2, &temp[Y][0])){
                    isParallel = Y;
                }
                getSidesZPlanes(sidesPlanes);
                if(!getIntersection(source.z, pixel.z, sidesPlanes, 2, &temp[Z][0])){
                    isParallel = Z;
                }

                aMin = getAMin(temp, isParallel);
                aMax = getAMax(temp, isParallel);

                if(aMin < aMax){
                    //computes Min-Max plane indexes 
                    Ranges indeces[3];
                    indeces[X] = getRangeOfIndex(source.x, pixel.x, isParallel, aMin, aMax, X);
                    indeces[Y] = getRangeOfIndex(source.y, pixel.y, isParallel, aMin, aMax, Y);
                    indeces[Z] = getRangeOfIndex(source.z, pixel.z, isParallel, aMin, aMax, Z);

                    //computes lenghts of the arrays containing parametric value of the intersection with each set of parallel planes
                    int lenX = indeces[X].maxIndx - indeces[X].minIndx;
                    int lenY = indeces[Y].maxIndx - indeces[Y].minIndx;
                    int lenZ = indeces[Z].maxIndx - indeces[Z].minIndx;
                    if(lenX < 0){
                        lenX = 0;
                    }
                    if(lenY < 0){
                        lenY = 0;
                    }
                    if(lenZ < 0){
                        lenZ = 0;
                    }
                    const int lenA = lenX + lenY + lenZ;

                    //computes ray-planes intersection Nx + Ny + Nz
                    getAllIntersections(source.x, pixel.x, indeces[X], aX, X);
                    getAllIntersections(source.y, pixel.y, indeces[Y], aY, Y);
                    getAllIntersections(source.z, pixel.z, indeces[Z], aZ, Z);

                    //computes segments Nx + Ny + Nz
                    merge3(aX, aY, aZ, lenX, lenY, lenZ, aMerged);

                    //associates each segment to the respective voxel Nx + Ny + Nz
                    const int pixelIndex = positionIndex * nSidePixels * nSidePixels + r *nSidePixels + c;
                    absorbment[pixelIndex] += computeAbsorption(source, pixel, positionIndex, aMerged, lenA, slice, f);
                    amax = fmax(amax, absorbment[pixelIndex]);
                    amin = fmin(amin, absorbment[pixelIndex]);

                }
            }
        }
    }
    *absMax = amax;
    *absMin = amin;
}

