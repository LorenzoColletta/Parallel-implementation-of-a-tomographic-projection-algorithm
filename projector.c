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
    int minIndx;
    int maxIndx;
};

double sin_table[1024], cos_table[1024];

int VOXEL_MAT;
int DETECTOR;
int DOD;
int DOS;

// number of voxel along X axis [0], Y axis [1], Z Axis [2]
int nVoxel[3];

//number of parallel planes alogn X axis [0], Y axis [1], Z axis [2]
int nPlanes[3];

//position of the top-left pixel of the detector relative to the center of the detector
double elementOffset;

//number of pixel along a detector's side
int nSidePixels;


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
 * Returns the minimum value between 'a' and 'b'.
 */
int min(int a, int b){
    return a < b ? a : b;
}

/**
 * Returns the minimum value between 'a', 'b' and 'c'
 */
int minABC(int a, int b, int c){
    return min(a,min(b,c));
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
#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, nVoxel, offset, diameter, VOXEL_MAT)
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

#pragma omp parallel for collapse(3) default(none) shared(f, nOfSlices, offset, sideLength, nVoxel, innerToOuterDiff, rightSide, sphereCenter, VOXEL_MAT)
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
int getIntersection(double a, double b, double *plane, int nPlanes, double *inters){
    if(a - b != 0){
        for(int i = 0; i < nPlanes; i++){
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
struct ranges getRangeOfIndex(const double source, const double pixel, int isParallel, double aMin, double aMax, enum axis ax){
    struct ranges idxs;
    double firstPlane, lastPlane;
    int voxelDim;

    if(ax == X){
        voxelDim = VOXEL_X;
        firstPlane = getXPlane(0);
        lastPlane = getXPlane(nPlanes[X] - 1);
    } else if( ax == Y){
        voxelDim = VOXEL_Y;
        firstPlane = getYPlane(0);
        lastPlane = getYPlane(nPlanes[Y] - 1);
    } else {
        voxelDim = VOXEL_Z;
        firstPlane = getZPlane(0);
        lastPlane = getZPlane(nPlanes[Z] - 1);
    }

    //gets range of indeces of XZ parallel planes
    if(isParallel != Y){
        if(pixel - source >= 0){
            idxs.minIndx = nPlanes[ax] - ceil((lastPlane - aMin * (pixel - source) - source) / voxelDim);
            idxs.maxIndx = 1 + floor((aMax * (pixel - source) + source - firstPlane) / voxelDim);
        } else {
            idxs.minIndx = nPlanes[ax] - ceil((lastPlane - aMax * (pixel - source) - source) / voxelDim);
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
void getAllIntersections(const double source, const double pixel, const struct ranges planeIndexRange, double *a, enum axis ax){
    int start = 0, end = 0;
    int direction = 1;
    double d;

    start = planeIndexRange.minIndx;
    end = planeIndexRange.maxIndx;
    double plane[end - start];
    if(ax == X){
        plane[0] = getXPlane(start);
        d = VOXEL_X;
        if(pixel - source < 0){
            plane[0] = getXPlane(end);
            d = -VOXEL_X;
        }
    } else if(ax == Y){
        plane[0] = getYPlane(start);
        d = VOXEL_Y;
        if(pixel - source < 0){
            plane[0] = getYPlane(end);
            d = -VOXEL_Y;
        }
    } else if(ax == Z){
        plane[0] = getZPlane(start);
        d = VOXEL_Z;
        if(pixel - source < 0){
            plane[0] = getZPlane(end);
            d = -VOXEL_Z;
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
 * Returns the planes of the object's sides orthogonal to the 'x' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesXPlanes(double *planes){
    planes[0] = getXPlane(0);
    planes[1] = getXPlane(nPlanes[X] - 1);
}

/**
 * Returns the planes of the object's sides orthogonal to the 'y' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesYPlanes(double *planes, int slice){
    planes[0] = getYPlane(slice);
    planes[1] = getYPlane(min(nPlanes[Y] - 1, OBJ_BUFFER + slice));
}

/**
 * Returns the planes of the object's sides orthogonal to the 'z' axis.
 * 'planes' is the pointer to an array of two elements.
 */
void getSidesZPlanes(double *planes){
    planes[0] = getZPlane(0);
    planes[1] = getZPlane(nPlanes[Z] - 1);
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
double computeAbsorption(struct point source, struct point pixel, int angle, double *a, int lenA, int slice, double *f){
    const double d12 = sqrt(pow(pixel.x - source.x, 2) + pow(pixel.y - source.y, 2) + pow(pixel.z - source.z, 2));
    
    double absorbment = 0.0;

    for(int i = 0; i < lenA - 1; i ++){
        const double segments = d12 * (a[i + 1] - a[i]);
        const double aMid = (a[i + 1] + a[i]) / 2;
        const int xRow = min((int)((source.x + aMid * (pixel.x - source.x) - getXPlane(0)) / VOXEL_X), nVoxel[X] - 1);
        const int yRow = minABC((int)((source.y + aMid * (pixel.y - source.y) - getYPlane(slice)) / VOXEL_Y), nVoxel[Y] - 1, slice + OBJ_BUFFER - 1);
        const int zRow = min((int)((source.z + aMid * (pixel.z - source.z) - getZPlane(0)) / VOXEL_Z), nVoxel[Z] - 1);

        absorbment += f[(yRow) * nVoxel[X] * nVoxel[Z] + zRow * nVoxel[Z] + xRow] * segments;
    }
    return absorbment;
}


/**
 * Computes the projection of a sub-section of the object onto the detector for each source position.
 * 'slice' is the index of the sub-section of the object.
 * 'f' is an array stores the coefficients of the voxels cointained in the sub-section.
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
                    pixel = getPixel(r,c, nTheta / 2);
                } else {
                    pixel = getPixel(r,c,angle);
                }


                //computes Min-Max parametric values
                double aMin, aMax;
                double temp2[2];
                int isParallel = -1;
                getSidesXPlanes(temp2);
                if(!getIntersection(source.x, pixel.x, temp2, 2, &temp[X][0])){
                    isParallel = X;
                }
                getSidesYPlanes(temp2, slice);
                if(!getIntersection(source.y, pixel.y, temp2, 2, &temp[Y][0])){
                    isParallel = Y;
                }
                getSidesZPlanes(temp2);
                if(!getIntersection(source.z, pixel.z, temp2, 2, &temp[Z][0])){
                    isParallel = Z;
                }

                aMin = getAMin(temp, isParallel);
                aMax = getAMax(temp, isParallel);

                if(aMin < aMax){
                    //computes Min-Max plane indexes 
                    struct ranges indeces[3];
                    indeces[X] = getRangeOfIndex(source.x, pixel.x, isParallel, aMin, aMax, X);
                    indeces[Y] = getRangeOfIndex(source.y, pixel.y, isParallel, aMin, aMax, Y);
                    indeces[Z] = getRangeOfIndex(source.z, pixel.z, isParallel, aMin, aMax, Z);

                    //computes lenghts of the arrays containing parametric value of the intersection with each set of parallel planes
                    int lenX = indeces[X].maxIndx - indeces[X].minIndx;
                    int lenY = indeces[Y].maxIndx - indeces[Y].minIndx;
                    int lenZ = indeces[Z].maxIndx - indeces[Z].minIndx;
                    lenX = lenX < 0 ? 0 : lenX;
                    lenY = lenY < 0 ? 0 : lenY;
                    lenZ = lenZ < 0 ? 0 : lenZ;
                    const int lenA = lenX + lenY + lenZ;

                    //computes ray-planes intersection Nx + Ny + Nz
                    getAllIntersections(source.x, pixel.x, indeces[X], aX, X);
                    getAllIntersections(source.y, pixel.y, indeces[Y], aY, Y);
                    getAllIntersections(source.z, pixel.z, indeces[Z], aZ, Z);

                    //computes segments Nx + Ny + Nz
                    mergeABC(aX, aY, aZ, lenX, lenY, lenZ, aMerged);

                    //associates each segment to the respective voxel Nx + Ny + Nz
                    const int pixelIndex = angle * nSidePixels * nSidePixels + r *nSidePixels + c;
                    absorbment[pixelIndex] += computeAbsorption(source, pixel, angle, aMerged, lenA, slice, f);
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


    int n = 2352;
    int objectType = 0;
    if(argc > 4){
        fprintf(stderr,"Usage: %s [0-1] [object Type]\n First parameter must be 0 or 1, it indicates whether to rotate the detector (0) or keep it stationary (1); object type must be 1,2 or 3.",argv[0]);
        return EXIT_FAILURE;
    }
    if(argc > 1){
        n = atoi(argv[1]);
    }
    if(argc > 2){
        stationaryDetector = atoi(argv[2]);
    }
    if(argc > 3){
        objectType = atoi(argv[3]);
    }

    VOXEL_MAT = n * VOXEL_X * 125 / 294;
    DETECTOR = n * PIXEL;
    DOD = 1.5 * VOXEL_MAT;
    DOS = 6 * VOXEL_MAT;

    nVoxel[X] = VOXEL_MAT / VOXEL_X;
    nVoxel[Y] = VOXEL_MAT / VOXEL_Y;
    nVoxel[Z] = VOXEL_MAT / VOXEL_Z;

    nPlanes[X] = (VOXEL_MAT / VOXEL_X) + 1;
    nPlanes[Y] = (VOXEL_MAT / VOXEL_Y) + 1;
    nPlanes[Z] = (VOXEL_MAT / VOXEL_Z) + 1;

    elementOffset = DETECTOR / 2 - PIXEL / 2;

    nSidePixels = DETECTOR / PIXEL;

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
