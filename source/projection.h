#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "common.h"

#define OBJ_BUFFER 1000


/**
 * Computes the projection of a sub-section of the object onto the detector for each source position.
 * 'slice' is the index of the sub-section of the object.
 * 'f' is an array stores the coefficients of the voxels cointained in the sub-section.
 * 'absorbment' is the resulting array, contains the value of absorbtion for each pixel.
 * 'absMax' is the maximum absorbtion computed.
 * 'absMax' is the minimum absorbtion computed.
*/
void computeProjections(int slice, double *f, double *absorbment, double *absMax, double *absMin);