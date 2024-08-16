/**
 * Generates a sub-section of a solid cubic object given its side length.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'sideLength' length of the side of the cubic object
*/
void generateCubeSlice(double *f, int nOfSlices, int offset, int sideLength);



/**
 * Generates a sub-section of a solid spherical object given its diameter.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'diameter' is the diameter of the sphere.
*/
void generateSphereSlice(double *f, int nOfSlices, int offset, int diameter);


/**
 * Generates a sub-section of a solid cubic object with an internal spherical cavity.
 * 'f' is the pointer to the array on which to store the sub-section.
 * 'nOfSlices' is the number of voxel along the Y axis.
 * 'offset' distance (in number of voxel) of the slices to be generated from the initial slice.
 * 'sideLength' lenght of the side of the object.
*/
void generateCubeWithSphereSlice(double *f, int nOfSlices, int offset, int sideLength);