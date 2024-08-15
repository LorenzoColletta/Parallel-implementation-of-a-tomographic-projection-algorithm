#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/**
 * The following constants are input parameters values to be considered as a reference.
 */
#define VOXEL_X_DIM 100                         // voxel side length along x-axis
#define VOXEL_Y_DIM 100                         // voxel side length along y-axis
#define VOXEL_Z_DIM 100                         // voxel side length along z-axis

#define PIXEL_DIM 85                            // side length of a pixel of the detector
#define ANGULAR_TRAJECTORY 90                   // total angular distance traveled by the source 
#define POSITIONS_ANGULAR_DISTANCE 15           // angular distance between each source position

#define OBJECT_SIDE_LENGTH                      // side length of the object
#define DETECTOR_SIDE_LENGTH                    // side lenth of the detector
#define DISTANCE_OBJECT_DETECTOR                // distance between the object's center and the detector
#define DISTANCE_OBJECT_SOURCE                  // distance between the object's center and the source position

/**
 * The following global variables represent the input parameters needed to set up the environment.
 * A reference to each variable value is provided by the previously defined constants.
 */
extern int gl_pixelDim;                         // side length of a pixel of the detector        
extern int gl_angularTrajectory;                // total angular distance traveled by the source 
extern int gl_positionsAngularDistance;         // angular distance between each source position
extern int gl_objectSideLenght;                 // side length of the object
extern int gl_detectorSideLength;               // side lenth of the detector                   
extern int gl_distanceDetectorObject;           // distance between the object's center and the detector           
extern int gl_distanceDetectorSource;           // distance between the object's center and the source position
extern int gl_voxelXDim;                        // voxel side along x-axis
extern int gl_voxelYDim;                        // voxel side along y-axis
extern int gl_voxelZDim;                        // voxel side along z-axis
extern int gl_nVoxel[3];                        // an array containing the number of voxel the object is composed of along the X axis [0], Y axis [1], Z axis [2]





