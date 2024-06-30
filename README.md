# 3D-image-reconstruction-based-on-tomographic-data
## Usage
### Compile  
    gcc -std=c99 -Wall -Wpedantic -fopenmp  projector.c -lm -o projector
### Run
    ./projector [integer] [0-1] [1-2-3] > image.pgm

The first parameter is the number of mesuring unit per side of the detector.

The second parameter must be either 0 or 1: 
* in case 1 is given, the projections are computed on a rotating detector;
* in case 0 is given, projections are computed on a stationary detector;
* in case no parameter is given, 0 is the default behavior.

The third parameter can be:
* in case 1 is given the computed object is a solid cubic object with an internal spherical cavity;
* in case 2 is given the computed object is a solid spherical object;
* in case no value is given or it is neither 1 nor 2, the computed object is a solid cubic object;

Example:

    ./projector 2352 0 1 > image.pgm
### Convert
    convert CubeWithSphere.pgm CubeWithSphere.jpeg
