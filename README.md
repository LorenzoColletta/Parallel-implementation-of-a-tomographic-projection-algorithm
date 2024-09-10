# Parallel implementation of a tomographic projection algorithm
The purpose of this repository is to provide a parallel implementation of a tomographic projection algorithm based on the Siddon algorithm and an implementation of an input generation algorithm for the projection.
## Generating the input
The input to the projection algorithm is a three-dimensional voxel grid mapped to a one-dimensional array of doubles.
The source code that implements an algorithm to generate the input to the projection algorithm is provided in ```inputgeneration.c```.

### Compile
    gcc -Wall -Wpedantic -std=c99 -fopenmp inputgeneration.c ./source/voxel.c -I./source/ -lm -o inputgeneration

This program can also be compiled so that the output file has no header (the 16 integer values described below). **This output file structure should not be used as input to the projection.c program**. Use the **-DRAW** argument when compiling:

    gcc -Wall -Wpedantic -std=c99 -fopenmp -DRAW inputgeneration.c ./source/voxel.c -I./source/ -lm -o inputgeneration

### Run:
    inputgeneration output.dat [object Type] [integer] 

* First parameter is the name of the file to store the output in;
* Second parameter is optional and can be: 1 (solid cube with spherical cavity), 2 (solid sphere) or 3 (solid cube), if not passed 3 (solid cube) is default;
* Third parameter is the number of pixel per side of the detector, every other parameter is set based to its value, if no value is given, default values are used;

### Output file structure:
The voxel (three-dimensional) grid is represented as a stack of two-dimensional grids. Considering a three-dimensional cartesian system where the x-axis is directed from left to right, the y-axis is directed upwards, and the z-axis is orthogonal to them, a two-dimensional grid can be viewed as a horizontal slice, orthogonal to the y-axis, of the object.

First a sequence of 16 integer values is given, representing in order:
* the side length of a pixel of the detector        
* total angular distance traveled by the source 
* angular distance between each source position
* side length of the object
* side lenth of the detector                   
* distance between the object's center and the detector           
* distance between the object's center and the source position
* voxel side lenght along x-axis
* voxel side lenght along y-axis
* voxel side lenght along z-axis
* the number $v1$ of voxel the object is composed of along the X axis  
* the number $v2$ of voxel the object is composed of along the Y axis 
* the number $v3$ of voxel the object is composed of along the Z axis
* the number of planes along the X axis
* the number of planes along the Y axis
* the number of planes along the Z axis 

Then, the values composing the voxel grid are given for a total of $v1\*v2\*v3$ (double) values.
Each sequence of length $v1*v3$ represents a horizontal slice of the object stored as a one-dimensional array of elements ordered first by the x coordinate and then by the z coordinate.
The first slice memorized is the bottom one, followed by the other slices in ascending order of the y coordinate.

## Computing the projections
### Compile
The projection algorithm can be compiled in two different ways: 
* **with -DBINARY argument:** returns the output in a binary file, in this case the values ​​are stored as type double values ​​for maximum accuracy;

    ```gcc -Wall -Wpedantic -std=c99 -fopenmp -DBINARY projector.c ./source/voxel.c -I./source/ -lm -o projector``` 
  
* **without -DBINARY argument:** it returns the output in a text file in pgm format, in this case the calculated values ​​are converted to a gray scale from 0 to 255.

    ```gcc -Wall -Wpedantic -std=c99 -fopenmp projector.c ./source/voxel.c -I./source/ -lm -o projector``` 

### Run

    projector input.dat output.dat/output.pgm

- First parameter is the name of the input file. The structure of the input file should follow the one described above.
- Second parameter is the name of a text or a binary file to store the output in.

### Binary file structure:
This file structure is used in case the projection algorithm is compiled with ```-DBINARY``` argument.
The computed output is a set of projection images each with a resolution of $n*n$.
The output file is structured as follows:
* the first value is the number of images produced (integer type)
* the second value is $n$ (integer type), the resolution of the image side
* the third value is the maximum value computed (type double)
* the third value is the minimum value computed (type double)
* each sequence of values ​​representing an image is preceded by a value (type double) which indicates the angle from which the image was computed
* an image is a sequence of $n*n$ double values; the image is stored as a one-dimensional array, sorted first by the x coordinate and then by the z coordinate (considering a three-dimensional Cartesian space with the x axis from left to right, the y axis oriented upwards and z perpendicular to them)

### Convert
In case the projection algorithm is compiled to store the output in a pgm file, you can convert the output to another image file format with the following command:
    convert CubeWithSphere.pgm CubeWithSphere.jpeg

## MAKEFILE
This repository provides a makefile to automate the compilation of both ```inputgeneration.c``` and ```projection.c``` source files. To compile the ```projection.c``` source code with the binary file structure, use the ```FILE_BINARY=yes``` argument when running make:

    make FILE_BINARY=yes

To compile ```inputgeneration.c``` so that the output file doesn't include a header use the ```VOXEL_MODEL_RAW=yes``` argument:

    make VOXEL_MODEL_RAW=yes

**In this case, the output file from ```inputgeneration``` should not be used as input to the ```projection``` program.**

Both ```FILE_BINARY=yes``` and ```VOXEL_MODEL_RAW=yes``` arguments can be used simultaneously.

    make FILE_BINARY=yes VOXEL_MODEL_RAW=yes
