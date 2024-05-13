# 3D-image-reconstruction-based-on-tomographic-data
## Usage
### Compile  
    gcc -std=c99 -Wall -Wpedantic -fopenmp  projector.c -lm -o projector
### Run
    ./projector 0 1 > CubeWithSphere.pgm
### Convert
    convert CubeWithSphere.pgm CubeWithSphere.jpeg