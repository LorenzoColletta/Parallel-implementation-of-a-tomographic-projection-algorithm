/****************************************************************************
 *
 * voxel.h
 *
 * Copyright (C) 2024 Lorenzo Colletta
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ****************************************************************************/

#ifndef VOXEL_H
#define VOXEL_H

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

#endif
