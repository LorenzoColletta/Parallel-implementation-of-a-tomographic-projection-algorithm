/****************************************************************************
 *
 * projection.h
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

#ifndef PROJECTION_H
#define PROJECTION_H

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

#endif
