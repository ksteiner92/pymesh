/**
 * This file is part of edge.
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
 * eigeh.h
 *
 *  Created on: 05.06.2017
 *      Author: Klaus Steiner
 */

#ifndef EDGE_EIGEN_H
#define EDGE_EIGEN_H

/**
 * Workaround for Bug:
 * http://eigen.tuxfamily.org/bz/show_bug.cgi?id=874
 */
#ifndef MKL_BLAS
#define MKL_BLAS -1
#endif
#define EIGEN_USE_MKL_ALL
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#endif //EDGE_EIGEN_H
