//
// Created by klaus on 2019-11-08.
//

#ifndef PYULB_TYPES_H
#define PYULB_TYPES_H

#include "eigen.h"

typedef unsigned long long int ullong;

typedef long ID;

using EigenDStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
template<typename MatrixType> using EigenDRef = Eigen::Ref<MatrixType, 0, EigenDStride>;
template<typename MatrixType> using EigenDMap = Eigen::Map<MatrixType, 0, EigenDStride>;
using MatrixXid = Eigen::Matrix<ID, Eigen::Dynamic, Eigen::Dynamic>;


#endif //PYULB_TYPES_H
