/* 
 * File:   fpotencia_libs.h
 * Author: Santiago Peñate Vera
 *
 * Created on 11 de noviembre de 2014, 22:43
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 Linux:
 * To satisfy the dependencies right away in Ubuntu (or other linux distro)
 * install:
 *  - libboost-graph-parallel-dev
 *  - libeigen3-dev
 */
#include<iostream>
#include <ctime>

//EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>


//BOOST: graphs
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/breadth_first_search.hpp>
//#include <boost/pending/indirect_cmp.hpp>
//#include <boost/range/irange.hpp>

//General math
#include <complex>
#include <cmath>

//Aseritions library
#include <assert.h> 

#include "enumaratons.h"

typedef std::complex<double> cx_double;

typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

typedef Eigen::VectorXcd cx_vec;
typedef Eigen::MatrixXcd cx_mat;

typedef Eigen::Matrix3cd cx_mat3; // 3x3 complex matrix

typedef Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> sp_cx_mat;
typedef Eigen::SparseMatrix<double> sp_mat;
typedef Eigen::SparseVector<double>  sp_vec;

typedef Eigen::DiagonalMatrix<std::complex<double>,Eigen::Dynamic> cx_diag_mat;


typedef unsigned int uint;


/* FPOTENCIA_LIBS_H */
