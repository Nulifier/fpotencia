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
#pragma once

#include <ctime>
#include <iostream>

// EIGEN
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>

// General math
#include <cmath>
#include <complex>
#include <numbers>

// Assertion library
#include <cassert>

#include "enumaratons.h"

namespace fPotencia {
	using cx_double = std::complex<double>;

	using mat = Eigen::MatrixXd;
	using vec = Eigen::VectorXd;

	using cx_vec = Eigen::VectorXcd;
	using cx_mat = Eigen::MatrixXcd;

	using cx_mat3 = Eigen::Matrix3cd; // 3x3 complex matrix

	using sp_cx_mat =
	    Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>;
	using sp_mat = Eigen::SparseMatrix<double>;
	using sp_vec = Eigen::SparseVector<double>;

	using cx_diag_mat =
	    Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic>;

	using uint = unsigned int;
}
