/*
 * File:   Solution.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 8 de agosto de 2014, 10:43
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <complex>
#include <iostream>
//#include <vector>

#include "Solution.h"
#include "fpotencia_libs.h"

namespace fPotencia {
	void cx_solution::resize(unsigned int n) {
		S.setZero(n);
		V.setZero(n);
	}

	void cx_solution::print(const std::string& title) const {
		std::cout << title << std::endl;
		std::cout << "S\t\tV" << std::endl;
		for (uint i = 0; i < length(); i++) {
			std::cout << S[i] << "\t\t" << V[i] << "->" << abs(V[i]) << "^"
			          << arg(V[i]) << std::endl;
		}
	}

	void solution::resize(unsigned int n) {
		P.setZero(n);
		Q.setZero(n);
		Vmag.setZero(n);
		Varg.setZero(n);
	}

	cx_vec solution::V() const {
		cx_vec T(length());
		for (unsigned int i = 0; i < length(); ++i) {
			T[i] = std::polar(Vmag[i], Varg[i]);
		}
		return T;
	}

	void solution::print(const std::string& title) const {
		std::cout << title << std::endl;
		std::cout << "P\tQ\tV\tD" << std::endl;
		for (uint i = 0; i < length(); i++) {
			std::cout << P[i] << "\t" << Q[i] << "\t" << Vmag[i] << "\t" << Varg[i]
			          << std::endl;
		}
	}

	cx_solution solution::toComplex() const {
		cx_solution sol;
		sol.resize(length());
		sol.S.real() = P;
		sol.S.imag() = Q;
		for (uint i = 0; i < length(); i++) {
			sol.V(i) = V(i);
		}
		return sol;
	}
}
