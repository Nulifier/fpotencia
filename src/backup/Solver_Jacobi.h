/*
 * File:   Solver_Jacobi.h
 * Author: Santiago Peñate Vera
 *
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include <cmath>
#include "Circuit.h"
#include "Solution.h"

namespace fPotencia {
	/*
	 * This class implemets the Jacobi iterative solver
	 * this method is very similar to the gauss seidel method but
	 * using the impedance matrix instead of the admittance matrix.
	 *
	 * This methos is the simplest available to solve power flows
	 */
	class Solver_Jacobi {
	public:

		Solver_Jacobi(Circuit model);

		virtual ~Solver_Jacobi();

		/*Properties*/
		Circuit Model;

		double EPS = 1e-9;

		int Iterations = 0;

		int Max_Iter = 100;

		Solver_State solve();

	private:

		std::vector<int> BUSES;

		cx_solution Sol;

	};
}
