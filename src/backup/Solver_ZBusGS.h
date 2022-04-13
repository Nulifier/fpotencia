/*
 * File:   Solver_ZBlockSubs.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de noviembre de 2014, 23:05
 * Copyright (C) 2014 Santiago Peñate Vera
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
	 * This class implemets the Z-Matrix block susbtitution
	 * this method is very simulat to the gauss seidel method but
	 * using the impedance matrix instead of the admittance matrix.
	 *
	 * The number of iterations can be compared to Newton Raphson, while the
	 * computation per iteration is really low in comparison.
	 * Due to the use of the Z matrix which is full unlike the Y matrix,
	 * this method enhances the convergence.
	 *
	 * On the other hand the Z-matrix is obtained as the inverse of the Y matrix,
	 * but this is only done once, so if many power flow simulations are expected,
	 * the calculation of Z will not vary if the circuit topology remains.
	 */
	class Solver_ZBusGS {
	public:

		Solver_ZBusGS(Circuit model);

		virtual ~Solver_ZBusGS();

		/*Properties*/
		Circuit Model;

		double EPS = 1e-9;

		int Iterations = 0;

		int Max_Iter = 100;


		Solver_State solve();


	private:

		vector<int> BUSES;

		cx_solution Sol;

		cx_vec E0; //voltage in the past iteration
		
		cx_vec E; //voltage in the current iteration

		cx_vec C; //Current vector in the past iteration

		cx_vec I; //Current vector in the current iteration
		
		void calculate_I(cx_vec * V);
		
		void calculate_C(cx_vec * V);

		cx_double calculate_I_k(uint k, cx_vec * V);

		cx_double calculate_V_new_k(uint k, bool pv);

		bool converged();

		bool converged_k(uint k);

		void calculate_slack_power();

		double calculate_Q_k(uint k);

		bool checks();
	};
}
