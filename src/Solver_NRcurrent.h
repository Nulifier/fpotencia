/*
 * File:   Solver_NRpolar.h
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
#include "Solver.h"
#include "Circuit.h"
#include "Solution.h"

namespace fPotencia {
	/*
	 * This class implements the Newton Raphson method in current equations
	 * Described at:
	 * Developments in the Newton Raphson Power Flow Formulation Based on 
	 * Current Injections
	 * By Vander Menengoy da Costa, Nelson Martins2 and Josk Luiz R. Pereira
	 * 1999
	 */
	class NRcurrentSolver : public Solver {
	public:
		NRcurrentSolver(Circuit& model);

		/**
		 * Construct the solver with an initial solution.
		 */
		NRcurrentSolver(Circuit& model, cx_solution sol_);

		/*Properties*/
		Circuit& Model;

		double tolerance = 1e-6;

		int Iterations = 0;

		unsigned int maxIterations = 2;

		Result solve(bool printIterations = false) override;

		bool canSolve() const override;

		void update_solution_power_from_circuit();

	private:
		std::vector<int> PQPV;

		cx_solution Sol;

		vec Pesp;

		vec Qesp;

		/// This function returns the calculated increments of y.
		void inc_y(vec &x, cx_solution &sol, uint N);

		/// Calculate the a, b, c & d parameters.
		void abcd(uint k, cx_solution &sol, double &a, double &b, double &c, double &d);

		/// Calculates the Jacobian.
		void Jacobian(mat &J, cx_solution &sol, uint N, bool updating);

		/// Check if the solution converged.
		bool converged(vec &X, uint Nj) const;

		/// Generates a solution object from a vector X.
		void update_solution(cx_solution &sol, vec &x, uint N);

		/// Calculate the slack bus power.
		void calculate_slack_power(); //calculate the slack bus power   

		/// Fills the initial power values.
		void fill_specified_values();
	};
}
