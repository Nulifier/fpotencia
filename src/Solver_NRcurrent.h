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
	class Solver_NRcurrent {
	public:

		Solver_NRcurrent(Circuit model);

		Solver_NRcurrent(Circuit model, cx_solution sol_);

		virtual ~Solver_NRcurrent();

		/*Properties*/
		Circuit Model;

		double EPS = 1e-6;

		int Iterations = 0;

		unsigned int maxIterations = 2;

		Solver_State solve(); //Solves the grid

		void update_solution_power_from_circuit();

	private:
		
		bool Debug = true;

		std::vector<int> PQPV;

		cx_solution Sol;

		vec Pesp;

		vec Qesp;

		/*
		 * This function returns the calculated increments of y
		 */
		void inc_y(vec &x, cx_solution &sol, uint N);


		/*Calculate the a, b, c & d parameters
		 */
		void abcd(uint k, cx_solution &sol, double &a, double &b, double &c, double &d);


		/*
		 * Calculates the Jacobian wich is passed by refference
		 */
		void Jacobian(mat &J, cx_solution &sol, uint N, bool updating);


		/*
		 */
		bool converged(vec &X, uint Nj); //check if the solution converged


		/*
		 */
		void update_solution(cx_solution &sol, vec &x, uint N);


		/*
		 */
		void calculate_slack_power(); //calculate the slack bus power   


		/*
		 */
		bool checks(); //check the solvability with this method

		/*
		 */
		void fill_especifyed_values();
	};
}
