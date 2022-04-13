/*
 * File:   Solver_NRpolar.h
 * Author: Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
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
	 * This class implements the Nerwton Raphson method
	 * in polar coordinates to solve the circuit.
	 */
	class Solver_NRrect {
	public:

		Solver_NRrect(Circuit model);

		Solver_NRrect(Circuit model, cx_solution sol_);

		virtual ~Solver_NRrect();

		/*Properties*/
		Circuit Model;

		double EPS = 1e-6;

		int Iterations = 0;

		int Max_Iter = 100;

		Solver_State solve(); //Solves the grid

		void update_solution_power_from_circuit();
		
	private:

		vector<int> BUSES;

		vector<int> PQPV;

		vector<int> PQPVPQPV;

		sp_vec Pesp;
		
		sp_vec Qesp;
		
		sp_vec V2esp;

		cx_solution Sol;

		void Jacobian(mat &J, uint npq, uint npv); //calculate the jacobian, J is passed by refference

		void get_mismatches(vec &inc, uint npq, uint npv); //inc is passed by refference

		void calculate_Q(uint npq, uint npv); //calculate the reative power at the PV buses

		double Q(uint k);

		double P(uint k);
		
		double c(uint i);
		
		double d(uint i);

		bool converged(vec PQinc, uint npqpvpq, double &error); //check if the solution converged

		void update_solution(vec X, uint npq, uint npv);

		void calculate_slack_power(); //calculate the slack bus power        

		bool checks(); //check the solvability with this method

		void fill_especifyed_values();
		
	};
}
