/*
 * File:   Solver_Jacobi.cpp
 * Author: Santiago Peñate Vera
 *
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "Solver_Jacobi.h"

namespace fPotencia {

	/*
	 * Class constructor
	 */
	Solver_Jacobi::Solver_Jacobi(Circuit model) {
		Model = model;

		Sol = Model.get_initial_cx_solution();
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_cx_solution();
		}
		
		BUSES.reserve(Model.loadBusIndices.size()); // preallocate memory
		BUSES.insert(BUSES.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		//BUSES.insert(BUSES.end(), Model.PV_list.begin(), Model.PV_list.end());

	}

	/*
	 * Class destructor
	 */
	Solver_Jacobi::~Solver_Jacobi() {
	}

	/*
	 * Solve the circuit using the Jacibi method
	 */
	Solver_State Solver_Jacobi::solve() {

		Solver_State state = Solver_State::Not_Converged;
		uint n = Model.buses.size();
		Iterations = 0;
		bool converged = false;
		double Vdiff_max;

		while (!converged && Iterations <= Max_Iter) {

			Vdiff_max = 0.0;

			for (uint k : BUSES) { //This loop can be executed in parallel
				//Calculate bus currents
				cx_double I = conj(Sol.S.coeff(k) / Sol.V.coeff(k));

				//Calculate other buses current contribution
				cx_double Iinj(0.0, 0.0);
				for (uint j = 0; j < n; j++)
					if (k != j)
						Iinj += Model.Y.coeff(k, j) * Sol.V.coeff(j);

				//Calculate the new voltage
				cx_double Vk_new = (1.0 / Model.Y.coeff(k, k)) * (I - Iinj);

				//Calculaton of the voltage difference for the convergence check
				double Vdiff = abs(Vk_new - Sol.V.coeff(k)); //might lead to not compare the voltage angles...check
				if (Vdiff > Vdiff_max)
					Vdiff_max = Vdiff;

				//Asign the new voltage
				Sol.V(k) = Vk_new;
			}

			//Convergence check
			//cout << iter << ":maxDiff: " << Vdiff_max << endl;
			if (Vdiff_max <= EPS)
				converged = true;

			Iterations++;
		}

		if (converged) {
			state = Solver_State::Converged;

			//Calculate the slack bus power
			for (uint k : Model.slackBusIndices) {
			cx_double I(0.0, 0.0);
			for (uint j = 0; j < Model.buses.size(); j++) {
				I += Model.Y.coeff(k, j) * Sol.V(j);
			}
			Sol.S(k) = Sol.V(k) * conj(I); //now this is the power
		}


			//Set the solution to the circuit
			Model.set_solution(Sol);
		} else if (Iterations == Max_Iter)
			state = Solver_State::Not_Converged;

		return state;
	}
}
