/*
 * File:   Solver_NRpolar.cpp
 * Author: Santiago Peñate Vera
 *
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "Solver_NRcurrent.h"

namespace fPotencia {

	/* Constructor
	 */
	Solver_NRcurrent::Solver_NRcurrent(Circuit model) {
		Model = model;

		Sol = Model.get_initial_cx_solution();
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_cx_solution();
		}

		fill_especifyed_values();

		if (!checks())
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
	}

	/* Constructor where the solution is given
	 */
	Solver_NRcurrent::Solver_NRcurrent(Circuit model, cx_solution sol_) {
		Model = model;

		Sol = sol_;

		fill_especifyed_values();

		if (!checks())
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
	}

	/* Destructor
	 */
	Solver_NRcurrent::~Solver_NRcurrent() {
	}


	//check the solvability with this method

	bool Solver_NRcurrent::checks() {
		bool val = true;

		//Only one slack bus allowed
		if (Model.slackBusIndices.size() > 1)
			val = false;

		return val;
	}

	/*Fills the initial power values
	 */
	void Solver_NRcurrent::fill_especifyed_values() {

		uint N = Model.loadBusIndices.size() + Model.generatorBusIndices.size();
		for (uint i = 0; i < Model.buses.size(); i++) {
			if (Model.buses[i].Type == BusType::PQ || Model.buses[i].Type == BusType::PV)
				PQPV.push_back(i);
		}


		Pesp = vec(Model.buses.size());
		Qesp = vec(Model.buses.size());

		for (uint k : PQPV) {
			Pesp(k) = Sol.Pi(k); //P at PQ and PV buses
			Qesp(k) = Sol.Qi(k); //Q at PQ buses
		}

	}

	/* This function returns the calculated increments of y
	 */
	void Solver_NRcurrent::inc_y(vec &x, cx_solution &sol, uint N) {
		if (Debug)
			std::cout << "\nΔY:" << std::endl;


		uint a = 0;
		uint k;
		double Ical_r, Ical_m, incP, incQ, Vk2;

		for (uint idx = 0; idx < N; idx++) { //filas
			k = PQPV[idx];
			//Increment of real current
			Ical_r = 0;
			Ical_m = 0;
			for (uint i = 0; i < Model.buses.size(); i++) {
				Ical_r += Model.G(k, i) * sol.Vr(i) - Model.B(k, i) * sol.Vi(i);
				Ical_m += Model.G(k, i) * sol.Vi(i) + Model.B(k, i) * sol.Vr(i);
			}

			incP = Pesp.coeff(k) - (sol.Vr(k) * Ical_r + sol.Vi(k) * Ical_m);
			incQ = Qesp.coeff(k) - (sol.Vi(k) * Ical_r - sol.Vr(k) * Ical_m);

			Vk2 = sol.Vr(k) * sol.Vr(k) + sol.Vi(k) * sol.Vi(k); //Square voltage

			if (Model.buses[k].Type == BusType::PQ) {
				//increment of imaginary current
				x(a) = (sol.Vi(k) * incP - sol.Vr(k) * incQ) / Vk2; //1

				//increment of real current
				x(a + 1) = (sol.Vr(k) * incP + sol.Vi(k) * incQ) / Vk2; //2

				if (Debug) {
					std::cout << "ΔIm" << k << "\t=\t" << x.coeff(a) << std::endl;
					std::cout << "ΔIr" << k << "\t=\t" << x.coeff(a + 1) << std::endl;
				}
			} else if (Model.buses[k].Type == BusType::PV) {

				//increment of imaginary current
				x(a) = sol.Vi(k) * incP / Vk2; //3

				//increment of real current
				x(a + 1) = sol.Vr(k) * incP / Vk2; //4

				if (Debug) {
					std::cout << "ΔIm*" << k << "\t=\t" << x.coeff(a) << std::endl;
					std::cout << "ΔIr*" << k << "\t=\t" << x.coeff(a + 1) << std::endl;
				}
			}

			a += 2;
		}

	}

	/*Calculates the abcd parameters
	 */
	void Solver_NRcurrent::abcd(uint k, cx_solution &sol, double &a, double &b, double &c, double &d) {
		double V4 = pow(sol.Vr(k), 4) + pow(sol.Vi(k), 4);

		a = (sol.Qi(k) * (sol.Vr(k) * sol.Vr(k) - sol.Vi(k) * sol.Vi(k))
				- 2.0 * sol.Vr(k) * sol.Vi(k) * sol.Pi(k))
				/ V4;

		d = a;

		b = (sol.Pi(k) * (sol.Vr(k) * sol.Vr(k) - sol.Vi(k) * sol.Vi(k))
				+ 2.0 * sol.Vr(k) * sol.Vi(k) * sol.Qi(k))
				/ V4;
		c = -b;
	}

	/*
	 * Calculates the Jacobian
	 */
	void Solver_NRcurrent::Jacobian(mat &J, cx_solution &sol, uint N, bool updating) {

		/* indices: x, y: conceptual bus index
		 *          i, k: real bus indices
		 *          a, b: jacobian indices
		 *
		 */
		uint Nj = 2 * N;


		//cout << "Updating:" << updating << endl;
		if (!updating)
			J.setZero(Nj, Nj);

		//Model.print_buses_state();

		//W
		uint a = 0;
		uint b, k, i;
		double ak, bk, ck, dk, Vk2, B1, B2, G1, G2;
		cx_double ZERO(0.0, 0.0);

		for (uint x = 0; x < N; x++) { //rows
			b = 0;
			k = PQPV[x];
			for (uint y = 0; y < N; y++) { //cols
				i = PQPV[y];

				if (Model.Y.coeff(k, i) != ZERO) {
					if (i == k) { //Diagonal sub-Jacobians

						abcd(k, sol, ak, bk, ck, dk); //always for the diagonal
						B1 = Model.B(k, k) - ak; //B'
						B2 = -1 * Model.B(k, k) - dk; //B''
						G1 = Model.G(k, k) - bk; // G'
						G2 = Model.G(k, k) - ck; //G''

						if (Model.buses[k].Type == BusType::PQ) { //Ykk*
							// always update
							J(a, b) = B1;
							J(a, b + 1) = G1;
							J(a + 1, b) = G2;
							J(a + 1, b + 1) = B2;

							//only to debug
							/*J(a, b) = 1.1;
							J(a, b + 1) = 1.2;
							J(a + 1, b) = 1.3;
							J(a + 1, b + 1) = 1.4;*/
						} else if (Model.buses[k].Type == BusType::PV) {//Ykk**
							// always update

							Vk2 = sol.Vr(k) * sol.Vr(k) - sol.Vi(k) * sol.Vi(k); //Square voltage

							J(a, b) = G1 - B1 * sol.Vi(k) / sol.Vr(k);
							J(a, b + 1) = sol.Vr(k) / Vk2;
							J(a + 1, b) = B2 - G2 * sol.Vi(k) / sol.Vr(k);
							J(a + 1, b + 1) = -1.0 * sol.Vi(k) / Vk2;

							//only to debug
							/*J(a, b) = 2.1;
							J(a, b + 1) = 2.2;
							J(a + 1, b) = 2.3;
							J(a + 1, b + 1) = 2.4;*/
						}

					} else { //Non diagonal sub-Jacobians

						if (Model.buses[k].Type == BusType::PQ
								&& Model.buses[i].Type == BusType::PQ) { //Ykm*
							//does not update
							if (!updating) {
								J(a, b) = Model.B(k, i);
								J(a, b + 1) = Model.G(k, i);
								J(a + 1, b) = Model.G(k, i);
								J(a + 1, b + 1) = -1 * Model.B(k, i);

								//only to debug
								/*J(a, b) = 2;
								J(a, b + 1) = 2;
								J(a + 1, b) = 2;
								J(a + 1, b + 1) = 2;*/
							}
						} else if (Model.buses[i].Type == BusType::PV) { //Ylk**
							// always update
							J(a, b) = Model.G(k, i) - Model.B(k, i) * sol.Vi(k) / sol.Vr(k);
							J(a, b + 1) = 0.0;
							J(a + 1, b) = -1.0 * Model.B(k, i) - Model.G(k, i) * sol.Vi(k) / sol.Vr(k);
							J(a + 1, b + 1) = 0.0;

							//only to debug
							/*J(a, b) = 3;
							J(a, b + 1) = 0;
							J(a + 1, b) = 3;
							J(a + 1, b + 1) = 0;*/
						} else if (Model.buses[k].Type == BusType::PV) {//Ykl** = Ykl*
							//does not update
							if (!updating) {
								J(a, b) = Model.B(k, i);
								J(a, b + 1) = Model.G(k, i);
								J(a + 1, b) = Model.G(k, i);
								J(a + 1, b + 1) = -1 * Model.B(k, i);

								//only to debug
								/*J(a, b) = 2;
								J(a, b + 1) = 2;
								J(a + 1, b) = 2;
								J(a + 1, b + 1) = 2;*/
							}
						}
					}
				}

				b += 2;
			}
			a += 2;
		}
	}

	/*check if the solution converged
	 */
	bool Solver_NRcurrent::converged(vec &X, uint Nj) {
		for (uint i = 0; i < Nj; i++)
			if (abs(X.coeff(i)) > EPS)
				return false;

		return true;
	}

	/* Generates a solution object from a vector X
	 */
	void Solver_NRcurrent::update_solution(cx_solution & sol, vec &x, uint N) {
		if (Debug)
			std::cout << "\nΔX:" << std::endl;
		uint a = 0;
		uint k;
		for (uint idx = 0; idx < N; idx++) { //filas
			k = PQPV[idx];
			if (Model.buses[k].Type == BusType::PQ) {

				//x[a] -> inc Vr
				//x[a+1] -> inc Vi
				sol.V(k) += cx_double(x.coeff(a), x.coeff(a + 1));

				if (Debug) {
					std::cout << "ΔVr*" << k << "\t=\t" << x.coeff(a) << std::endl;
					std::cout << "ΔVm*" << k << "\t=\t" << x.coeff(a + 1) << std::endl;
				}

			} else if (Model.buses[k].Type == BusType::PV) {
				//x[a] -> inc Vi
				//x[a+1] -> inc Q                
				sol.V(k) = cx_double(sol.Vr(k), sol.Vi(k) + x.coeff(a));
				sol.S(k) = cx_double(sol.Pi(k), sol.Qi(k) + x.coeff(a + 1));

				//PQ PV controller
				if (sol.Qi(k) < Model.buses[k].min_q) {

					std::cout << "PV-> PQ Bus " << k << ":: Q=" << sol.Qi(k) << ", qmin:" << Model.buses[k].min_q << ", qmax:" << Model.buses[k].max_q << std::endl;

					Model.buses[k].Type = BusType::PQ;
					sol.S(k) = cx_double(sol.Pi(k), Model.buses[k].min_q);

				} else if (sol.Qi(k) > Model.buses[k].max_q) {

					std::cout << "PV-> PQ Bus " << k << ":: Q=" << sol.Qi(k) << ", qmin:" << Model.buses[k].min_q << ", qmax:" << Model.buses[k].max_q << std::endl;

					Model.buses[k].Type = BusType::PQ;
					sol.S(k) = cx_double(sol.Pi(k), Model.buses[k].max_q);
				}

				if (Debug) {
					std::cout << "ΔVm*" << k << "\t=\t" << x.coeff(a) << std::endl;
					std::cout << "ΔQ*" << k << "\t=\t" << x.coeff(a + 1) << std::endl;
				}
			}

			a += 2;
		}
	}

	/*calculate the slack bus power 
	 */
	void Solver_NRcurrent::calculate_slack_power() {
		for (uint k : Model.slackBusIndices) {
			cx_double I(0.0, 0.0);
			for (uint j = 0; j < Model.buses.size(); j++) {
				I += Model.Y.coeff(k, j) * Sol.V(j);
			}
			I = Sol.V(k) * conj(I); //now this is the power
			Sol.S(k) = I;
		}
	}

	/*Solves the grid
	 */
	Solver_State Solver_NRcurrent::solve() {

		Sol.print("Initial solution:");

		uint npv = Model.generatorBusIndices.size();
		uint npq = Model.loadBusIndices.size();
		uint N = npq + npv;
		uint Nj = 2 * N;

		//Declare vectors
		mat J(Nj, Nj);
		vec incX(Nj);
		vec incY(Nj);

		//generate vector of mismatches
		inc_y(incY, Sol, N);

		//check convergence
		bool converged_ = converged(incY, Nj);
		bool updating = false;

		Iterations = 0;
		while (!converged_ && Iterations < maxIterations) {

			//prints
			std::cout << "\n\nIter: " << Iterations << "\n" << std::endl;

			//Update Jacobian
			Jacobian(J, Sol, N, updating);
			Eigen::FullPivLU<mat>LU(J); //Full pivot LU 
			//cout << "\nJ:\n" << J << endl;
			if (!updating)
				updating = true; //the Jacobian has been created, now only update it.

			//Solve the linear system to get the increments 
			incX = LU.solve(incY);

			//Update the solution object with the solution vector: recalculate
			// the voltages with the increments
			update_solution(Sol, incX, N);
			Sol.print("solution:");

			//Update mismatches
			inc_y(incY, Sol, N); //generate vector of mismatches

			//check convergence
			converged_ = converged(incY, Nj);

			Iterations++;
		}

		//Sol.print("NR current solution:");


		if (converged_) {
			calculate_slack_power();
			Model.set_solution(Sol);
			return Solver_State::Converged;
		} else {
			Sol.print("Solution so far:");
			return Solver_State::Not_Converged;
		}
	}

	/*
	 */
	void Solver_NRcurrent::update_solution_power_from_circuit() {

	}
}
