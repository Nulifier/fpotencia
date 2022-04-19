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
	NRcurrentSolver::NRcurrentSolver(Circuit& model): Model(model) {
		Sol = Model.get_initial_cx_solution();
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_cx_solution();
		}

		fill_specified_values();

		if (!canSolve()) {
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
		}
	}

	NRcurrentSolver::NRcurrentSolver(Circuit& model, cx_solution sol_): Model(model) {
		Sol = sol_;

		fill_specified_values();

		if (!canSolve()) {
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
		}
	}

	bool NRcurrentSolver::canSolve() const {
		return Model.slackBusIndices.size() <= 1;
	}

	void NRcurrentSolver::fill_specified_values() {

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

	void NRcurrentSolver::inc_y(vec &x, cx_solution &sol, uint N) {
		uint a = 0;
		uint k;
		double Ical_r, Ical_m, incP, incQ, Vk2;

		for (uint idx = 0; idx < N; idx++) { //filas
			k = PQPV[idx];
			//Increment of real current
			Ical_r = 0;
			Ical_m = 0;
			for (uint i = 0; i < Model.buses.size(); ++i) {
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
			}
			else if (Model.buses[k].Type == BusType::PV) {
				//increment of imaginary current
				x(a) = sol.Vi(k) * incP / Vk2; //3

				//increment of real current
				x(a + 1) = sol.Vr(k) * incP / Vk2; //4
			}

			a += 2;
		}
	}

	void NRcurrentSolver::abcd(uint k, cx_solution &sol, double &a, double &b, double &c, double &d) {
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

	void NRcurrentSolver::Jacobian(mat &J, cx_solution &sol, uint N, bool updating) {
		/* indices: x, y: conceptual bus index
		 *          i, k: real bus indices
		 *          a, b: jacobian indices
		 */
		uint Nj = 2 * N;

		if (!updating) {
			J.setZero(Nj, Nj);
		}

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
						}
						else if (Model.buses[k].Type == BusType::PV) {//Ykk**
							// always update

							Vk2 = sol.Vr(k) * sol.Vr(k) - sol.Vi(k) * sol.Vi(k); //Square voltage

							J(a, b) = G1 - B1 * sol.Vi(k) / sol.Vr(k);
							J(a, b + 1) = sol.Vr(k) / Vk2;
							J(a + 1, b) = B2 - G2 * sol.Vi(k) / sol.Vr(k);
							J(a + 1, b + 1) = -1.0 * sol.Vi(k) / Vk2;
						}
					}
					else { //Non diagonal sub-Jacobians
						if (Model.buses[k].Type == BusType::PQ
								&& Model.buses[i].Type == BusType::PQ) { //Ykm*
							//does not update
							if (!updating) {
								J(a, b) = Model.B(k, i);
								J(a, b + 1) = Model.G(k, i);
								J(a + 1, b) = Model.G(k, i);
								J(a + 1, b + 1) = -1 * Model.B(k, i);
							}
						}
						else if (Model.buses[i].Type == BusType::PV) { //Ylk**
							// always update
							J(a, b) = Model.G(k, i) - Model.B(k, i) * sol.Vi(k) / sol.Vr(k);
							J(a, b + 1) = 0.0;
							J(a + 1, b) = -1.0 * Model.B(k, i) - Model.G(k, i) * sol.Vi(k) / sol.Vr(k);
							J(a + 1, b + 1) = 0.0;
						}
						else if (Model.buses[k].Type == BusType::PV) {//Ykl** = Ykl*
							//does not update
							if (!updating) {
								J(a, b) = Model.B(k, i);
								J(a, b + 1) = Model.G(k, i);
								J(a + 1, b) = Model.G(k, i);
								J(a + 1, b + 1) = -1 * Model.B(k, i);
							}
						}
					}
				}

				b += 2;
			}
			a += 2;
		}
	}

	bool NRcurrentSolver::converged(vec &X, uint Nj) const {
		for (uint i = 0; i < Nj; i++) {
			if (abs(X.coeff(i)) > tolerance) {
				return false;
			}
		}

		return true;
	}

	void NRcurrentSolver::update_solution(cx_solution & sol, vec &x, uint N) {
		uint a = 0;
		uint k;
		for (uint idx = 0; idx < N; idx++) { //filas
			k = PQPV[idx];
			if (Model.buses[k].Type == BusType::PQ) {

				//x[a] -> inc Vr
				//x[a+1] -> inc Vi
				sol.V(k) += cx_double(x.coeff(a), x.coeff(a + 1));

			}
			else if (Model.buses[k].Type == BusType::PV) {
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
			}

			a += 2;
		}
	}

	void NRcurrentSolver::calculate_slack_power() {
		for (unsigned int k : Model.slackBusIndices) {
			cx_double I(0.0, 0.0);
			for (uint j = 0; j < Model.buses.size(); j++) {
				I += Model.Y.coeff(k, j) * Sol.V(j);
			}
			I = Sol.V(k) * conj(I); //now this is the power
			Sol.S(k) = I;
		}
	}

	Solver::Result NRcurrentSolver::solve() {
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
		bool didConverge = converged(incY, Nj);
		bool updating = false;

		for (unsigned int i = 0; i < maxIterations && !didConverge; ++i) {
			//Update Jacobian
			Jacobian(J, Sol, N, updating);
			
			Eigen::FullPivLU<mat> LU(J); //Full pivot LU 
			
			if (!updating) {
				updating = true; //the Jacobian has been created, now only update it.
			}

			//Solve the linear system to get the increments 
			incX = LU.solve(incY);

			//Update the solution object with the solution vector: recalculate
			// the voltages with the increments
			update_solution(Sol, incX, N);

			//Update mismatches
			inc_y(incY, Sol, N); //generate vector of mismatches

			//check convergence
			didConverge = converged(incY, Nj);
		}

		if (didConverge) {
			calculate_slack_power();
			Model.set_solution(Sol);
			return Solver::Solved;
		} else {
			return Solver::NotSolved;
		}
	}

	void NRcurrentSolver::update_solution_power_from_circuit() {

	}
}
