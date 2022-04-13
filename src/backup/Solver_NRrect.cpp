/*
 * File:   Solver_NRpolar.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Solver_NRrect.h"
#include "Circuit.h"

namespace fPotencia {

	/*
	 * constructor
	 */
	Solver_NRrect::Solver_NRrect(Circuit model) {

		Model = model;

		Sol = Model.get_initial_cx_solution();
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_cx_solution();
		}

		BUSES.reserve(Model.loadBusIndices.size() + Model.generatorBusIndices.size()); // preallocate memory
		BUSES.insert(BUSES.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		BUSES.insert(BUSES.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		PQPV.reserve(2 * Model.loadBusIndices.size() + Model.generatorBusIndices.size()); // preallocate memory
		PQPV.insert(PQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPV.insert(PQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		PQPVPQPV.reserve(2 * (Model.loadBusIndices.size() + Model.generatorBusIndices.size())); // preallocate memory
		PQPVPQPV.insert(PQPVPQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		fill_especifyed_values();

		if (!checks())
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
	}

	/*
	 * constructor, given previous solution
	 */
	Solver_NRrect::Solver_NRrect(Circuit model, cx_solution sol_) {

		Model = model;

		Sol = sol_;

		BUSES.reserve(Model.loadBusIndices.size() + Model.generatorBusIndices.size()); // preallocate memory
		BUSES.insert(BUSES.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		BUSES.insert(BUSES.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		PQPV.reserve(2 * Model.loadBusIndices.size() + Model.generatorBusIndices.size()); // preallocate memory
		PQPV.insert(PQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPV.insert(PQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		PQPVPQPV.reserve(2 * (Model.loadBusIndices.size() + Model.generatorBusIndices.size())); // preallocate memory
		PQPVPQPV.insert(PQPVPQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		PQPVPQPV.insert(PQPVPQPV.end(), Model.generatorBusIndices.begin(), Model.generatorBusIndices.end());

		fill_especifyed_values();

		if (!checks())
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
	}

	/*
	 * destructor
	 */
	Solver_NRrect::~Solver_NRrect() {
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * checks wether if the grid can be solved with this method
	 */
	bool Solver_NRrect::checks() {
		bool val = true;

		//Only one slack bus allowed
		if (Model.slackBusIndices.size() > 1)
			val = false;

		return val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the slack bus power
	 */
	void Solver_NRrect::calculate_slack_power() {
		for (uint k : Model.slackBusIndices) {
			cx_double I(0.0, 0.0);
			for (uint j = 0; j < Model.buses.size(); j++)
				I += Model.Y.coeff(k, j) * Sol.V(j);

			Sol.S(k) = Sol.V(k) * conj(I); //now this is the power
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the reactive power of the bus k (usefull for PV uses)
	 */
	void Solver_NRrect::calculate_Q(uint npq, uint npv) {
		double val;
		uint k;
		for (uint i = npq - 1; i < npq + npv; i++) {
			k = PQPV[i];
			Sol.S(k) = cx_double(Sol.S.coeff(k).real(), Q(k));
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the active power at a bus
	 */
	double Solver_NRrect::P(uint i) {
		double val = 0.0;
		for (uint k = 0; k < Model.buses.size(); k++) {
			val += Sol.Vr(i)* (Model.G(i, k) * Sol.Vr(k) - Model.B(i, k) * Sol.Vi(k)) + Sol.Vi(i)* (Model.G(i, k) * Sol.Vi(k) + Model.B(i, k) * Sol.Vr(k));
		}
		//return Model.G(i, i) * Sol.Vr(i) * Sol.Vr(i) + Model.G(i, i) * Sol.Vi(i) * Sol.Vi(i) + val;
		return val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the reactive power at a bus
	 */
	double Solver_NRrect::Q(uint i) {
		double val = 0.0;
		for (uint k = 0; k < Model.buses.size(); k++) {
			val += Sol.Vi(i)* (Model.G(i, k) * Sol.Vr(k) - Model.B(i, k) * Sol.Vi(k)) - Sol.Vr(i)* (Model.G(i, k) * Sol.Vi(k) + Model.B(i, k) * Sol.Vr(k));
		}
		//return Model.B(i, i) * Sol.Vr(i) * Sol.Vr(i) + Model.B(i, i) * Sol.Vi(i) * Sol.Vi(i) + val;
		return val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Real Current injection at bus i
	 */
	double Solver_NRrect::c(uint i) {
		double val = 0.0;
		for (uint k = 0; k < Model.buses.size(); k++)
			if (k != i)
				val += Model.G(i, k) * Sol.Vr(k) - Model.B(i, k) * Sol.Vi(k);

		return Model.G(i, i) * Sol.Vr(i) - Model.B(i, i) * Sol.Vi(i) + val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Real Current injection at bus i
	 */
	double Solver_NRrect::d(uint i) {
		double val = 0.0;
		for (uint k = 0; k < Model.buses.size(); k++)
			if (k != i)
				val += Model.G(i, k) * Sol.Vi(k) + Model.B(i, k) * Sol.Vr(k);

		return Model.G(i, i) * Sol.Vi(i) + Model.B(i, i) * Sol.Vr(i) + val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * This function gets the especified values correctly
	 */
	void Solver_NRrect::fill_especifyed_values() {

		Pesp = sp_vec(Model.buses.size());
		Qesp = sp_vec(Model.buses.size());
		V2esp = sp_vec(Model.buses.size());

		for (uint k : PQPV)
			Pesp.coeffRef(k) = Sol.S.coeff(k).real(); //P at PQ and PV buses

		for (uint k : Model.loadBusIndices)
			Qesp.coeffRef(k) = Sol.S.coeff(k).imag(); //Q at PQ buses

		for (uint k : Model.generatorBusIndices)
			V2esp.coeffRef(k) = pow(abs(Sol.V.coeff(k)), 2.0); //V at PV buses ^ 2
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the jacobian of the circuit
	 */
	void Solver_NRrect::Jacobian(mat &J, uint npq, uint npv) {
		//matrix(rows, cols)
		uint npqpv = npq + npv;
		uint i, k;
		mat J1(npqpv, npqpv);
		mat J2(npqpv, npqpv);
		mat J3(npq, npqpv);
		mat J4(npq, npqpv);
		mat J5(npv, npqpv);
		mat J6(npv, npqpv);
		double val;

		//cout << "J1 and 2" << endl;
		//J1 and J2
		for (uint a = 0; a < npqpv; a++) {
			i = PQPV[a];
			for (uint b = 0; b < npqpv; b++) {
				k = PQPV[b];
				if (i == k) {
					J1(a, b) = Model.G(i, i) * Sol.Vr(i) + Model.B(i, i) * Sol.Vi(i) + c(i);
					J2(a, b) = -Model.B(i, i) * Sol.Vr(i) + Model.G(i, i) * Sol.Vi(i) + d(i);
				} else {
					J1(a, b) = Model.G(i, k) * Sol.Vr(i) + Model.B(i, k) * Sol.Vi(i);
					J2(a, b) = Model.G(i, k) * Sol.Vi(i) - Model.B(i, k) * Sol.Vr(i);
				}
			}
		}

		//cout << "J3 and 4" << endl;
		//J3 and J4
		for (uint a = 0; a < npq; a++) {
			i = Model.loadBusIndices[a];
			for (uint b = 0; b < npqpv; b++) {
				k = PQPV[b];
				if (i == k) {
					J3(a, b) = Model.G(i, i) * Sol.Vi(i) - Model.B(i, i) * Sol.Vr(i) - d(i);
					J4(a, b) = -Model.G(i, i) * Sol.Vr(i) - Model.B(i, i) * Sol.Vi(i) + c(i);
				} else {
					J3(a, b) = Model.G(i, k) * Sol.Vi(i) - Model.B(i, k) * Sol.Vr(i);
					J4(a, b) = -Model.G(i, k) * Sol.Vr(i) - Model.B(i, k) * Sol.Vi(i);
				}
			}
		}

		//cout << "J5 and 6" << endl;
		//J5 and J6
		for (uint a = 0; a < npv; a++) {
			i = Model.generatorBusIndices[a];
			J5(a, a) = 2.0 * Sol.Vr(i);
			J6(a, a) = 2.0 * Sol.Vi(i);
		}

		J << J1, J2,
				J3, J4,
				J5, J6;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the power increments
	 */
	void Solver_NRrect::get_mismatches(vec &inc, uint npq, uint npv) {

		uint npqpv, npqpvpq, k;
		npqpv = npq + npv;
		npqpvpq = npq + npv + npq;

		for (uint i = 0; i < 2 * npqpv; i++) {
			k = PQPVPQPV[i];
			if (i < npqpv) {
				//cout << i << " incP:" << Pesp.coeff(k) << " - " << P(k) << endl;
				inc(i) = Pesp.coeff(k) - P(k); //P at PQ and PV buses
			} else if (i >= npqpv && i < npqpvpq) {
				//cout << i << " incQ:" << Qesp.coeff(k) << " - " << Q(k) << endl;
				inc(i) = Qesp.coeff(k) - Q(k); //Q at PQ buses
			} else {
				//cout << i << " incV:" << V2esp.coeff(k) << " - " << pow(abs(Sol.V.coeff(k)), 2.0) << endl;
				inc(i) = V2esp.coeff(k) - pow(abs(Sol.V.coeff(k)), 2.0); //V at PV buses ^ 2
			}
		}

	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Check the convergence
	 * Error gets defined as max(abs(inc))
	 */
	bool Solver_NRrect::converged(vec inc, uint npqpvpq, double &error) {
		double err;
		error = 0.0;        
		for (uint k = 0; k < npqpvpq; k++) {
			err = abs(inc.coeff(k));
			if (err > error)
				error = err;
		}

		if (error > EPS)
			return false;

		return true;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Check the convergence
	 */
	void Solver_NRrect::update_solution(vec X, uint npq, uint npv) {

		uint npqpv, npqpvpq, k;
		npqpv = npq + npv;
		npqpvpq = npq + npv + npq;

		double d;
		for (uint i = 0; i < npqpv; i++) {
			k = PQPV[i];
			// update V (complex) at bus k
			Sol.V(k) += cx_double(X.coeff(i), X.coeff(i + npqpv));

			if (i >= npq) { //Control the PV voltages
				d = atan(Sol.Vi(k) / Sol.Vr(k));
				Sol.V(k) = cx_double(V2esp.coeff(k) * cos(d), V2esp.coeff(k) * sin(d)); //update V at PV buses ^ 2
			}
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Run the solve process
	 */
	Solver_State Solver_NRrect::solve() {

		uint npq = Model.loadBusIndices.size();
		uint npv = Model.generatorBusIndices.size();
		uint n2pqpv = 2 * (npq + npv); //size of the arrays

		//System : J*X = K
		mat J(n2pqpv, n2pqpv);
		vec increments(n2pqpv);
		vec new_voltages(n2pqpv);

		bool conv;
		double error;
		double alpha = 1.0;
		double error_old;
		Iterations = 0;

		//first calculations
		get_mismatches(increments, npq, npv);
		//First convergence check initializes the old error
		conv = converged(increments, n2pqpv, error_old);

		//std::cout << "K:\n" << K << std::endl;
		//std::cout << "Converged: " << conv << std::endl;

		while (!conv && Iterations < Max_Iter) {
			//std::cout << "-----------------------Iter: " << Iterations << std::endl;
			//std::cout << "K:\n" << K << std::endl;
			//Calculate the jacobian
			Jacobian(J, npq, npv);
			//std::cout << "J:\n" << J << std::endl;

			Eigen::FullPivLU<mat>lu(J); //Full pivot LU
			new_voltages = lu.solve(alpha * increments);
			//std::cout << "X:\n" << X << std::endl;

			//update the solution
			update_solution(new_voltages, npq, npv);

			calculate_Q(npq, npv); //Calculate the reactive power for the PV buses

			//Calculate the increment of power for the new iteration
			get_mismatches(increments, npq, npv);
			//std::cout << "K:\n" << K << std::endl;

			//Check convergency
			conv = converged(increments, n2pqpv, error);

			//Robustness
			cout << alpha << " : " << error << " <> " << error_old << endl;
			if (error > error_old) {
				alpha *= 0.5;
				cout << "alpha:" << alpha << endl;
				if (alpha < EPS) {
					return Solver_State::Not_Converged;
				}
			} else {
				alpha = 1.0;
			}


			//Sol.print("Soltution:");
			error_old = error;
			Iterations++;
		}

		if (!conv || Iterations == Max_Iter)
			return Solver_State::Not_Converged;
		else {
			//std::cout << "Converged in " << Iterations << " iterations." << std::endl;
			calculate_slack_power();
			Model.set_solution(Sol);
			//Sol.print("Final Soltution:");
			return Solver_State::Converged;
		}
	}

	/*This function updates the solver solution object power values using the
	 * circuit's own solution power values. this is specially usefull when updating
	 * the circuit power values while keeping the previous voltage solution
	 */
	void Solver_NRrect::update_solution_power_from_circuit() {
		Sol.S = Model.get_initial_cx_solution().S;
		fill_especifyed_values();
	}

}