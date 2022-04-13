/*
 * File:   Solver_ZBlockSubs.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de noviembre de 2014, 23:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Circuit.h"
#include "Bus.h"
#include "Solver_ZBusGS.h"

namespace fPotencia {

	/*//////////////////////////////////////////////////////////////////////////
	 * solver object constructor given the circuit model
	 */
	Solver_ZBusGS::Solver_ZBusGS(Circuit model) {

		Model = model;

		Sol = Model.get_initial_cx_solution();
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_cx_solution();
		}

		BUSES.reserve(Model.loadBusIndices.size()); // preallocate memory
		BUSES.insert(BUSES.end(), Model.loadBusIndices.begin(), Model.loadBusIndices.end());
		//BUSES.insert(BUSES.end(), Model.PV_list.begin(), Model.PV_list.end());

		E0.resize(Model.buses.size());
		E.resize(Model.buses.size());
		I.resize(Model.buses.size());
		C.resize(Model.buses.size());


		for (uint i = 0; i < Model.buses.size(); i++) {
			E0(i) = Sol.V.coeff(i);
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * solver object destructor
	 */
	Solver_ZBusGS::~Solver_ZBusGS() {
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * checks wether if the grid can be solved
	 */
	bool Solver_ZBusGS::checks() {
		bool val = true;

		//The Zbus only works with PQ buses
		if (Model.generatorBusIndices.size() > 0)
			val = false;

		//Only one slack bus allowed
		if (Model.slackBusIndices.size() > 1)
			val = false;

		return val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculates the current for the bus k
	 */
	cx_double Solver_ZBusGS::calculate_I_k(uint k, cx_vec * V) {
		return (conj(Sol.S.coeff(k)) / conj((*V).coeff(k)));
		//return (conj(Sol.S.coeff(k)) / conj((*V).coeff(k))) -Model.Y.coeff(k, k) * (*V).coeff(k);
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculates the current for all the bus bars
	 */
	void Solver_ZBusGS::calculate_I(cx_vec * V) {
		for (int k : BUSES) {
			I(k) = calculate_I_k(k, V);
			//cout << "I_" << k << ": S*/I*:" << conj(S_esp(k)) << "/" << conj(Sol.V(k)) << "=" << I.coeff(k) << endl;
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculates the C constant for all the bus bars
	 */
	void Solver_ZBusGS::calculate_C(cx_vec * V) {

		uint s = Model.slackBusIndices[0]; //slack bus index
		cx_double V1 = (*V).coeff(s); //slack bus voltage

		for (int k : BUSES) {
			C(k) = cx_double(0.0, 0.0);
			for (int j : BUSES) {
				C(k) += Model.Zred.coeff(k, j) * Model.Y.coeff(j, s) * V1;
			}
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculates the new voltage for the busbar k
	 */
	cx_double Solver_ZBusGS::calculate_V_new_k(uint i, bool pv) {
		cx_double v(0.0, 0.0);

		for (uint j : BUSES) {
			v += Model.Zred.coeff(j, i) * I.coeff(i);
			//cout << "V_" << k << "{" << k << "," << j << "}+=" << Model.Zred.coeff(k, j) << "*" << I[j] << "-" << C[k] << endl;
		}

		v -= C.coeff(i);

		if (pv)
			/*
			 * This is the voltage correction for the PV buses based on the
			 * connected generator voltage set point
			 */
			v = v * Model.buses[i].v_set_point / abs(v);

		return v;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Checks if the element k of the solution converged
	 */
	bool Solver_ZBusGS::converged_k(uint k) {
		return ( abs(E0.coeff(k) - E.coeff(k)) <= EPS);
		//too simple convergence criteria
		/*
				cx_double sum(0.0, 0.0);

				for (int j = 0; j < Model.buses.size(); j++)
					sum += Model.Y.coeff(k, j) * Sol.V.coeff(j);

				cx_double cdelta = Sol.S.coeff(k) - conj(Sol.V.coeff(k) * sum);
				double delta = abs(cdelta);
				return (abs(cdelta) <= EPS);
		 */
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Checks if all the solution converged
	 */
	bool Solver_ZBusGS::converged() {

		for (uint k : BUSES) {
			//cout << "Convergence:" << abs(newSol.V[k]) << " - " << abs(Sol.V[k]) << " = " << abs(abs(newSol.V[k]) - abs(Sol.V[k])) << endl;
			if (!converged_k(k)) {
				//cout << "converged false" << endl;
				return false;
			}
		}
		cout << "converged true" << endl;
		return true;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculate the reactive power of the bus k (usefull for PV uses)
	 */
	double Solver_ZBusGS::calculate_Q_k(uint k) {

		double val = 0.0;

		for (uint j = 0; j < Model.buses.size(); j++) {
			val += Sol.Vi(j)*(Model.G(k, j) * Sol.Vr(j) - Model.B(k, j) * Sol.Vi(j))
					- Sol.Vr(j)*(Model.G(k, j) * Sol.Vi(j) + Model.B(k, j) * Sol.Vr(j));
		}

		return val;
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Calculates the power for the slack bus
	 */
	void Solver_ZBusGS::calculate_slack_power() {
		//VD (or slack) bus
		cx_double I(0.0, 0.0);
		for (int k : Model.slackBusIndices) {
			for (uint j = 0; j < Model.buses.size(); j++) {
				I += Model.Y.coeff(k, j) * Sol.V.coeff(j);
			}
			Sol.S(k) = Sol.V.coeff(k) * conj(I);
		}
	}

	/*//////////////////////////////////////////////////////////////////////////
	 * Solve using the Z-bus method with the Gauss-Seidel algorithm
	 */
	Solver_State Solver_ZBusGS::solve() {

		if (!checks()) {
			std::cout << "\nThe grid cannot be solved using Zbus Gauss Seidel\n" << std::endl;
			return Solver_State::Not_Solvable_with_Method;
		}

		std::cout << "E0:\n" << E0 << std::endl;
		calculate_C(&E0);
		calculate_I(&E0);
		std::cout << "I:\n" << I << std::endl;
		std::cout << "C:\n" << C << std::endl;

		bool has_converged = false;

		Iterations = 0;

		int k;
		double P, Q;

		while (!has_converged && Iterations < Max_Iter) {

			has_converged = true;

			for (uint k : BUSES) {
				//Calculate the voltage
				E(k) = calculate_V_new_k(k, false);
				I(k) = calculate_I_k(k, &E);
				
				//Check the convergence
				if (!converged_k(k))
					has_converged = false;
			}

			for (uint i = 0; i < Model.buses.size(); i++) {
				E0(i) = E.coeff(i);
			}

			std::cout << "Iter:" << Iterations << std::endl;
			std::cout << "E:\n" << E << std::endl;
			std::cout << "E0:\n" << E0 << std::endl;
			std::cout << "I:\n" << I << std::endl;
			Iterations++;
		}

		for (uint i = 0; i < Model.buses.size(); i++) {
			Sol.V(i) = E.coeff(i);
		}

		cout << "iterations:" << Iterations << endl;

		if (has_converged) {
			calculate_slack_power(); //calculate the power for the slack buses
			Model.set_solution(Sol); //sets the solution to the circuit and calculates the flows
			Sol.print("Solution:");
			return Solver_State::Converged;
		} else
			return Solver_State::Not_Converged;
	}


}