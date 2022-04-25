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

#include "Solver_NRpolar.h"
#include "Circuit.h"

namespace fPotencia {
	NRpolarSolver::NRpolarSolver(Circuit& model):
			Model(model),
			Sol(Model.get_initial_solution())
	{
		if (!Sol.initialized) {
			Model.compile(false);
			Sol = Model.get_initial_solution();
		}

		PQPV.reserve(
				2 * Model.loadBusIndices.size()
					+ Model.generatorBusIndices.size());
		PQPV.insert(
				PQPV.end(),
				Model.loadBusIndices.begin(),
				Model.loadBusIndices.end());
		PQPV.insert(
				PQPV.end(),
				Model.generatorBusIndices.begin(),
				Model.generatorBusIndices.end());

		LastPQ.reserve(Model.loadBusIndices.size());
		LastPQ.insert(
				LastPQ.end(),
				Model.loadBusIndices.begin(),
				Model.loadBusIndices.end());

		LastPV.reserve(Model.generatorBusIndices.size());
		LastPV.insert(
				LastPV.end(),
				Model.generatorBusIndices.begin(),
				Model.generatorBusIndices.end());

		Pesp = Sol.P;
		Qesp = Sol.Q;
		
		if (!canSolve()) {
			throw std::invalid_argument("The circuit failed the solver compatibility test.");
		}
	}

	NRpolarSolver::NRpolarSolver(Circuit& model, const solution& sol_):
		NRpolarSolver(model)
	{
		Sol = sol_;
	}

	bool NRpolarSolver::canSolve() const {
		return Model.slackBusIndices.size() <= 1;
	}

	void NRpolarSolver::calculate_slack_power() {        
		for (unsigned int k: Model.slackBusIndices) {
			cx_double I(0.0, 0.0);
			for (uint j = 0; j < Model.buses.size(); j++) {
				I += Model.Y.coeff(k, j) * Sol.V(j);
			}
			I = Sol.V(k) * conj(I); //now this is the power
			Sol.P(k) = I.real();
			Sol.Q(k) = I.imag();
		}
	}

	void NRpolarSolver::calculate_Q(uint npq, uint npv) {
		double val;
		uint k;
		for (uint i = npq - 1; i < npq + npv; i++) {
			k = PQPV[i];
			val = Q(k);
			Sol.Q(k) = val;
		}
	}

	double NRpolarSolver::P(uint k) {
		double val = 0.0;
		for (uint j = 0; j < Model.buses.size(); j++) {
			val += Sol.Vmag.coeff(j)
					*(Model.G(k, j) * cos(Sol.Varg.coeff(k) - Sol.Varg.coeff(j))
					+ Model.B(k, j) * sin(Sol.Varg.coeff(k) - Sol.Varg.coeff(j)));
		}
		return Sol.Vmag.coeff(k) * val;
	}

	double NRpolarSolver::Q(uint k) {
		double val = 0.0;
		for (uint j = 0; j < Model.buses.size(); j++) {
			val += Sol.Vmag.coeff(j)
					*(Model.G(k, j) * sin(Sol.Varg.coeff(k) - Sol.Varg.coeff(j))
					- Model.B(k, j) * cos(Sol.Varg.coeff(k) - Sol.Varg.coeff(j)));
		}
		return Sol.Vmag.coeff(k) * val;
	}

	void NRpolarSolver::correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X) {
		// Find the PV buses that violate the reative power limit
		std::vector<int> lst;
		uint k;
		for (uint i = npq; i < npq + npv; i++) {
			k = PQPV[i];
			if (Model.buses[k].max_q != 0.0) {
				if (Sol.Q(k) >= Model.buses[k].max_q) {
					Sol.Q(k) = Model.buses[k].max_q; //truncate Q to the limit
					lst.push_back(k); //add the PV bus to the list to be treated as a PQ bus
				}
			}
			else {
				//std::cout << "Probably invalid reactie power value at bus " << Model.buses[k].Name << std::endl;
			}
		}

		// Change the lists and arrays size to accomodate the new situation
		npq += lst.size();
		npv -= lst.size();
		uint npqpvpq = 2 * npq + npv; //size of the arrays

		// Resize the linear system, since if npq and npv vary their size vary as well
		J = mat(npqpvpq, npqpvpq);
		K = vec(npqpvpq);
		X = vec(npqpvpq);

		// add the PV buses indices from lst to the PQ list
		LastPQ.insert(LastPQ.end(), lst.begin(), lst.end());

		// Remove the same lst indices fromthe PV list
		for (uint k : lst)
			LastPV.erase(std::remove(LastPV.begin(), LastPV.end(), k), LastPV.end());

		PQPV.clear();
		PQPV.reserve(LastPQ.size() + LastPV.size()); // preallocate memory
		PQPV.insert(PQPV.end(), LastPQ.begin(), LastPQ.end());
		PQPV.insert(PQPV.end(), LastPV.begin(), LastPV.end());
	}

	void NRpolarSolver::Jacobian(mat &J, vec &V, vec &D, uint npq, uint npv) {
		//matrix(rows, cols)
		uint npqpv = npq + npv;
		double val;
		uint k, j;

		J.setZero();

		//J1 
		for (uint a = 0; a < npqpv; a++) { //rows
			k = PQPV[a];
			//diagonal
			J(a, a) = -Q(k) - Model.B(k, k) * V.coeff(k) * V.coeff(k);

			//non diagonal elements
			for (uint b = 0; b < npqpv; b++) {
				if (b != a) {
					j = PQPV[b];
					val = V.coeff(k) * V.coeff(j)
							*(Model.G(k, j) * sin(D.coeff(k) - D.coeff(j))
							- Model.B(k, j) * cos(D.coeff(k) - D.coeff(j)));
					J(a, b) = val;
				}
			}
		}

		//J2
		unsigned int da = 0;
		unsigned int db = npqpv;
		for (uint a = 0; a < npqpv; a++) { //rows
			k = PQPV[a];
			//diagonal
			if (a < npq)
				J(a + da, a + db) = P(k) + Model.G(k, k) * V.coeff(k) * V.coeff(k);

			//non diagonal elements
			for (uint b = 0; b < npq; b++) {
				if (b != a) {
					j = PQPV[b];
					val = V.coeff(k) * V.coeff(j)
							*(Model.G(k, j) * cos(Sol.Varg.coeff(k) - Sol.Varg.coeff(j))
							+ Model.B(k, j) * sin(Sol.Varg.coeff(k) - Sol.Varg.coeff(j)));
					J(a + da, b + db) = val;
				}
			}
		}


		//J3
		da = npqpv;
		db = 0;
		for (uint a = 0; a < npq; a++) { //rows
			k = PQPV[a];
			//diagonal
			J(a + da, a + db) = P(k) - Model.G(k, k) * V.coeff(k) * V.coeff(k);

			//non diagonal elements
			for (uint b = 0; b < npqpv; b++) {
				if (b != a) {
					j = PQPV[b];
					val = V.coeff(k) * V.coeff(j)
							*(Model.G(k, j) * cos(D.coeff(k) - D.coeff(j))
							+ Model.B(k, j) * sin(D.coeff(k) - D.coeff(j)));
					J(a + da, b + db) = -val;
				}
			}
		}

		//J4
		da = npqpv;
		db = npqpv;
		for (uint a = 0; a < npq; a++) { //rows
			k = PQPV[a];
			//diagonal
			J(a + da, a + db) = Q(k) - Model.B(k, k) * V.coeff(k) * V.coeff(k);

			//non diagonal elements
			for (uint b = 0; b < npq; b++) {
				if (b != a) {
					j = PQPV[b];
					val = V.coeff(k) * V.coeff(j)
							*(Model.G(k, j) * sin(D.coeff(k) - D.coeff(j))
							- Model.B(k, j) * cos(D.coeff(k) - D.coeff(j)));
					if (val != 0.0) {
						J(a + da, b + db) = val;
					}
				}
			}
		}


	}
	
	
	/*
	 
	 def mu(Ybus, J, F, dV, dx, pvpq, pq):
	"""
	Calculate the Iwamoto acceleration parameter as described in:
	"A Load Flow Calculation Method for Ill-Conditioned Power Systems" by Iwamoto, S. and Tamura, Y.
	Args:
		Ybus: Admittance matrix
		J: Jacobian matrix
		F: mismatch vector
		dV: voltage increment (in complex form)
		dx: solution vector as calculated dx = solve(J, F)
		pvpq: array of the pq and pv indices
		pq: array of the pq indices

	Returns:
		the Iwamoto's optimal multiplier for ill conditioned systems
	"""
	# evaluate the Jacobian of the voltage derivative
	dS_dVm, dS_dVa = dSbus_dV(Ybus, dV)  # compute the derivatives

	J11 = dS_dVa[array([pvpq]).T, pvpq].real
	J12 = dS_dVm[array([pvpq]).T, pq].real
	J21 = dS_dVa[array([pq]).T, pvpq].imag
	J22 = dS_dVm[array([pq]).T, pq].imag

	# theoretically this is the second derivative matrix
	# since the Jacobian has been calculated with dV instead of V
	J2 = vstack([
			hstack([J11, J12]),
			hstack([J21, J22])
			], format="csr")

	a = F
	b = J * dx
	c = 0.5 * dx * J2 * dx

	g0 = -a.dot(b)
	g1 = b.dot(b) + 2 * a.dot(c)
	g2 = -3.0 * b.dot(c)
	g3 = 2.0 * c.dot(c)

	roots = np.roots([g3, g2, g1, g0])
	# three solutions are provided, the first two are complex, only the real solution is valid
	return roots[2].real
	 
	 */


	double NRpolarSolver::solve3rdDegreePolynomial(
			double d,
			double c,
			double b,
			double a,
			double x)
			const
	{
		double fx = a * x * x * x + b * x * x + c * x + d;
		double fxd = 3.0 * a * x * x + 2.0 * b * x + c;
		double incx = fx / fxd;

		while (abs(incx) > tolerance) {
			x -= incx;
			fx = a * x * x * x + b * x * x + c * x + d;
			fxd = 3.0 * a * x * x + 2.0 * b * x + c;
			incx = fx / fxd;
		}

		return x;
	}
	
	double NRpolarSolver::mu(mat &J, mat &J2, vec &F, vec &dV, vec &dD, vec & dx, uint npq, uint npv){
		Jacobian(J2, dV, dD, npq, npv);
		
		vec a = F;
		vec b = J * (dx);
		
		vec c(2*npq+npv); //= dx. * b * 0.5;
		for (uint i=0;i<(2*npq+npv); i++) //this loop is because EIGEN does not want to perform this simple element wise vector multiplication...
			c(i) = dx.coeff(i) * b.coeff(i) * 0.5;

		double g0 = -1* a.dot(b);
		double g1 = b.dot(b) + 2 * a.dot(c);
		double g2 = -3.0 * b.dot(c);
		double g3 = 2.0 * c.dot(c);
		   
		double sol = solve3rdDegreePolynomial(g3, g2, g1, g0, 1.0);
		return sol;         
	}

	void NRpolarSolver::get_power_inc(vec& PQinc, uint npq, uint npv) {

		uint npqpv = npq + npv;
		uint k;
		PQinc.setZero();

		for (uint a = 0; a < npqpv; a++) {
			//For PQ and PV buses; calculate incP
			k = PQPV[a];
			PQinc(a) = Pesp.coeff(k) - P(k);

			if (a < npq) //only for PQ buses, calculate incQ
				PQinc(a + npqpv) = Qesp.coeff(k) - Q(k);
		}
	}

	bool NRpolarSolver::converged(const vec& PQinc, uint npqpvpq) const
	{
		for (uint k = 0; k < npqpvpq; k++)
			if (abs(PQinc.coeff(k)) > tolerance)
				return false;

		return true;
	}
	
	void NRpolarSolver::get_increments(vec X, vec &incV, vec &incD, uint npq, uint npv){
	
		uint npqpv = npq + npv;
		uint k;

		for (uint a = 0; a < npqpv; a++) {
			k = PQPV[a];
			incD(k) = X.coeff(a);

			if (a < npq) {
				incV(k) = X.coeff(a + npqpv);
			}
		}
	}

	void NRpolarSolver::update_solution(vec X, uint npq, uint npv) {

		uint npqpv = npq + npv;
		uint k;

		for (uint a = 0; a < npqpv; a++) {
			k = PQPV[a];
			Sol.Varg(k) += X.coeff(a);

			if (a < npq) {
				Sol.Vmag(k) = Sol.Vmag.coeff(k) * (1.0 + X.coeff(a + npqpv));
			}
		}

		// Correct for PV buses
		for (uint i = npq; i < npq + npv; i++) {
			k = PQPV[i];
			cx_double v = Sol.V(k);
			v *= Model.buses[k].v_set_point / abs(v);
			Sol.Vmag(k) = abs(v);
			Sol.Varg(k) = arg(v);
		}
	}


	Solver::Result NRpolarSolver::solve(bool printIterations) {
		uint npq = Model.loadBusIndices.size();
		uint npv = Model.generatorBusIndices.size();
		uint npqpvpq = 2 * npq + npv;

		//System : J*X = K
		mat J(npqpvpq, npqpvpq);
		mat J2(npqpvpq, npqpvpq);
		vec K(npqpvpq);
		vec X(npqpvpq);
		vec incV(Sol.length());
		vec incD(Sol.length());
		
		// First shot: Perhaps the model already converged?

		get_power_inc(K, npq, npv);
		bool didConverge = converged(K, npqpvpq);

		for (unsigned int i = 0; (i < maxIterations) && !didConverge; ++i) {
			Jacobian(J, Sol.Vmag, Sol.Varg, npq, npv);

			Eigen::FullPivLU<mat> lu(J); //Full pivot LU
			X = lu.solve(K);
			
			get_increments(X, incV, incD, npq, npv);
			
			auto mu_ = mu(J, J2, K, incV, incD, X, npq, npv);
			
			//upgrade the solution
			update_solution(X * mu_, npq, npv);

			//Calculate the increment of power for the new iteration
			get_power_inc(K, npq, npv);

			didConverge = converged(K, npqpvpq);

			if (printIterations) {
				std::cout << "Iteration " << i << std::endl;
			}
		}
		
		//Calculate the reactive power for the PV buses:
		calculate_Q(npq, npv);

		if (didConverge) {
			calculate_slack_power();
			Model.set_solution(Sol.toComplex());
			return Solver::Solved;
		}
		else {
			return Solver::NotSolved;
		}
	}
	
	/**
	 * This function updates the solver solution object power values using the
	 * circuit's own solution power values. this is specially usefull when updating
	 * the circuit power values while keeping the previous voltage solution
	 */
	void NRpolarSolver::update_solution_power_from_circuit(){    
		Sol.P = Model.get_initial_solution().P;
		Sol.Q = Model.get_initial_solution().Q;
		Pesp = Sol.P;
		Qesp = Sol.Q;
	}
}
