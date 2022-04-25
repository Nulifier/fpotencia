/*!
 * \file Solver_NRpolar.h
 * \author Santiago Peñate Vera
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

#include "Solver.h"
#include "Circuit.h"
#include "Solution.h"

namespace fPotencia {
	/**
	 * This class implements the Newton-Raphson method of load flow analysis using polar coordinates.
	 */
	class NRpolarSolver: public Solver {
	public:
		/**
		 * Creates a new Solver.
		 *
		 * \param[in] model The circuit the solver should perform load flow
		 *  analysis on
		 */
		NRpolarSolver(Circuit& model);

		/**
		 * Constructs a new solver instance that works off an initial solution.
		 *
		 * \param[in] model The circuit the solver should perform load flow
		 *  analysis on
		 *
		 * \param[in] sol_ The initial solution the solver should start
		 *  working with
		 */
		NRpolarSolver(Circuit& model, const solution& sol_);

		/// The circuit model the solver analyses.
		Circuit& Model;

		/**
		 * Allowable tolerance of the solver instance.
		 */
		double tolerance;

		/**
		 * Maximum number of iterations.
		 */
		unsigned maxIterations;

		/**
		 * Solves a polynomial of 3rd degree.
		 *
		 * This method solves a polynomial defined by the coeffients
		 * g0, g1, g3 and g3 such that $d + c*x + b*x^2 + a*x^3 = 0$.
		 *
		 * Provides the real solution using the Newton-Raphson technique
		 */
		double solve3rdDegreePolynomial(
				double d,
				double c,
				double b,
				double a,
				double x) const;

		/// Checks whether a particular solution converged.
		bool converged(vec const& PQinc, uint npqpvpq) const;

		/*!
		 * \brief Solves the grid
		 *
		 * \return Solver_State::converged if the grid was solved
		 *
		 * \sa Solver_State
		 */
		Solver::Result solve(bool printIterations = false) override;

		[[nodiscard]] bool canSolve() const override;

		void update_solution_power_from_circuit();
		
	private:
		/// Vector of indices of the pq and pv busses.
		std::vector<int> PQPV;

		/// List of PQ bus indices in the previous iteration.
		std::vector<int> LastPQ;
		
		/// List of PV bus indices in the previous iteration.
		std::vector<int> LastPV;

		/// Vector of specified active power.
		vec Pesp;

		/// Vector of specified reactive power.
		vec Qesp;

		solution Sol;

		/// Calculate the Jacobian of the circuit.
		void Jacobian(mat &J, vec &V, vec &D, uint npq, uint npv);
		
		double mu(const mat &J, mat &J2, const vec &F, vec &dV, vec &dD, vec & dx, uint npq, uint npv);

		/// Calculate the power increments.
		void get_power_inc(vec &PQinc, uint npq, uint npv); //PQinc is passed by refference

		/// Calculate the reactive power of the bus k (usefull for PV uses).
		void calculate_Q(uint npq, uint npv); //calculate the reative power at the PV buses

		/// Calculate the active power at a bus.
		double P(uint k);

		/// Calculate the reactive power at a bus.
		double Q(uint k);

		void update_solution(vec X, uint npq, uint npv);
		
		void get_increments(vec X, vec &incV, vec &incD, uint npq, uint npv);

		/// Calculate the slack bus power
		void calculate_slack_power();

		/// This function corects the PV buses that exeed the reative power limit.
		void correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X);
	};
}
