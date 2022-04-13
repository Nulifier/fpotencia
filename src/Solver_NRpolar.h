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

	/*!
	 * \brief This class implements the Newton-Raphson method of load flow
	 *  analysis using polar coordinates.
	 */
	class NRpolarSolver: public Solver
	{
	public:


		/*!
		 * \brief Creates a new solver
		 *
		 * \param[in] model The circuit the solver should perform load flow
		 *  analysis on
		 */
		NRpolarSolver(Circuit const& model);


		/*!
		 * \brief Constructs a new solver instance that works off an initial
		 *  solution
		 *
		 * \param[in] model The circuit the solver should perform load flow
		 *  analysis on
		 *
		 * \param[in] sol_ The initial solution the solver should start
		 *  working with
		 */
		NRpolarSolver(Circuit const& model, solution const& sol_);


		virtual ~NRpolarSolver() noexcept;


		//!  \brief The circuit model the solver analyses
		Circuit Model;


		/*!
		 * \brief Allowable tolerance of the solver instance
		 *
		 * \sa DEFAULT_SOLUTION_TOLERANCE
		 */
		double tolerance;


		/*!
		 * \brief Maximum number of iterations
		 *
		 * \sa DEFAULT_MAX_ITERATIONS
		 */
		unsigned maxIterations;


		/*!
		 * \brief Solves a polynomial of 3rd degree
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


		//! \brief Checks whether a particular solution converged
		bool converged(vec const& PQinc, uint npqpvpq) const;


		/*!
		 * \brief Solves the grid
		 *
		 * \return Solver_State::converged if the grid was solved
		 *
		 * \sa Solver_State
		 */
		virtual Solver::Result powerFlow(Circuit& grid) override;


		void update_solution_power_from_circuit();
		
	private:

		std::vector<int> BUSES;

		std::vector<int> PQPV;

		std::vector<int> LastPQ;
		
		std::vector<int> LastPV;

		vec Pesp;

		vec Qesp;

		solution Sol;

		void Jacobian(mat &J, vec &V, vec &D, uint npq, uint npv); //calculate the jacobian, J is passed by refference
		
		double mu(mat &J, mat &J2, vec &F, vec &dV, vec &dD, vec & dx, uint npq, uint npv);
		void get_power_inc(vec &PQinc, uint npq, uint npv); //PQinc is passed by refference

		void calculate_Q(uint npq, uint npv); //calculate the reative power at the PV buses

		double Q(uint k);

		double P(uint k);

		void update_solution(vec X, uint npq, uint npv);
		
		void get_increments(vec X, vec &incV, vec &incD, uint npq, uint npv);

		void calculate_slack_power(); //calculate the slack bus power        


		//! \brief Checks whether the solver can work on the given model
		bool checks() const;


		void correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X);
		
	};
}
