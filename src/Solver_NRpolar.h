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

#ifndef SOLVER_NRPOLAR_H
#define	SOLVER_NRPOLAR_H


#include <cmath>

#include "Solver.h"
#include "Circuit.h"
#include "Solution.h"


using namespace std;


namespace fPotencia {

    /*!
     * \brief This class implements the Nerwton-Raphson method of load flow
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
         * \brief Creates the Jacobian
         *
         * \param J
         * \param V
         * \param D
         * \param npq
         * \param npv
         */
        void jacobian(mat &J, vec &V, vec &D, uint npq, uint npv);


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
                double x)
                const;


        //! \brief Compute the acceleration factor
        double mu(
                mat const& jacobian,
                vec const& mismatches,
                vec& dx,
                size_t numLoads,
                size_t numGenerators);


        /*!
         * \brief Calculates the values for \f$ \Delta P \f$
         *  and \f$ \Delta Q \f$, i.e., the values by which each load and
         *  generator power flow value is adjusted.
         *
         * These values are the convergence criterium: If each of them is
         * below the tolerance set, the algorithm has converged and the solver
         * is finished.
         *
         * \param[in] numLoads Number of loads in the circuit
         *
         * \param[in] numGenerators Number of generators in the circuit
         *
         * \return The vector containing all delta P and delta Q values
         */
        void adjustDeltaPQ(vec& pqDeltas, uint numLoads, uint numGenerators);


        /*!
         * \brief Checks whether a particular solution converged
         *
         * \param[in] PQinc The vector containing the increments (deltas) of
         *  real and reactive power values for the next step of the iteration
         */
        bool converged(vec const& PQinc) const;


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

        vector<int> BUSES;

        vector<int> PQPV;

        vector<int> LastPQ;

        vector<int> LastPV;

        vec Pesp;

        vec Qesp;

        solution Sol;



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

#endif	/* SOLVER_NRPOLAR_H */

}
