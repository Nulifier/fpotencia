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

    /*
     * constructor
     */
    NRpolarSolver::NRpolarSolver(Circuit const& model):
            Model(model),
            tolerance(DEFAULT_SOLUTION_TOLERANCE),
            maxIterations(DEFAULT_MAX_ITERATIONS),
            Sol(Model.get_initial_solution())
    {
        if (!Sol.initialized) {
            Model.compile(false);
            Sol = Model.get_initial_solution();
        }


        BUSES.reserve(
                Model.loadBusIndices.size()
                    + Model.generatorBusIndices.size());
        BUSES.insert(
                BUSES.end(),
                Model.loadBusIndices.begin(),
                Model.loadBusIndices.end());
        BUSES.insert(
                BUSES.end(),
                Model.generatorBusIndices.begin(),
                Model.generatorBusIndices.end());

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

        if (!checks()) {
            throw std::invalid_argument(
                    "The circuit failed the solver compatibility test.");
        }
    }


    NRpolarSolver::NRpolarSolver(
            const Circuit& model,
            const solution& sol_):
                NRpolarSolver(model)
    {
        Sol = sol_;
    }


    NRpolarSolver::~NRpolarSolver() noexcept
    {
    }


    bool NRpolarSolver::checks() const
    {
        return Model.slackBusIndices.size() <= 1;
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the slack bus power
     */
    void NRpolarSolver::calculate_slack_power() {
        for (auto k: Model.slackBusIndices) {
            cx_double I(0.0, 0.0);
            for (uint j = 0; j < Model.buses.size(); j++) {
                I += Model.Y.coeff(k, j) * Sol.Vcx(j);
            }
            I = Sol.Vcx(k) * conj(I); //now this is the power
            Sol.P(k) = I.real();
            Sol.Q(k) = I.imag();
        }
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the reactive power of the bus k (usefull for PV uses)
     */
    void NRpolarSolver::calculate_Q(uint npq, uint npv) {
        double val;
        uint k;
        for (uint i = npq - 1; i < npq + npv; i++) {
            k = PQPV[i];
            val = Q(k);
            Sol.Q(k) = val;
        }
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the active power at a bus
     */
    double NRpolarSolver::P(uint k) {
        double val = 0.0;
        for (uint j = 0; j < Model.buses.size(); j++) {
            val += Sol.V.coeff(j)
                    *(Model.G(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j))
                    + Model.B(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j)));
        }
        return Sol.V.coeff(k) * val;
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the reactive power at a bus
     */
    double NRpolarSolver::Q(uint k) {
        double val = 0.0;
        for (uint j = 0; j < Model.buses.size(); j++) {
            val += Sol.V.coeff(j)
                    *(Model.G(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j))
                    - Model.B(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j)));
        }
        return Sol.V.coeff(k) * val;
    }


    void NRpolarSolver::jacobian(mat &J, vec &V, vec &D, uint npq, uint npv)
    {
        uint npqpv = npq + npv;
        double val;
        uint k, j;
        uint da, db;

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
                    //if (val != 0.0)
                    J(a, b) = val;
                }
            }
        }

        //J2
        da = 0;
        db = npqpv;
        for (uint a = 0; a < npqpv; a++) { //rows
            k = PQPV[a];
            //diagonal
            //std::cout << "J2D:" << (a + da) << "," << (a + db) << std::endl;
            if (a < npq)
                J(a + da, a + db) = P(k) + Model.G(k, k) * V.coeff(k) * V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npq; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = V.coeff(k) * V.coeff(j)
                            *(Model.G(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j))
                            + Model.B(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j)));
                    //if (val != 0.0)
                    //std::cout << "J2ij:" << (a + da) << "," << (b + db) << std::endl;
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
            //std::cout << "J3:" << (a + da) << "," << (a + db) << std::endl;
            J(a + da, a + db) = P(k) - Model.G(k, k) * V.coeff(k) * V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npqpv; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = V.coeff(k) * V.coeff(j)
                            *(Model.G(k, j) * cos(D.coeff(k) - D.coeff(j))
                            + Model.B(k, j) * sin(D.coeff(k) - D.coeff(j)));
                    //if (val != 0.0)
                    //std::cout << "J3:" << (a + da) << "," << (b + db) << std::endl;
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
            //std::cout << "J4:" << (a + da) << "," << (a + db) << std::endl;
            J(a + da, a + db) = Q(k) - Model.B(k, k) * V.coeff(k) * V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npq; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = V.coeff(k) * V.coeff(j)
                            *(Model.G(k, j) * sin(D.coeff(k) - D.coeff(j))
                            - Model.B(k, j) * cos(D.coeff(k) - D.coeff(j)));
                    if (val != 0.0) {
                        //std::cout << "J4:" << (a + da) << "," << (b + db) << std::endl;
                        J(a + da, b + db) = val;
                    }
                }
            }
        }


    }


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


    double NRpolarSolver::mu(
            mat const& jacobian,
            vec const& mismatches,
            vec& dx,
            size_t numLoads,
            size_t numGenerators)
    {
        vec a = mismatches;
        vec b = jacobian * (dx);
        vec c(dx.size());
        assert(static_cast<size_t>(c.size()) == 2*numLoads + numGenerators);
        assert(dx.size() == b.size());

        for (vec::Index i = 0; i < dx.size(); i++) {
            c(i) = dx.coeff(i) * b.coeff(i) * 0.5;
        }

        double g0 = -1.0 * a.dot(b);
        double g1 = b.dot(b) + 2.0 * a.dot(c);
        double g2 = -3.0 * b.dot(c);
        double g3 = 2.0 * c.dot(c);

        return solve3rdDegreePolynomial(g3, g2, g1, g0, 1.0);
    }


    void NRpolarSolver::adjustDeltaPQ(
            vec& pqDeltas,
            uint numLoads,
            uint numGenerators)
    {
        assert(pqDeltas.size() == 2 * numLoads + numGenerators);
        auto npqpv = numLoads + numGenerators;
        pqDeltas.setZero();

        for (vec::Index i = 0; i < npqpv; i++) {
            // For both buses, calculate the delta P:

            auto k = PQPV[i];
            pqDeltas(i) = Pesp.coeff(k) - P(k);

            if (i < numLoads) { // Only for PQ buses, also calculate Q:
                pqDeltas(i + npqpv) = Qesp.coeff(k) - Q(k);
            }
        }
    }


    bool NRpolarSolver::converged(const vec& PQinc) const
    {
        return PQinc.lpNorm<Infinity>() < tolerance;
    }


    void NRpolarSolver::get_increments(vec X, vec &incV, vec &incD, uint npq, uint npv){

        uint npqpv = npq + npv;

        for (uint a = 0; a < npqpv; a++) {
            auto k = PQPV[a];
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
            Sol.D(k) += X.coeff(a);

            if (a < npq)
                Sol.V(k) = Sol.V.coeff(k) * (1.0 + X.coeff(a + npqpv));
        }

        //Correct for PV buses
        for (uint i = npq; i < npq + npv; i++) {
            k = PQPV[i];
            cx_double v = Sol.Vcx(k);
            v *= Model.buses[k].v_set_point / abs(v);
            Sol.V(k) = abs(v);
            Sol.D(k) = arg(v);
        }
    }


    Solver::Result NRpolarSolver::powerFlow(Circuit& grid)
    {
        uint npq = Model.loadBusIndices.size();
        uint npv = Model.generatorBusIndices.size();
        uint npqpvpq = 2 * npq + npv;

        //System : J*X = K
        mat J(npqpvpq, npqpvpq);
        mat J2(npqpvpq, npqpvpq);
        vec X(npqpvpq);
        vec K(npqpvpq);
        vec incV(Sol.Lenght);
        vec incD(Sol.Lenght);

        // First shot: Perhaps the model already converged?

        adjustDeltaPQ(K, npq, npv);
        auto didConverge = converged(K);

        for (unsigned i = 0; i < maxIterations && ! didConverge; ++i) {
            jacobian(J, Sol.V, Sol.D, npq, npv);

            Eigen::FullPivLU<mat>lu(J); //Full pivot LU
            X = lu.solve(K);

            get_increments(X, incV, incD, npq, npv);

            // Update solution, use optimization factor µ:

            update_solution(X * mu(J, K, X, npq, npv), npq, npv);

            //Calculate the delta P/delta Q values:

            adjustDeltaPQ(K, npq, npv);

            didConverge = converged(K);
        }

        //Calculate the reactive power for the PV buses:

        calculate_Q(npq, npv);

        if (! didConverge) {
            return Solver::NotSolved;
        } else {
            calculate_slack_power();
            grid.set_solution(Sol.get_cx());
            return Solver::Solved;
        }
    }
}
