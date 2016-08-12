/*
 * File:   Circuit.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:06
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <complex>
#include <algorithm>

#include "Bus.h"
#include "Load.h"
#include "Line.h"
#include "Transformer.h"
#include "Generator.h"
#include "Shunt.h"

#include "Circuit.h"


using namespace Eigen;


namespace fPotencia {

    /*
     * Circuit class default onstructor
     */
    Circuit::Circuit(): default_voltage(1.0, 0.0)
    {
    }


    Circuit::~Circuit() noexcept
    {
    }

    /*
     * Add bus to the class.
     *
     * This function adds a bus object and assigns an ordinal number to the bus
     * object incomming to the function, so that object, even externally will
     * be refferenced to this circuit.
     *
     * I want to force all the buses to have different names.
     *
     * This is usefull to tell the branch elements the indices of the buses
     * where it is connected to.
     */
    void Circuit::add_Bus(Bus &bus)
    {
        if (buses.empty()) {
            bus.index = 0;
        } else {
            bus.index = buses.back().index + 1;
        }

        buses.emplace_back(bus);
    }


    /*
     * This function composes the circuit admittance matrix
     *
     * Each circuit branch element has a function to compose its own
     * admittance matrix. As each branch element "knows" the indices of the
     * busbars where it is connected, it can create an admittance matrix of the
     * dimension of the crcuit admittance matrix. If those elemenr Y matrices
     * are Sparse, then the mis use of memory is minimal ans the composition
     * of the circuit Y matriz becomes super simple: the sum of all the
     * branc elements Y_element matrices created as sparse matrices.
     */
    void Circuit::compose_Y() {
        auto n = Circuit::buses.size();

        if (n > 0) {
            Y = sp_cx_mat(n, n);
            double vbase;
            //Ymod = sp_mat(n, n);
            //Yang = sp_mat(n, n);
            //Add the shunt admittances
            for (uint i = 0; i < shunts.size(); i++) {
                shunts[i].get_element_Y(n, Y);
            }

            //Add the Lines admittances
            for (uint i = 0; i < lines.size(); i++) {
                //if the line values are not in per units, then make them per unit
                if (lines[i].valueType() == fPotencia::absolute) {
                    //do the per unit base change
                    vbase = buses[lines[i].buses().first].nominal_voltage;
                    lines[i].impedanceBase_ = (vbase * vbase) / Sbase;
                }

                assert(n > lines[i].buses().first);
                assert(n > lines[i].buses().second);
                lines[i].placeAdmittance(Y);
            }

            //Add the transforers admittance matrices.
            //The transformer model is formulated in such way that the values
            //come already in per unit
            for (uint i = 0; i < transformers.size(); i++) {
                transformers[i].get_element_Y(n, Y);
                //std::cout << "\n\n" << std::endl;
            }

        } else {
            throw std::invalid_argument("There are no buses");
        }
    }

    /*
     * Removes a column of zeros at the index k
     * Removes a row of zeros to the index k
     */
    void Circuit::removeCross(cx_mat& matrix, uint k) {

        //assumed that matrix is square as Y and Z are.
        uint n = matrix.cols() - 1;
        cx_mat m(n, n);
        uint a, b, i, j;
        for (i = 0; i < n; i++) {
            if (i < k)
                a = i;
            else
                a = i + 1;

            for (j = 0; j < n; j++) {

                if (j < k)
                    b = j;
                else
                    b = j + 1;

                m(i, j) = matrix(a, b);
            }
        }

        matrix = m;
    }

    /*
     * Adds a column of zeros at the index k
     * Adds a row of zeros to the index k
     */
    void Circuit::expandOnPoint(cx_mat& matrix, uint k) {

        //assumed that matrix is square as Y and Z are.
        uint n = matrix.cols();
        cx_mat m(n + 1, n + 1);
        uint a, b, i, j;
        for (i = 0; i < n; i++) {
            if (i < k)
                a = i;
            else
                a = i + 1;

            for (j = 0; j < n; j++) {

                if (j < k)
                    b = j;
                else
                    b = j + 1;

                m(a, b) = matrix(i, j);
            }
        }

        cx_double zero(0, 0);
        for (i = 0; i < n + 1; i++) {
            m(k, i) = zero;
            m(i, k) = zero;
        }

        matrix = m;
    }

    /*
     * This function creates the reduced impedance matrix of the circuit.
     *
     * This matrix contains a null column and row in the indices of the slack
     * (or VD) buses. This is usefull for the numerical algorithms to reach
     * a good solution
     *
     * Only applicable after Y is created
     */
    void Circuit::compose_Z() {
        //This is the reduced version of Z.
        //The Z-bus methods require that the slack bus column and row are removed

        //Create a copy of the full admittance matrix
        cx_mat Yred(Y);

        //remove a column and a row in those indices matching the slack buses indices
        for (uint i = 0; i < slackBusIndices.size(); i++)
            removeCross(Yred, slackBusIndices[i]);

        // perform a fast LU inverse of the reduced admittance matrix
        //this leads to the reduced impedance matrix (of the wrong dimensions)
        Eigen::FullPivLU<cx_mat> lu(Yred);
        Zred = lu.inverse();

        //to make it of the correct dimensions, lets add zero rows and columns
        //where there were removed at the beginning
        for (uint i = 0; i < slackBusIndices.size(); i++)
            expandOnPoint(Zred, slackBusIndices[i]);


        Eigen::FullPivLU<cx_mat> lu2(Y);
        Z = lu2.inverse();
    }

    /* Given the circuir objects, it build the relations among them: Impedance
     * and admittance matrices.
     *
     * This should be called every time the circuit topology changes
     */
    void Circuit::compile(bool guess_angles)
    {
        // Make sure we do not double-add:

        loadBusIndices.clear();
        slackBusIndices.clear();
        generatorBusIndices.clear();

        double maxpower = realPowerMaximum();
        Sbase = pow(10, 1 + floor(log10(maxpower)));

        // Set bus types:


        // All external grids make their respecitve buses VD buses:

        for (const auto& externalGrid: externalGrids) {
            auto& bus = buses.at(externalGrid.bus);
            if (bus.Type == undefined_bus_type) {
                bus.Type = VD;
            }
        }

        /* Iterate over all generator buses and make them PV buses if the
         * generator is indeed supplying voltage; otherwise, the bus is a PQ
         * bus.
         */

        for (auto const& generator: generators) {
            auto& bus = buses.at(generator.busIndex());

            if (undefined_bus_type == bus.Type) {
                bus.Type = PV;
            }

            // Apply generator limits conversion to the buses:

            bus.min_q = generator.reactivePowerLimits().first / Sbase;
            bus.max_q = generator.reactivePowerLimits().second / Sbase;

            if (Generator::pu == generator.voltageType()) {
                bus.v_set_point = generator.voltage();
            } else { // absolute volts:
                bus.v_set_point = generator.voltage() / bus.nominal_voltage;
            }
        }

        // All buses which carry loads and are not yet classified, are PQs:

        for (auto const& load: loads) {
            auto& bus = buses.at(load.busIndex());
            if (bus.Type == undefined_bus_type) {
                bus.Type = PQ;
            }
        }

        // ... finally, all other, non-qualified buses are PQ buses:

        Vbase = 0.0;
        for (auto& bus: buses) {
            if (bus.Type == undefined_bus_type) {
                bus.Type = PQ;
            }


            // Rise Vbase, if necessary:

            if (bus.nominal_voltage > Vbase) {
                Vbase = bus.nominal_voltage;
            }


            // Build the bus lists:

            switch (bus.Type) {
            case VD:
                slackBusIndices.push_back(bus.index);
                break;
            case PQ:
                loadBusIndices.push_back(bus.index);
                break;
            case PV:
                generatorBusIndices.push_back(bus.index);
                break;
            default:
                throw std::invalid_argument("Unknown bus type");
            }
        }

        Zbase = Vbase * Vbase / Sbase;

        Ybase = 1.0 / Zbase;

        //Calculate the power connected to each bus according to the type
        generate_initial_solution();

        //create the admittance matrix
        compose_Y();

        //create the admittance matrix
        compose_Z();

        if (guess_angles)
            //Calculate the power connected to each bus according to the type (this time with the new voltage)
            correct_initial_solution();
    }


    double Circuit::realPowerMaximum() const
    {
        double maximum = 0.0;

        for (const auto& generator: generators) {
            if (std::abs(generator.power()) > maximum) {
                maximum = std::abs(generator.realPower());
            }
        }

        for (const auto& load: loads) {
            if (std::abs(load.power()) > maximum) {
                maximum = std::abs(load.realPower());
            }
        }

        return maximum;
    }


    /*
     * Generates the Circuit initial solution.
     *
     * The circuit object contains two solutions that are the same;
     * - sol: the polar solution, where the voltage is in polar mode and the
     *        power is in cartesian mode
     * - cx_sol: where the voltage and the power are in cartesian mode
     *
     * Both solutions contain the circuit initial guess of the solution
     * in per unit values of the voltage and the power. This translates to
     * - Voltage = 1 +0j
     * - Power = sum of the bus connected load and generation power divided by
     *           the base power Sbase.
     *
     * Sbase is calculated to match the order of the circuit power values.
     */
    void Circuit::generate_initial_solution(bool keep_last_solution)
    {
        // Reset all powers:

        for (auto& bus: buses) {
            bus.connected_power = cx_double(0.0, 0.0);
        }

        // Calculate the bus connected generation and load:

        for (auto const& generator: generators) {
            auto& bus = buses.at(generator.busIndex());
            bus.connected_power += generator.power();
        }

        for (const auto& load: loads) {
            auto& bus = buses.at(load.busIndex());
            bus.connected_power -= load.power();
        }

        sol.resize(buses.size());
        cx_sol.resize(buses.size());

        for (uint i = 0; i < buses.size(); i++) {

            if (buses[i].Type == PQ) {

                sol.P(i) = buses[i].connected_power.real() / Sbase;
                sol.Q(i) = buses[i].connected_power.imag() / Sbase;
                if (!keep_last_solution) {
                    sol.V(i) = default_voltage.real();
                    sol.D(i) = default_voltage.imag();
                }

            } else if (buses[i].Type == PV) {

                sol.P(i) = buses[i].connected_power.real() / Sbase;
                sol.V(i) = default_voltage.real();
                if (!keep_last_solution) {
                    sol.Q(i) = buses[i].connected_power.imag() / Sbase;
                    sol.D(i) = default_voltage.imag();
                }

            } else if (buses[i].Type == VD) {
                //The slack values must always be initilaized to keep the
                //same voltage refference for the solvers
                sol.P(i) = 0.0;
                sol.Q(i) = 0.0;
                sol.V(i) = default_voltage.real();
                sol.D(i) = default_voltage.imag();
            }

            cx_sol.S(i) = cx_double(sol.P[i], sol.Q[i]);
            cx_sol.V(i) = cx_double(sol.V[i], sol.D[i]);
        }

        sol.initialized = true;
        cx_sol.initialized = true;
    }

    /*
     * This class generates an initial estimate of the voltage angles by
     * running a DC Power Flow simulation. This simulation is reduced
     * to the multiplication of the Z matrix of the circuit by the active
     * power vector due to a number of assumntions (quite unrealistic)
     * This voltage angles solution is quite usefull to initialize the
     * solution to be sent to bigger AC solvers.
     */
    void Circuit::correct_initial_solution() {

        //this generates the DC-Solution angles
        dc_angles = Zred.imag() * cx_sol.getP();
        cout << "Initial D:\n" << dc_angles << endl;

        for (uint i = 0; i < buses.size(); i++) {
            sol.D(i) = dc_angles(i);

            cx_sol.V(i) = cx_double(sol.V[i], sol.D[i]);
        }

        sol.initialized = true;
        cx_sol.initialized = true;
        sol.print("Initial Solution:");
        cx_sol.print("Initial CX Solution:");
    }

    /*
     * Return the circuit initial solution in polar mode
     */
    solution Circuit::get_initial_solution() {
        return sol;
    }

    /*
     * Return the circuit initial solution in complex mode
     */
    cx_solution Circuit::get_initial_cx_solution() {
        return cx_sol;
    }

    /*
     * Calls to calculate the branch elements current and power flow.
     * Each branch element has a function that obtains the current and power
     * flows given the terminals volatge and the element admittance matrix.
     * The admittance matrix is known from the previous step of generating
     * the circuit admittance matrix
     */
    void Circuit::calculate_flows(cx_solution sol_) {

        for (uint i = 0; i < lines.size(); i++) {
            lines[i].setPower(sol_);
        }

        for (uint i = 0; i < transformers.size(); i++) {
            transformers[i].calculate_current(sol_);
        }

        for (uint i = 0; i < shunts.size(); i++) {
            shunts[i].calculate_current(sol_);
        }
    }

    /*
     * Sets the circuit solution and applies the solution to the
     * circuit objects. For example, it sets the busbars voltages
     * and calculates the lines power flows, which are accesible at
     * the line objects.
     */
    void Circuit::set_solution(cx_solution sol_) {

        //Set the main circuit solution
        cx_sol = sol_;

        //Undo the power base change
        cx_double s_base(Sbase, 0);
        for (uint i = 0; i < sol_.Lenght; i++)
            sol_.S(i) *= s_base;

        //Undo the voltage change and assign the buses power and voltage
        for (uint i = 0; i < sol_.Lenght; i++) {
            buses[i].voltage_pu = sol_.V.coeff(i);
            sol_.V(i) *= cx_double(buses[i].nominal_voltage, 0);
            buses[i].voltage = sol_.V.coeff(i);
            buses[i].power = sol_.S.coeff(i);

            //copy cx_sol to sol
            sol.P(i) = cx_sol.S.coeff(i).real();
            sol.Q(i) = cx_sol.S.coeff(i).imag();
            sol.V(i) = abs(cx_sol.V(i));
            sol.D(i) = arg(cx_sol.V(i));
        }

        calculate_flows(sol_);
    }

    /*
     * Gets the real part of a circuit admittance matrix element
     * at row i and column j
     */
    double Circuit::G(int i, int j) {
        //cx_double com = (cx_double) (Y.coeff(i, j));
        return Y.coeff(i, j).real();
    }

    /*
     * Gets the imaginary part of a circuit admittance matrix element
     * at row i and column j
     */
    double Circuit::B(int i, int j) {
        //cx_double com = (cx_double) (Y.coeff(i, j));
        return Y.coeff(i, j).imag();
    }

    /*
     * Prints the circuit element power flows and voltages
     */
    void Circuit::print() {

        cout << "Sbase = " << Sbase << endl;

        for (uint i = 0; i < transformers.size(); i++)
            transformers[i].print();

        for (uint i = 0; i < shunts.size(); i++)
            shunts[i].print();

        for (uint i = 0; i < buses.size(); i++)
            buses[i].print();
    }

    /**/
    void Circuit::print_buses_state() {

        string type;
        for (uint i = 0; i < buses.size(); i++) {
            type = BusType_name[buses[i].Type];
            cout << i << ": " + type << endl;
        }
    }

    /**/
    void Circuit::printCXMat(cx_mat m, string header) {
        std::cout << header << std::endl;
        for (uint i = 0; i < buses.size(); i++)
            for (uint j = 0; j < buses.size(); j++)
                if (m.coeff(i, j) != cx_double(0, 0))
                    std::cout << "(" << i << "," << j << ") = " << m.coeff(i, j) << std::endl;

    }

    void Circuit::printMat(sp_mat m, string header) {
        std::cout << header << std::endl;
        for (uint i = 0; i < buses.size(); i++)
            for (uint j = 0; j < buses.size(); j++)
                if (m.coeff(i, j) != 0.0)
                    std::cout << "(" << i << "," << j << ") = " << m.coeff(i, j) << std::endl;

    }

    /*
     * Checks how much does the solution diverge
     *
     * The vector S has to be equal to Vx(YxV)* for a perfect solution.
     * The vector S is calculated by any of the AC or DC solvers
     *
     * S = VxI* is the most fundametal equation in the power flow simulation
     * if it is fullfilled the circuit has reached a valuable solution,
     * however, because of the use of numerical algorithms, the equality is
     * never reached and therefoe we need to stand a certain threshold of
     * inequality: S = threshold * (VxI*) where I* = (YxV)*
     */
    void Circuit::check_solution() {

        uint n = buses.size();
        cx_vec V = cx_sol.V;
        cx_vec S = cx_sol.S;
        cx_mat YVconj(n, 1);
        sp_cx_mat Vdiag(n, n);

        cx_mat A = Y*V;

        double r, im;
        for (uint i = 0; i < n; i++) {
            r = A.coeff(i, 0).real();
            im = A.coeff(i, 0).imag();
            YVconj(i, 0) = cx_double(r, -im);
            Vdiag.insert(i, i) = V(i, 0);
        }

        cx_mat delta = S - Vdiag*YVconj; //if delta is zero for all the values the solution is perfect

        double mismatch = 0;
        for (uint i = 0; i < n; i++)
            mismatch += abs(delta(i, 0));

        mismatch /= n;
    }
}
