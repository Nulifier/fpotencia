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

#include "Circuit.h"
#include <complex>

#include "Bus.h"
#include "Generator.h"
#include "Line.h"
#include "Load.h"
#include "Shunt.h"
#include "Transformer.h"
#include "enumaratons.h"

namespace fPotencia {
	void Circuit::add_Bus(Bus& bus) {
		if (buses.empty()) {
			bus.index = 0;
		}
		else {
			bus.index = buses.back().index + 1; // Sequencial bus numbering
		}

		buses.push_back(bus);
	}

	void Circuit::remove_Bus(Bus bus) {
		for (Buses::const_iterator itr = buses.begin(); itr != buses.end();
		     ++itr) {
			if (itr->index == bus.index) {
				buses.erase(itr);
				return;
			}
		}
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
		int n = static_cast<int>(Circuit::buses.size());

		if (n > 0) {
			Y = sp_cx_mat(n, n);
			// Add the shunt admittances
			for (auto& shunt : shunts) {
				shunt.get_element_Y(n, Y);
			}

			// Add the Lines admittances
			for (auto& line : lines) {
				// if the line values are not in per units, then
				// make them per unit
				if (!line.values_in_per_unit) {
					// do the per unit base change
					double vbase = buses[line.bus1].nominal_voltage;
					line.Zbase = (vbase * vbase) / Sbase;
				}

				line.get_element_Y(n, Y);
			}

			// Add the transforers admittance matrices.
			// The transformer model is formulated in such way that the
			// values come already in per unit
			for (auto& transformer : transformers) {
				transformer.get_element_Y(n, Y);
			}
		}
		else {
			throw std::invalid_argument("There are no buses");
		}
	}

	/*
	 * Removes a column of zeros at the index k
	 * Removes a row of zeros to the index k
	 */
	void Circuit::removeCross(cx_mat& matrix, uint k) {
		// assumed that matrix is square as Y and Z are.
		uint n = matrix.cols() - 1;
		cx_mat m(n, n);
		for (unsigned int i = 0; i < n; i++) {
			unsigned int a = i < k ? i : i + 1;

			for (unsigned int j = 0; j < n; j++) {
				unsigned int b = j < k ? j : j + 1;

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
		// assumed that matrix is square as Y and Z are.
		uint n = matrix.cols();
		cx_mat m(n + 1, n + 1);

		for (unsigned int i = 0; i < n; ++i) {
			unsigned int a = i < k ? i : i + 1;

			for (unsigned int j = 0; j < n; ++j) {
				unsigned int b = j < k ? j : j + 1;

				m(a, b) = matrix(i, j);
			}
		}

		cx_double zero(0, 0);
		for (unsigned int i = 0; i < n + 1; ++i) {
			m(k, i) = zero;
			m(i, k) = zero;
		}

		matrix = m;
	}

	void Circuit::compose_Zred() {
		// This is the reduced version of Z.
		// The Z-bus methods require that the slack bus column and row are
		// removed

		// Create a copy of the full admittance matrix
		cx_mat Yred(Y);

		// remove a column and a row in those indices matching the slack
		// buses indices
		for (unsigned int slackBusIndice : slackBusIndices) {
			removeCross(Yred, slackBusIndice);
		}

		// perform a fast LU inverse of the reduced admittance matrix
		// this leads to the reduced impedance matrix (of the wrong dimensions)
		Eigen::FullPivLU<cx_mat> lu(Yred);
		m_Zred = lu.inverse();

		// to make it of the correct dimensions, lets add zero rows and
		// columns where there were removed at the beginning
		for (unsigned int slackBusIndice : slackBusIndices) {
			expandOnPoint(*m_Zred, slackBusIndice);
		}
	}

	/*
	 * Only applicable after Y is created
	 */
	void Circuit::compose_Z() {
		Eigen::FullPivLU<cx_mat> lu2(Y);
		m_Z = lu2.inverse();
	}

	/* Given the circuir objects, it build the relations among them: Impedance
	 * and admittance matrices.
	 *
	 * This should be called every time the circuit topology changes
	 */
	void Circuit::compile(bool guess_angles) {
		// Make sure we do not double-add:

		loadBusIndices.clear();
		slackBusIndices.clear();
		generatorBusIndices.clear();

		// Calculate the base power in a comparable scale to the power in the
		// grid
		double maxpower = get_max_power();
		// cout << "max power:" << maxpower << endl;
		Sbase = pow(10, 1 + floor(log10(maxpower)));

		// Calculate the bus types

		// the presence of an external grid makes the bus VD (or slack)
		for (const auto& externalGrid : externalGrids) {
			if (buses[externalGrid.bus].Type == BusType::undefined_bus_type) {
				buses[externalGrid.bus].Type = BusType::VD;
			}
		}

		/*the presence of an generator makes the bus PV (if it is voltage
		 * controlled) or PQ otherwise
		 */
		for (const auto& generator : generators) {
			if (buses[generator.bus].Type == BusType::undefined_bus_type) {
				if (generator.voltage_controlled) {
					/* if the generator is set to contol the voltage
					 * then, make the bus a PV bus, otherwise leve
					 * the bus fo the last loop to make it a PQ bus
					 */
					buses[generator.bus].Type = BusType::PV;
				}
			}
			// apply generator limits conversion to the buses
			buses[generator.bus].min_q = generator.min_Q / Sbase;
			buses[generator.bus].max_q = generator.max_Q / Sbase;
			if (generator.Vset_in_per_unit) {
				buses[generator.bus].v_set_point = generator.voltage_set_point;
			}
			else {
				buses[generator.bus].v_set_point =
				    generator.voltage_set_point /
				    buses[generator.bus].nominal_voltage;
			}
		}

		// the presence of a load makes the bus PQ if it has not ben
		// classified
		for (auto& load : loads) {
			if (buses[load.bus].Type == BusType::undefined_bus_type) {
				buses[load.bus].Type = BusType::PQ;
				// PQ_list.push_back(buses[loads[i].bus].index);
			}
		}

		// not clasified buses are set to PQ mode
		Vbase = 0.0;
		for (auto& buse : buses) { // final check
			if (buse.Type == BusType::undefined_bus_type) {
				buse.Type = BusType::PQ;
			}

			if (buse.nominal_voltage > Vbase) {
				Vbase = buse.nominal_voltage; // set Vbase as the
				                              // biggest voltage
			}

			// set the bus types lists
			switch (buse.Type) {
			case BusType::VD:
				slackBusIndices.push_back(buse.index);
				break;
			case BusType::PQ:
				loadBusIndices.push_back(buse.index);
				break;
			case BusType::PV:
				generatorBusIndices.push_back(buse.index);
				break;
			default:
				throw std::invalid_argument("Unknown bus type");
			}
		}

		Zbase = Vbase * Vbase / Sbase;

		Ybase = 1.0 / Zbase;

		// Calculate the power connected to each bus according to the type
		generate_initial_solution();

		// create the admittance matrix
		compose_Y();

		// create the admittance matrix
		// compose_Z();

		if (guess_angles) {
			// Calculate the power connected to each bus according to the
			// type (this time with the new voltage)
			correct_initial_solution();
		}
	}

	double Circuit::get_max_power() const {
		double mx = 0;

		// Calculate the bus connected generation and load
		for (const auto& generator : generators) {
			if (abs(generator.power) > mx) {
				mx = abs(generator.power.real());
			}
		}

		for (const auto& load : loads) {
			if (abs(load.power) > mx) {
				mx = abs(load.power.real());
			}
		}

		return mx;
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
	void Circuit::generate_initial_solution(bool keep_last_solution) {

		// Initialize all powers
		for (auto& buse : buses) {
			buse.connected_power = cx_double(0.0, 0.0);
		}

		// Calculate the bus connected generation and load
		for (auto& generator : generators) {
			buses[generator.bus].connected_power += generator.power;
		}

		for (auto& load : loads) {
			buses[load.bus].connected_power -= load.power;
		}

		sol.resize(buses.size());
		cx_sol.resize(buses.size());

		for (uint i = 0; i < buses.size(); i++) {

			if (buses[i].Type == BusType::PQ) {

				sol.P(i) = buses[i].connected_power.real() / Sbase;
				sol.Q(i) = buses[i].connected_power.imag() / Sbase;
				if (!keep_last_solution) {
					sol.Vmag(i) = default_voltage.real();
					sol.Varg(i) = default_voltage.imag();
				}
			}
			else if (buses[i].Type == BusType::PV) {

				sol.P(i) = buses[i].connected_power.real() / Sbase;
				sol.Vmag(i) = default_voltage.real();
				if (!keep_last_solution) {
					sol.Q(i) = buses[i].connected_power.imag() / Sbase;
					sol.Varg(i) = default_voltage.imag();
				}
			}
			else if (buses[i].Type == BusType::VD) {
				// The slack values must always be initilaized to keep
				// the same voltage refference for the solvers
				sol.P(i) = 0.0;
				sol.Q(i) = 0.0;
				sol.Vmag(i) = default_voltage.real();
				sol.Varg(i) = default_voltage.imag();
			}

			cx_sol.S(i) = cx_double(sol.P[i], sol.Q[i]);
			cx_sol.V(i) = cx_double(sol.Vmag[i], sol.Varg[i]);
		}

		sol.initialized = true;
		cx_sol.initialized = true;
	}

	/*
	 * This class generates an initial estimate of the voltage angles by
	 * running a DC Power Flow simulation. This simulation is reduced
	 * to the multiplication of the Z matrix of the circuit by the active
	 * power vector due to a number of assumptions (quite unrealistic)
	 * This voltage angles solution is quite usefull to initialize the
	 * solution to be sent to bigger AC solvers.
	 */
	void Circuit::correct_initial_solution() {
		// this generates the DC-Solution angles
		dc_angles = Zred().imag() * cx_sol.P();

		for (uint i = 0; i < buses.size(); i++) {
			sol.Varg(i) = dc_angles(i);

			cx_sol.V(i) = cx_double(sol.Vmag[i], sol.Varg[i]);
		}

		sol.initialized = true;
		cx_sol.initialized = true;
	}

	/*
	 * Return the circuit initial solution in polar mode
	 */
	solution Circuit::get_initial_solution() { return sol; }

	/*
	 * Return the circuit initial solution in complex mode
	 */
	cx_solution Circuit::get_initial_cx_solution() { return cx_sol; }

	/*
	 * Calls to calculate the branch elements current and power flow.
	 * Each branch element has a function that obtains the current and power
	 * flows given the terminals volatge and the element admittance matrix.
	 * The admittance matrix is known from the previous step of generating
	 * the circuit admittance matrix
	 */
	void Circuit::calculate_flows(const cx_solution& sol_) {
		for (auto& line : lines) {
			line.calculate_current(sol_);
		}

		for (auto& transformer : transformers) {
			transformer.calculate_current(sol_);
		}

		for (auto& shunt : shunts) {
			shunt.calculate_current(sol_);
		}
	}

	/*
	 * Sets the circuit solution and applies the solution to the
	 * circuit objects. For example, it sets the busbars voltages
	 * and calculates the lines power flows, which are accesible at
	 * the line objects.
	 */
	void Circuit::set_solution(cx_solution sol_) {

		// Set the main circuit solution
		cx_sol = sol_;

		// Undo the power base change
		cx_double s_base(Sbase, 0);
		for (uint i = 0; i < sol_.length(); ++i) {
			sol_.S(i) *= s_base;
		}

		// Undo the voltage change and assign the buses power and voltage
		for (uint i = 0; i < sol_.length(); ++i) {
			buses[i].voltage_pu = sol_.V.coeff(i);
			sol_.V(i) *= cx_double(buses[i].nominal_voltage, 0);
			buses[i].voltage = sol_.V.coeff(i);
			buses[i].power = sol_.S.coeff(i);

			// copy cx_sol to sol
			sol.P(i) = cx_sol.S.coeff(i).real();
			sol.Q(i) = cx_sol.S.coeff(i).imag();
			sol.Vmag(i) = abs(cx_sol.V(i));
			sol.Varg(i) = arg(cx_sol.V(i));
		}

		calculate_flows(sol_);
	}

	void Circuit::print() const {
		std::cout << "Sbase = " << Sbase << std::endl;

		for (const auto& line : lines) {
			line.print();
		}

		for (const auto& transformer : transformers) {
			transformer.print();
		}

		for (const auto& shunt : shunts) {
			shunt.print();
		}

		for (const auto& buse : buses) {
			buse.print();
		}
	}

	void Circuit::print_buses_state() {
		for (uint i = 0; i < buses.size(); i++) {
			std::cout << i << ": "
			          << BusType_name.at(static_cast<size_t>(buses[i].Type))
			          << std::endl;
		}
	}

	void Circuit::printCXMat(const cx_mat& m, const std::string& header) const {
		std::cout << header << std::endl;
		for (uint i = 0; i < buses.size(); i++) {
			for (uint j = 0; j < buses.size(); j++) {
				if (m.coeff(i, j) != cx_double(0, 0)) {
					std::cout << "(" << i << "," << j << ") = " << m.coeff(i, j)
					          << std::endl;
				}
			}
		}
	}

	void Circuit::printMat(const sp_mat& m, const std::string& header) const {
		std::cout << header << std::endl;
		for (uint i = 0; i < buses.size(); i++) {
			for (uint j = 0; j < buses.size(); j++) {
				if (m.coeff(i, j) != 0.0) {
					std::cout << "(" << i << "," << j << ") = " << m.coeff(i, j)
					          << std::endl;
				}
			}
		}
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

		cx_mat A = Y * V;

		double r, im;
		for (uint i = 0; i < n; i++) {
			r = A.coeff(i, 0).real();
			im = A.coeff(i, 0).imag();
			YVconj(i, 0) = cx_double(r, -im);
			Vdiag.insert(i, i) = V(i, 0);
		}

		cx_mat delta = S - Vdiag * YVconj; // if delta is zero for all the
		                                   // values the solution is perfect

		double mismatch = 0;
		for (uint i = 0; i < n; i++)
			mismatch += abs(delta(i, 0));

		mismatch /= n;
	}

	/*
	 * This function sets the load and generation power values (in actual values
	 * not in p.u) and updates the solution objects
	 */
	/*void Circuit::setPowerValues(double loadP[], double loadQ[], double
	genP[], double genQ[]) {

	                                if (sizeof (loadP) != sizeof (loadQ))
	                                                                throw
	std::invalid_argument("setPowerValues: The size of the load vectors are
	different");

	                                if (sizeof (loadP) != loads.size())
	                                                                throw
	std::invalid_argument("setPowerValues: The size of the load vectors does not
	match the loads size");

	                                if (sizeof (genP) != sizeof (genQ))
	                                                                throw
	std::invalid_argument("setPowerValues: The size of the generator vectors are
	different");

	                                if (sizeof (genP) != generators.size())
	                                                                throw
	std::invalid_argument("setPowerValues: The size of the load vectors does not
	match the loads size");

	                                //Set the load values
	                                for (uint i = 0; i < loads.size(); i++)
	                                                                loads[i].power
	= cx_double(loadP[i], loadQ[i]);

	                                //Set the generation values
	                                for (uint i = 0; i < generators.size(); i++)
	                                                                generators[i].power
	= cx_double(genP[i], genQ[i]);

	                                //Update the solution but keeping the
	previous values of voltage (and Q, D in case of PV buses)
	generate_initial_solution(true);
	}*/
}
