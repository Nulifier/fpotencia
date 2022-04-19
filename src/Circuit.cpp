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
#include "Circuit.h"

#include "Bus.h"
#include "Load.h"
#include "Line.h"
#include "Transformer.h"
#include "Generator.h"
#include "Shunt.h"

using namespace Eigen;

namespace fPotencia {
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
	void Circuit::add_Bus(Bus &bus) {
		if (find_bus(bus.Name) == -1) {

			if (buses.empty())
				bus.index = 0;
			else
				bus.index = buses.back().index + 1; //Sequencial bus numbering

			buses.push_back(bus);
		} else {
			std::cout << "The bus " << bus.Name << " has been already added, therefore this object is not being added" << std::endl;
		}
	}

	/*
	 * Remove bus by object: Removes a bus from the list, but all the buses
	 * indices remain the same
	 */
	void Circuit::remove_Bus(Bus bus) {
		int index = find_bus(bus.Name);
		buses.erase(buses.begin() + index);
	}

	/*
	 * Remove bus by name: Removes a bus from the list, but all the buses
	 * indices remain the same
	 */
	void Circuit::remove_Bus(std::string bus_name) {
		int index = find_bus(bus_name);
		buses.erase(buses.begin() + index);
	}

	/*Get bus index given the bus name*/
	int Circuit::find_bus(std::string bus_name) {
		bool found = false;
		for (uint i = 0; i < buses.size(); i++)
			if (!bus_name.compare(buses[i].Name))
				return i;

		return -1;
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
		int n = Circuit::buses.size();

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
				if (lines[i].values_in_per_unit == false) {
					//do the per unit base change
					vbase = buses[lines[i].bus1].nominal_voltage;
					lines[i].Zbase = (vbase * vbase) / Sbase;
				}

				lines[i].get_element_Y(n, Y);
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
	void Circuit::compile(bool guess_angles) {
		// Make sure we do not double-add:

		loadBusIndices.clear();
		slackBusIndices.clear();
		generatorBusIndices.clear();

		//Calculate the base power in a comparable scale to the power in the grid
		double maxpower = get_max_power();
		//cout << "max power:" << maxpower << endl;
		Sbase = pow(10, 1 + floor(log10(maxpower)));

		//Calculate the bus types

		//the presence of an external grid makes the bus VD (or slack)
		for (uint i = 0; i < externalGrids.size(); ++i) {
			if (buses[externalGrids[i].bus].Type == undefined_bus_type) {
				buses[externalGrids[i].bus].Type = VD;
			}
		}

		/*the presence of an generator makes the bus PV (if it is voltage
		 * controlled) or PQ otherwise
		 */
		for (uint i = 0; i < generators.size(); i++) {
			if (buses[generators[i].bus].Type == undefined_bus_type) {
				if (generators[i].voltage_controlled) {
					/* if the generator is set to contol the voltage
					 * then, make the bus a PV bus, otherwise leve the bus fo the
					 * last loop to make it a PQ bus
					 */
					buses[generators[i].bus].Type = PV;
				}
			}
			//apply generator limits conversion to the buses
			buses[generators[i].bus].min_q = generators[i].min_Q / Sbase;
			buses[generators[i].bus].max_q = generators[i].max_Q / Sbase;
			if (generators[i].Vset_in_per_unit)
				buses[generators[i].bus].v_set_point = generators[i].voltage_set_point;
			else
				buses[generators[i].bus].v_set_point = generators[i].voltage_set_point / buses[generators[i].bus].nominal_voltage;
		}

		//the presence of a load makes the bus PQ if it has not ben classified
		for (uint i = 0; i < loads.size(); i++)
			if (buses[loads[i].bus].Type == undefined_bus_type) {
				buses[loads[i].bus].Type = PQ;
				//PQ_list.push_back(buses[loads[i].bus].index);
			}

		//not clasified buses are set to PQ mode
		Vbase = 0.0;
		for (uint i = 0; i < buses.size(); i++) {//final check
			if (buses[i].Type == undefined_bus_type) {
				buses[i].Type = PQ;
			}

			if (buses[i].nominal_voltage > Vbase)
				Vbase = buses[i].nominal_voltage; //set Vbase as the biggest voltage

			//set the bus types lists
			switch (buses[i].Type) {
			case VD:
				slackBusIndices.push_back(buses[i].index);
				break;
			case PQ:
				loadBusIndices.push_back(buses[i].index);
				break;
			case PV:
				generatorBusIndices.push_back(buses[i].index);
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

		if (guess_angles) {
			//Calculate the power connected to each bus according to the type (this time with the new voltage)
			correct_initial_solution();
		}
	}

	/*
	 * Examines the load and generatio objects and returns the maximum value
	 * of the load or generation power in absolute value
	 */
	double Circuit::get_max_power() {
		double mx = 0;

		//Calculate the bus connected generation and load
		for (uint i = 0; i < generators.size(); i++)
			if (abs(generators[i].power) > mx)
				mx = abs(generators[i].power.real());

		for (uint i = 0; i < loads.size(); i++)
			if (abs(loads[i].power) > mx)
				mx = abs(loads[i].power.real());

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

		//Initialize all powers
		for (uint i = 0; i < buses.size(); i++) {
			buses[i].connected_power = cx_double(0.0, 0.0);
		}

		//Calculate the bus connected generation and load
		for (uint i = 0; i < generators.size(); i++)
			buses[generators[i].bus].connected_power += generators[i].power;

		for (uint i = 0; i < loads.size(); i++)
			buses[loads[i].bus].connected_power -= loads[i].power;

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
		std::cout << "Initial D:\n" << dc_angles << std::endl;

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
			lines[i].calculate_current(sol_);
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
		for (uint i = 0; i < sol_.Length; i++)
			sol_.S(i) *= s_base;

		//Undo the voltage change and assign the buses power and voltage
		for (uint i = 0; i < sol_.Length; i++) {
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

		std::cout << "Sbase = " << Sbase << std::endl;

		for (uint i = 0; i < lines.size(); i++)
			lines[i].print();

		for (uint i = 0; i < transformers.size(); i++)
			transformers[i].print();

		for (uint i = 0; i < shunts.size(); i++)
			shunts[i].print();

		for (uint i = 0; i < buses.size(); i++)
			buses[i].print();
	}

	/**/
	void Circuit::print_buses_state() {
		std::string type;
		for (uint i = 0; i < buses.size(); i++) {
			type = BusType_name[buses[i].Type];
			std::cout << i << ": " + type << std::endl;
		}
	}

	/**/
	void Circuit::printCXMat(cx_mat m, std::string header) {
		std::cout << header << std::endl;
		for (uint i = 0; i < buses.size(); i++)
			for (uint j = 0; j < buses.size(); j++)
				if (m.coeff(i, j) != cx_double(0, 0))
					std::cout << "(" << i << "," << j << ") = " << m.coeff(i, j) << std::endl;

	}

	void Circuit::printMat(sp_mat m, std::string header) {
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

	/*
	 * This function sets the load and generation power values (in actual values not in p.u) and updates the solution objects
	 */
	/*void Circuit::setPowerValues(double loadP[], double loadQ[], double genP[], double genQ[]) {

		if (sizeof (loadP) != sizeof (loadQ))
			throw std::invalid_argument("setPowerValues: The size of the load vectors are different");

		if (sizeof (loadP) != loads.size())
			throw std::invalid_argument("setPowerValues: The size of the load vectors does not match the loads size");

		if (sizeof (genP) != sizeof (genQ))
			throw std::invalid_argument("setPowerValues: The size of the generator vectors are different");

		if (sizeof (genP) != generators.size())
			throw std::invalid_argument("setPowerValues: The size of the load vectors does not match the loads size");

		//Set the load values
		for (uint i = 0; i < loads.size(); i++)
			loads[i].power = cx_double(loadP[i], loadQ[i]);

		//Set the generation values
		for (uint i = 0; i < generators.size(); i++)
			generators[i].power = cx_double(genP[i], genQ[i]);

		//Update the solution but keeping the previous values of voltage (and Q, D in case of PV buses)
		generate_initial_solution(true);
	}*/
}
