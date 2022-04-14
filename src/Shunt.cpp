/* 
 * File:   Shunt.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 6 de agosto de 2014, 10:06
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Shunt.h"

namespace fPotencia {
	Shunt::Shunt(std::string name, int bus, double R, double X) {
		Name = name;
		impedance = cx_double(R, X);
		bus1 = bus;
	}

	/*
	 * This function calculates the shunt admmitance and returns it in the form
	 * of a sparse matrix to make the circuit admittance atrix construction
	 * staight forward
	 */
	void Shunt::get_element_Y(int n, sp_cx_mat &Yret) {

		//dimension check
		if (bus1 > (n - 1)) {
			std::cout << "Shunt>>" << Name << ": Wrong Y dimension: " << n << std::endl;
			std::cout << Yret << std::endl;
			return;
		}

		//Fill the internal matrix
		Y_element = cx_double(1, 0) / impedance;

		//Fil the circuit admittance matrix value
		Yret.coeffRef(bus1, bus1) += Y_element;

	}

	/*
	 * This function calculates the amount of current going through a shunt
	 * element given a circuit solution
	 */
	void Shunt::calculate_current(cx_solution sol) {
		cx_double voltage;

		voltage = sol.V[bus1];

		current = Y_element*voltage;

		power = voltage * conj(current);

		power_losses = power;

		//cout << "Power at " + Name + ": " << power << endl;
		//cout << "\t Losses: " << power_losses << endl;
	}

	/*
	 * This function prints the shunt element calculated parameters
	 */
	void Shunt::print() {
		std::cout << Name << std::endl;
		std::cout << "\tPower: " << power << std::endl;

		std::cout << "\tCurrent: " << current << std::endl;

		std::cout << "\tLosses: " << power_losses << std::endl;
	}

}