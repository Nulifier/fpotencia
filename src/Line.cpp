/* 
 * File:   Line.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Line.h"

namespace fPotencia {
	LineType::LineType(const std::string& name, double r, double x, double b, bool per_unit_values) {
		Name = name;

		if (b == 0)
			b = 1e-9;

		impedance = cx_double(r, x);
		shunt_admittance = cx_double(0.0, b);
		values_in_per_unit = per_unit_values;
	}

	Line::Line(const std::string& name, int connection_bus1, int connection_bus2, LineType line_type, double line_length) {
		Name = name;
		bus1 = connection_bus1;
		bus2 = connection_bus2;
		length = line_length;

		SetType(line_type);
	}

	/*
	 * This function calculates the line impedance and admittance from a line type
	 * model
	 */
	void Line::SetType(LineType line_type) {
		shunt_admittance = line_type.shunt_admittance * length;
		impedance = line_type.impedance * length;
		values_in_per_unit = line_type.values_in_per_unit;

		//create the element admittance matrix
		cx_double y = cx_double(1, 0) / impedance;
		cx_double ys = shunt_admittance / cx_double(2, 0);

		//Fill the internal matrix
		Y_element = cx_mat(2, 2);
		Y_element(0, 0) = y + ys;
		Y_element(0, 1) = -y;
		Y_element(1, 1) = y + ys;
		Y_element(1, 0) = -y;

		//check for inf or nan
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				if (std::isinf(Y_element.coeff(i, j).real()) || std::isinf(Y_element.coeff(i, j).imag())) {
					std::cout << Y_element << std::endl;
					std::stringstream ss;
					ss << "Line>>" << Name << ": infinite or nan values in the element Y at: " << i << "," << j;
					throw std::invalid_argument(ss.str());
				}
	}

	void Line::get_element_Y(int n, sp_cx_mat &Yret) {
		//dimension check
		if (bus1 > (n - 1) || bus2 > (n - 1)) {
			std::stringstream ss;
			ss << "Line>>" << Name << ": Wrong Y dimension: " << n;
			throw std::invalid_argument(ss.str());
			return;
		}

		//set the circuit matrix values
		if (values_in_per_unit) {
			Yret.coeffRef(bus1, bus1) += Y_element.coeff(0, 0);
			Yret.coeffRef(bus1, bus2) += Y_element.coeff(0, 1);
			Yret.coeffRef(bus2, bus2) += Y_element.coeff(1, 1);
			Yret.coeffRef(bus2, bus1) += Y_element.coeff(1, 0);
		} else {
			Yret.coeffRef(bus1, bus1) += Y_element.coeff(0, 0) * Zbase;
			Yret.coeffRef(bus1, bus2) += Y_element.coeff(0, 1) * Zbase;
			Yret.coeffRef(bus2, bus2) += Y_element.coeff(1, 1) * Zbase;
			Yret.coeffRef(bus2, bus1) += Y_element.coeff(1, 0) * Zbase;
		}
	}

	void Line::calculate_current(cx_solution sol) {
		cx_mat voltage(2, 1);
		cx_mat current(2, 1);
		cx_mat power(2, 1);

		voltage(0, 0) = sol.V[bus1];
		voltage(1, 0) = sol.V[bus2];

		current = Y_element * voltage;

		power(0, 0) = voltage(0, 0) * conj(current(0, 0));
		power(1, 0) = voltage(1, 0) * conj(current(1, 0));

		current_bus1_to_bus2 = current(0, 0);
		current_bus2_to_bus1 = current(1, 0);

		power_bus1_to_bus2 = power(0, 0);
		power_bus2_to_bus1 = power(1, 0);

		power_losses = power_bus1_to_bus2 + power_bus2_to_bus1;
		/*if (power_bus1_to_bus2.real() > power_bus2_to_bus1.real())
				power_losses = power_bus1_to_bus2 + power_bus2_to_bus1;
				else
				power_losses = power_bus2_to_bus1 - power_bus1_to_bus2;*/
	}

	void Line::print() {
		std::cout << Name << std::endl;
		std::cout << "\t r:" << impedance.real() << ", x:" << impedance.imag() << ", c: " << shunt_admittance.imag() << std::endl;
		std::cout << "\tPower" << std::endl;
		std::cout << "\t bus 1 to 2: " << power_bus1_to_bus2 << std::endl;
		std::cout << "\t bus 2 to 1: " << power_bus2_to_bus1 << std::endl;

		std::cout << "\t Losses: " << power_losses << std::endl;

		std::cout << "\tCurrent" << std::endl;
		std::cout << "\t bus 1 to 2: " << current_bus1_to_bus2 << std::endl;
		std::cout << "\t bus 2 to 1: " << current_bus2_to_bus1 << std::endl;
	}
}
