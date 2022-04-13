/* 
 * File:   Line.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "Solution.h"
#include "fpotencia_libs.h"

namespace fPotencia {
	/*******************************************************************************
	 *LineType class definition
	 ******************************************************************************/
	class LineType final {
	public:
		LineType(std::string name, double r, double x, double b, bool per_unit_values);

		std::string Name;

		cx_double impedance;

		cx_double shunt_admittance;

		bool values_in_per_unit = false;

	private:

		cx_mat Y;
	};

	/*******************************************************************************
	 *Line class definition
	 ******************************************************************************/
	class Line final {
	public:
		Line(std::string name, int connection_bus1, int connection_bus2, LineType line_type, double line_lenght);

		void SetType(LineType line_type);

		std::string Name;

		int bus1 = 0;

		int bus2 = 0;

		double lenght = 0;

		bool values_in_per_unit;

		void get_element_Y(int n, sp_cx_mat &Yret);

		void calculate_current(cx_solution sol);

		void print();


		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double current_bus1_to_bus2;

		cx_double current_bus2_to_bus1;

		cx_double power_bus1_to_bus2;

		cx_double power_bus2_to_bus1;

		cx_double power_losses;

		double Zbase;

	private:

		cx_double impedance;

		cx_double shunt_admittance;

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/
		cx_mat Y_element; //calculated element admittance matrix (2x2)
	};
}
