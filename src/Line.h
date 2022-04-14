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
	class LineType final {
	public:
		/**
		 * Constructs a LineType.
		 * @param name The name of this line type.
		 * @param r The resistance per distance unit.
		 * @param x The reactance per distance unit.
		 * @param b The susceptance per distance unit.
		 */
		LineType(const std::string& name, double r, double x, double b, bool per_unit_values);

		std::string Name;

		cx_double impedance;

		cx_double shunt_admittance;

		bool values_in_per_unit = false;
	};

	class Line final {
	public:
		Line(const std::string& name, int connection_bus1, int connection_bus2, LineType line_type, double line_length);

		void SetType(LineType line_type);

		std::string Name;

		/** The bus id that this line connects from */
		int bus1;

		/** The bus id that this line connects to */
		int bus2;

		/** The length of this line in distance units */
		double length;

		bool values_in_per_unit;

		/*
		 * Returns the component admittance matrix in sparse format, this way
		 * the composition of the circuit admittance matrix is straight forward
		 */
		void get_element_Y(int n, sp_cx_mat &Yret);

		/**
		 * This function calculates the amount of current going through the line
	 	 * given a circuit solution
	 	 */
		void calculate_current(cx_solution sol);

		/** Prints all the calculated line parameters. */
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
