/* 
 * File:   Shunt.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:06
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "fpotencia_libs.h"
#include "Solution.h"

namespace fPotencia {
	class Shunt final {
	public:                
	    Shunt(std::string name, int bus, double R, double X);

		void get_element_Y(int n, sp_cx_mat &Yret);

		void calculate_current(cx_solution sol);

		void print() const;

		/*************************************************************************
		 * Properties
		 *************************************************************************/
		std::string Name;

		int bus1 = 0;

		cx_double impedance;

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double current;

		cx_double power;

		cx_double power_losses;

		private:

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double Y_element;

	};
}
