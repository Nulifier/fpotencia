/*
 * File:   Transformer.h
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

#include <math.h>
#include "fpotencia_libs.h"
#include "Solution.h"
#include "Bus.h"

namespace fPotencia {
	/*******************************************************************************
	 *TransformerType class definition
	 ******************************************************************************/
	class TransformerType final {
	public:

		TransformerType(std::string name, cx_double leakage_z, cx_double magnetizing_z);

		/*
		 * Name: Name of the type
		 * Yabc: Transformer admittance matrix -> can be created with 
		 *       TransformerConstructors.h
		 */
		TransformerType(std::string name, cx_mat Yabc);

		// Type name    
		std::string Name;

		// Tap position [Value arround 1]    
		double tap = 1.0;

		double phase_shift = PI / 6.0; // 30 deg by default (always in radians; 30 deg = pi/6 rad)

		//Calculated parameters in per unit values
		cx_double leakage_impedance; //r + j*x

		cx_double magnetizing_impedance; //rfe + j*xm
	};

	/*******************************************************************************
	 *Transformer class definition
	 ******************************************************************************/
	class Transformer {
	public:

		Transformer(std::string name, int connection_busHV, int connection_busLV, TransformerType transformer_type);

		virtual ~Transformer();

		void SetType(TransformerType transformer_type);

		void get_element_Y(int n, sp_cx_mat &Yret);

		void calculate_current(cx_solution sol);

		void print();

		//properties
		std::string Name;

		int HV_bus_index = 0;

		int LV_bus_index = 0;

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double current_primary_to_secondary;

		cx_double current_secondary_to_primary;

		cx_double power_primary_to_secondary;

		cx_double power_secondary_to_primary;

		cx_double power_losses;

	private:

		double tap;

		double pha_shift;

		cx_double leakage_impedance; //r + j*x

		cx_double magnetizing_impedance; //rfe + j*xm

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_mat Y_element;
	};
}
