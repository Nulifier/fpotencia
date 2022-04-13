/*
 * File:   Transformer.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Transformer.h"
#include "TransformerConstructors.h"

namespace fPotencia {

	/*******************************************************************************
	 *TransformerType class implementation
	 ******************************************************************************/

	/*
	 * Transformer type object initialization with calculated parameters
	 */
	TransformerType::TransformerType(std::string name, cx_double leakage_z, cx_double magnetizing_z) {
		//Calculated parameters in per unit
		Name = name;
		leakage_impedance = leakage_z;
		magnetizing_impedance = magnetizing_z;
		tap = 1.0;
		phase_shift = PI / 6.0;
	}

	/*
	 * Transformer type object destructor
	 */
	TransformerType::~TransformerType() {
	}



	/***************************************************************************
	 * Transformer class implementation
	 ***************************************************************************/

	/*
	 * Transformer object constructor
	 */
	Transformer::Transformer(std::string name, int connection_busHV, int connection_busLV, TransformerType transformer_type) {
		Name = name;
		HV_bus_index = connection_busHV;
		LV_bus_index = connection_busLV;
		SetType(transformer_type);
	}

	/*
	 * transformer object destructor
	 */
	Transformer::~Transformer() {
	}

	/*
	 * This function calculates the transformer impedance and admittance from a
	 * transformer type model and generates the transformer admmitance matrix
	 */
	void Transformer::SetType(TransformerType transformer_type) {
		tap = transformer_type.tap;
		pha_shift = transformer_type.phase_shift;
		leakage_impedance = transformer_type.leakage_impedance;
		magnetizing_impedance = transformer_type.magnetizing_impedance;

		/*Calculate the transformer admittance matrix*/
		cx_double Yl = cx_double(1, 0) / leakage_impedance;
		cx_double Ym = cx_double(1, 0) / magnetizing_impedance;
		double sin_phase = sin(pha_shift);
		double cos_phase = cos(pha_shift);
		cx_double tap_ = cx_double(tap * cos_phase, tap * sin_phase);

		Y_element = cx_mat(2, 2);
		Y_element(0, 0) = (Yl + Ym) / (tap_ * tap_); //primary-primary    
		Y_element(0, 1) = -(Ym + Yl) / conj(tap_); //primary-secondary
		Y_element(1, 0) = -Yl / cx_double(tap, 0.0); //secondary-primary
		Y_element(1, 1) = Yl; //secondary-secondary

		//check for inf or nan
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				if (std::isinf(Y_element.coeff(i, j).real()) || std::isinf(Y_element.coeff(i, j).imag())) {
					std::cout << Y_element << std::endl;
					std::cout << "Zm:" << magnetizing_impedance << std::endl;
					std::cout << "Zl:" << leakage_impedance << std::endl;
					std::cout << "tap:" << tap << std::endl;
					std::stringstream ss;
					ss << "Transformer>>" << Name << ": infinite or nan values in the element Y at: " << i << "," << j;
					throw std::invalid_argument(ss.str());
				}
	}

	/*
	 * This function modifies the circuit admittance matrix Yret with the
	 * previousy calculated transformer admittance matrix
	 */
	void Transformer::get_element_Y(int n, sp_cx_mat & Yret) {

		//dimension check
		if (HV_bus_index > (n - 1) || LV_bus_index > (n - 1)) {
			std::stringstream ss;
			ss << "Transformer>>" << Name << ": Wrong Y dimension: " << n;
			throw std::invalid_argument(ss.str());
			//std::cout << "Transformer>>" << Name << ": Wrong Y dimension: " << n << endl;
			return;
		}

		//yff
		Yret.coeffRef(HV_bus_index, HV_bus_index) += Y_element.coeff(0, 0); //cx_double(tap*tap, 0); //primary-primary        
		//yft
		Yret.coeffRef(HV_bus_index, LV_bus_index) += Y_element.coeff(0, 1); //primary-secondary
		//ytf
		Yret.coeffRef(LV_bus_index, HV_bus_index) += Y_element.coeff(1, 0); //secondary-primary
		//ytt
		Yret.coeffRef(LV_bus_index, LV_bus_index) += Y_element.coeff(1, 1); //secondary-secondary
	}

	/*
	 * This function calculates the current and powers flowing into the 
	 * transformer given a circuit solution
	 */
	void Transformer::calculate_current(cx_solution sol) {
		cx_mat voltage(2, 1);
		cx_mat current(2, 1);
		cx_mat power(2, 1);

		voltage(0, 0) = sol.V[HV_bus_index];
		voltage(1, 0) = sol.V[LV_bus_index];

		current = Y_element*voltage;

		power(0, 0) = voltage(0, 0) * conj(current(0, 0));
		power(1, 0) = voltage(1, 0) * conj(current(1, 0));

		current_primary_to_secondary = current(0, 0);
		current_secondary_to_primary = current(1, 0);

		power_primary_to_secondary = power(0, 0);
		power_secondary_to_primary = power(1, 0);

		power_losses = power_primary_to_secondary + power_secondary_to_primary;

		//cout << "Power at " + Name + "\n" << power << endl;
		//cout << "\t Losses: " << power_losses << endl;
	}

	/*
	 * This function  prints all the calculated values of the transformer
	 */
	void Transformer::print() {
		std::cout << Name << std::endl;
		std::cout << "\tPower" << std::endl;
		std::cout << "\t bus 1 to 2: " << power_primary_to_secondary << std::endl;
		std::cout << "\t bus 2 to 1: " << power_secondary_to_primary << std::endl;

		std::cout << "\t Losses: " << power_losses << std::endl;

		std::cout << "\tCurrent" << std::endl;
		std::cout << "\t bus 1 to 2: " << current_primary_to_secondary << std::endl;
		std::cout << "\t bus 2 to 1: " << current_secondary_to_primary << std::endl;
	}
}
