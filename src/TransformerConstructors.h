/* 
 * File:   LineConstructors.h
 * Author: Santiago Peñate Vera
 *
 * Created on 15 de March de 2015, 10:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "Solution.h"
#include "fpotencia_libs.h"

/*
 * This file contains the neccesary classes to model 3-phase overhead lines
 * and cables for distribution.
 * 
 * the refference book is:
 * Distribution systems modelling and analisys (3rd Edition)
 * by William H Kersting, CRC press.
 * 
 * The magnitudes are as used in the book, that means a quite heretic mix of 
 * units like they like in the US.
 * 
 * Units converters from and to the international metric system are provided
 */

namespace fPotencia {
	/*
	 * 
	 */
	class Transformer_3phase_Yabc {
	public:
		Transformer_3phase_Yabc();

		virtual ~Transformer_3phase_Yabc();

		cx_mat Yabc(cx_double leakage_impedance,
				TransformerConnexionType connexionType,
				double tap_hv,
				double tap_lv);
	};

	/* 
	 * Computes the transformer parameters from its short circuit test values
	 */
	class Transformer_ShortCircuit_Constructor {
	public:


		Transformer_ShortCircuit_Constructor(std::string name,
				TransformerConnexionType connexionType,
				double HV_nominal_voltage,
				double LV_nominal_voltage, double Nominal_power,
				double Copper_losses, double Iron_losses,
				double No_load_current, double Short_circuit_voltage,
				double GX_HV1,
				double GR_HV1);

		/**/
		virtual~Transformer_ShortCircuit_Constructor();

		cx_mat Y_abc; //Calculates admittance matrix from the constructor

		//Calculated parameters in per unit values
		cx_double leakage_impedance; //r + j*x

		cx_double magnetizing_impedance; //rfe + j*xm   

		// Tap angle phase shift [rad]    
		double phase_shift = PI / 6.0;
		
	private:

		void calculate_model();

		// Type name    
		std::string Name;

		// HV side nominal voltage [kV]    
		double Uhv;

		// LV side nominal voltage [kV]    
		double Ulv;

		// Nominal power [MVA]    
		double Sn;

		// Copper losses [kW]    
		double Pcu;

		// Iron losses [kW]    
		double Pfe;

		// No load current [%]    
		double I0;

		// Short circuit voltage [%]    
		double Usc;

		// Reactance contribution to the HV side [from 0 to 1]    
		double GX_hv1;

		// Resistance contribution to the HV side [from 0 to 1]    
		double GR_hv1;

		TransformerConnexionType connType;
	};
}
