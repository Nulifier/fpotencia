/* 
 * File:   Bus.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 9:54
 * 
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "fpotencia_libs.h"

namespace fPotencia {
	class Bus final {
	public:
		Bus(std::string name, BusType type, double Bus_Nominal_Voltage);

		int index = -1;

		//std::string Name;

		BusType Type;

		/** The nominal voltage of this bus (ie. 480V) */
		double nominal_voltage;

		/** The calculated voltage of this bus in volts */
		cx_double voltage;

		/** The calculated voltage of this bus in nominal voltage units */
		cx_double voltage_pu;

		/** The calculated power at this bus */
		cx_double power;

		/**
		 * The sum of the connected loads and generators.
		 * @note This is updated as part of calculating the initial solution.
		 */
		cx_double connected_power = cx_double(0, 0);

		/*-----Only for PV buses--------------------------------------*/

		double min_q = 0; //minimum reactive power per unit (changed in the circuit class)

		double max_q = 0; //maximum reactive power per unit (changed in the circuit class)

		double v_set_point = 1.0; //only used if the bus is PV
		/*------------------------------------------------------------*/

		/** Print the bus results. */
		void print();
	};
}
