/* 
 * File:   Bus.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 6 de agosto de 2014, 9:55
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Bus.h"

namespace fPotencia {
	Bus::Bus(std::string name, BusType type, double Bus_Nominal_Voltage) {
		//Name = name;
		Type = type;
		nominal_voltage = Bus_Nominal_Voltage;
	}

	void Bus::print() {
		//std::cout << Name + " -> " + BusType_name[Type] << std::endl;
		std::cout << "\tPower: " << power << std::endl;

		std::cout << "\tVoltage: " << voltage << "\tVoltage p.u.: " << voltage_pu << std::endl;
	}
}
