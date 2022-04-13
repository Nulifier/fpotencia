/* 
 * File:   Generator.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Generator.h"

namespace fPotencia {

	/*
	 * Constructor for non voltage controlled generators
	 */
	Generator::Generator(std::string name, uint connection_bus, double P, double Q) {
		Name = name;
		bus = connection_bus;
		power = cx_double(P, Q);
		Vset_in_per_unit = true;
		voltage_set_point = 1.0;
		voltage_controlled = false;
	}

	/*
	 * Class constructor for voltage controlled generators
	 */
	Generator::Generator(std::string name, uint connection_bus, double P, double Vset, double Qmin, double Qmax, bool Vset_per_unit) {
		Name = name;
		bus = connection_bus;
		power = cx_double(P, 0.0);
		voltage_set_point = Vset;
		voltage_controlled = true;
		min_Q = Qmin;
		max_Q = Qmax;
		Vset_in_per_unit = Vset_per_unit;
	}
}
