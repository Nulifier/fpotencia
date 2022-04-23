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
	Generator::Generator(const std::string& name, uint connection_bus, double P,
	                     double Q)
	    : Name(name), bus(connection_bus), voltage_set_point(1.0),
	      voltage_controlled(false), Vset_in_per_unit(true), min_Q(0),
	      max_Q(0) {

		power = cx_double(P, Q);
	}

	Generator::Generator(const std::string& name, uint connection_bus, double P,
	                     double Vset, double Qmin, double Qmax,
	                     bool Vset_per_unit)
	    : Name(name), bus(connection_bus), min_Q(Qmin), max_Q(Qmax),
	      voltage_set_point(Vset), voltage_controlled(true),
	      Vset_in_per_unit(Vset_per_unit) {

		power = cx_double(P, 0.0);
	}
}
