/*
 * File:   Load.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Load.h"

namespace fPotencia {
	Load::Load(const std::string& name, int connection_bus, double P, double Q)
	    : Name(name), bus(connection_bus) {

		power = cx_double(P, Q);
	}
}