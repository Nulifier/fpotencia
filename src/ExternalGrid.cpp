/* 
 * File:   ExternalGrid.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 8 de agosto de 2014, 14:45
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "ExternalGrid.h"

namespace fPotencia {
	/*
	 * External grid object constructor
	 */
	ExternalGrid::ExternalGrid(std::string name, int connection_bus) {
		Name = name;
		bus = connection_bus;
	}
}
