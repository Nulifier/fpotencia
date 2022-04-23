/* 
 * File:   Load.h
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

#include "fpotencia_libs.h"

namespace fPotencia {
	/** A load connected to a bus */
	class Load final {
	public:
		Load(const std::string& name, int connection_bus, double P, double Q);

		std::string Name;

		/** The bus this load is connected to */
		int bus;

		/* The power consumed by this load */
		cx_double power;
	};
}
