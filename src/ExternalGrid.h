/* 
 * File:   ExternalGrid.h
 * Author: Santiago Peñate Vera
 *
 * Created on 8 de agosto de 2014, 14:45
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#include "fpotencia_libs.h"

namespace fPotencia {
	class ExternalGrid final {
		public:
		ExternalGrid(std::string name, int connection_bus);

		std::string Name;

		/** The bus where this ExternalGrid is connected */
		int bus;

		/** The calculated power at this node */
		cx_double power;
	};
}
