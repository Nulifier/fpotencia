/*
 * File:   Solver_def.h
 * Author: Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#pragma once

#define PI 3.1415926535897932384626433832795

//Define unit conversions
#define mile2km 1.609344
#define ft2in   12.0
#define ft2m    0.305
#define in2m    0.025

#define m2ft 3.2808399
#define m2in 39.3700787
#define km2mile 0.621371192

namespace fPotencia {

	enum Solver_State {
		Converged,
		Not_Converged,
		Not_Solvable,
		Not_Solvable_with_Method
	};

	enum BusType {
		PQ,
		PV,
		VD, //Same as slack 
		undefined_bus_type
	};

	enum Cable_Type {
		ConcentricNeutral,
		TapeShield
	};

	enum Conductor_Type {
		Cable,
		Tape_Shield,
		Void
	};

	enum ConnectionPhase {
		A,
		B,
		C
	};

	enum TransformerConnexionType {
		Yg_Yg,
		Yg_Y,
		Yg_D,
		Y_Yg,
		Y_Y,
		Y_D,
		D_Yg,
		D_Y,
		D_D
	};

	enum Initialization_Mode {
		Positive_Sequence,
		Three_Phase
	};
	
	enum Units_mode{
		US, //meaning a weird mix of metric and imperial, read where aplicable
		Metric //metric system compatible m, km, etc.
	};

	static const char * BusType_name[] = {
		"PQ",
		"PV",
		"VD",
		"undefined bus type"
	};
}
