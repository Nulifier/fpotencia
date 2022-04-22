#pragma once

#include <Circuit.h>

class TestModels {
public:
	static fPotencia::Circuit ieee14Model();
	static fPotencia::Circuit lynnPowerllWithGenerator();
	static fPotencia::Circuit lynnPowerllWithoutGenerator();
	static fPotencia::Circuit ieee300Model();
};
