#pragma once

#include <gtest/gtest.h>
#include "Circuit.h"

class SolverTest: public ::testing::Test
{
public:
	SolverTest();
	virtual ~SolverTest() noexcept;
	fPotencia::Circuit generateIeee14Bus() const;
	fPotencia::Circuit generateLynnPowellWithGenerator() const;
	fPotencia::Circuit generateLynnPowellWithoutGenerator() const;
};
