#include <gtest/gtest.h>

#include "fpotencia.h"
#include "SolverTest.h"
#include "ZBusGSSolverTest.h"


TEST_F(ZBusGSSolverTest, ConvergesOnIeee14BusTest)
{
	auto model = generateIeee14Bus();

	bool estimate_angles = false;
	model.compile(estimate_angles);

	fPotencia::Solver_ZBusGS NRcs(model);
	NRcs.Max_Iter = 6;
	NRcs.EPS = 1e-9;

	auto state = NRcs.solve();

	ASSERT_EQ(fPotencia::Solver_State::Converged, state);
}


TEST_F(ZBusGSSolverTest, SolvesLynnPowellBusWithoutGeneratorTest)
{
	auto model = generateLynnPowellWithoutGenerator();
	model.compile(false);

	fPotencia::Solver_ZBusGS solver(model);
	solver.EPS = 1e-12;
	solver.Max_Iter = 100;
	auto state = solver.solve();

	ASSERT_EQ(fPotencia::Solver_State::Converged, state);

	for (auto& bus: solver.Model.buses) {
	bus.print();
	}

	static const double maxError = 1e-4;

	ASSERT_NEAR(
		1.0,
		solver.Model.buses.at(0).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.954483,
		solver.Model.buses.at(1).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.954023,
		solver.Model.buses.at(2).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.931462,
		solver.Model.buses.at(3).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.952366,
		solver.Model.buses.at(4).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.0,
		solver.Model.buses.at(0).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.040076,
		solver.Model.buses.at(1).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.039373,
		solver.Model.buses.at(2).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.059405,
		solver.Model.buses.at(3).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.044717,
		solver.Model.buses.at(4).voltage_pu.imag(),
		maxError);

	ASSERT_NEAR(
		159.8,
		solver.Model.buses.at(0).power.real(),
		0.1);
}


TEST_F(ZBusGSSolverTest, SolvesLynnPowellBusWithGeneratorTest)
{
	auto model = generateLynnPowellWithGenerator();
	model.compile(false);

	fPotencia::Solver_ZBusGS solver(model);
	solver.EPS = 1e-12;
	solver.Max_Iter = 100;
	auto state = solver.solve();

	ASSERT_EQ(fPotencia::Solver_State::Converged, state);

	for (auto& bus: solver.Model.buses) {
	bus.print();
	}


	static const double maxError = 1e-4;

	ASSERT_NEAR(
		1.0,
		solver.Model.buses.at(0).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.974707,
		solver.Model.buses.at(1).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.981296,
		solver.Model.buses.at(2).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.999916,
		solver.Model.buses.at(3).voltage_pu.real(),
		maxError);
	ASSERT_NEAR(
		0.980726,
		solver.Model.buses.at(4).voltage_pu.real(),
		maxError);

	ASSERT_NEAR(
		0.0,
		solver.Model.buses.at(0).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.026648,
		solver.Model.buses.at(1).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.021173,
		solver.Model.buses.at(2).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.012996,
		solver.Model.buses.at(3).voltage_pu.imag(),
		maxError);
	ASSERT_NEAR(
		-0.025049,
		solver.Model.buses.at(4).voltage_pu.imag(),
		maxError);

	ASSERT_NEAR(
		86.5,
		solver.Model.buses.at(0).power.real(),
		0.1);
}
