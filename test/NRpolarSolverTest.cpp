#include <gtest/gtest.h>

#include "fpotencia.h"
#include "SolverTest.h"
#include "NRpolarSolverTest.h"


TEST_F(NRpolarSolverTest, ConvergesOnIeee14BusTest)
{
	auto model = generateIeee14Bus();

	bool estimate_angles = false;
	model.compile(estimate_angles);

	fPotencia::NRpolarSolver NRcs(model);
	NRcs.maxIterations = 6;
	NRcs.tolerance = 1e-9;

	auto state = NRcs.powerFlow(model);

	ASSERT_EQ(fPotencia::Solver::Solved, state);
}


TEST_F(NRpolarSolverTest, SolvesLynnPowellBusWithoutGeneratorTest)
{
	auto model = generateLynnPowellWithoutGenerator();
	model.compile(false);

	fPotencia::NRpolarSolver solver(model);
	solver.tolerance = 1e-12;
	solver.maxIterations = 100;
	auto state = solver.powerFlow(model);

	ASSERT_EQ(fPotencia::Solver::Solved, state);


	static const double maxError = 1e-4;

	ASSERT_NEAR(
			1.0,
			model.buses.at(0).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.954483,
			model.buses.at(1).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.954023,
			model.buses.at(2).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.931462,
			model.buses.at(3).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.952366,
			model.buses.at(4).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.0,
			model.buses.at(0).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.040076,
			model.buses.at(1).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.039373,
			model.buses.at(2).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.059405,
			model.buses.at(3).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.044717,
			model.buses.at(4).voltage_pu.imag(),
			maxError);

	ASSERT_NEAR(
			159.8,
			model.buses.at(0).power.real(),
			0.1);
}


TEST_F(NRpolarSolverTest, SolvesLynnPowellBusWithGeneratorTest)
{
	auto model = generateLynnPowellWithGenerator();
	model.compile(false);

	fPotencia::NRpolarSolver solver(model);
	solver.tolerance = 1e-12;
	solver.maxIterations = 100;
	auto state = solver.powerFlow(model);

	ASSERT_EQ(fPotencia::Solver::Solved, state);

	static const double maxError = 1e-4;

	ASSERT_NEAR(
			1.0,
			model.buses.at(0).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.974707,
			model.buses.at(1).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.981296,
			model.buses.at(2).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.999916,
			model.buses.at(3).voltage_pu.real(),
			maxError);
	ASSERT_NEAR(
			0.980726,
			model.buses.at(4).voltage_pu.real(),
			maxError);

	ASSERT_NEAR(
			0.0,
			model.buses.at(0).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.026648,
			model.buses.at(1).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.021173,
			model.buses.at(2).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.012996,
			model.buses.at(3).voltage_pu.imag(),
			maxError);
	ASSERT_NEAR(
			-0.025049,
			model.buses.at(4).voltage_pu.imag(),
			maxError);

	ASSERT_NEAR(
			86.5,
			model.buses.at(0).power.real(),
			0.1);
}
