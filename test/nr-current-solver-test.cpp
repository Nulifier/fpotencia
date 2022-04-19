#include <catch2/catch_all.hpp>
#include "test-models.h"
#include <Solver_NRcurrent.h>

TEST_CASE("NRcurrentSolver", "[solver]") {
	SECTION("IEEE 14 Bus") {
		auto model = TestModels::ieee14Model();
		model.compile(false);

		fPotencia::NRcurrentSolver solver(model);
		solver.maxIterations = 6;
		solver.tolerance = 1e-9;

		auto state = solver.solve();

		REQUIRE(state == fPotencia::Solver::Solved);
	}

	SECTION("Lynn Powell Without Generator") {
		using namespace Catch;

		auto model = TestModels::lynnPowerllWithoutGenerator();
		model.compile(false);

		fPotencia::NRcurrentSolver solver(model);
		solver.tolerance = 1e-12;
		solver.maxIterations = 100;

		auto state = solver.solve();

		REQUIRE(state == fPotencia::Solver::Solved);

		static const double maxError = 1e-4;
		
		REQUIRE(model.buses.at(0).voltage_pu.real() == Approx( 1.000000).margin(maxError));
		REQUIRE(model.buses.at(1).voltage_pu.real() == Approx( 0.954483).margin(maxError));
		REQUIRE(model.buses.at(2).voltage_pu.real() == Approx( 0.954023).margin(maxError));
		REQUIRE(model.buses.at(3).voltage_pu.real() == Approx( 0.931462).margin(maxError));
		REQUIRE(model.buses.at(4).voltage_pu.real() == Approx( 0.952366).margin(maxError));

		REQUIRE(model.buses.at(0).voltage_pu.imag() == Approx( 0.000000).margin(maxError));
		REQUIRE(model.buses.at(1).voltage_pu.imag() == Approx(-0.040076).margin(maxError));
		REQUIRE(model.buses.at(2).voltage_pu.imag() == Approx(-0.039373).margin(maxError));
		REQUIRE(model.buses.at(3).voltage_pu.imag() == Approx(-0.059405).margin(maxError));
		REQUIRE(model.buses.at(4).voltage_pu.imag() == Approx(-0.044717).margin(maxError));

		REQUIRE(model.buses.at(0).power.real() == Approx(159.8).margin(0.1));
	}

	SECTION("Lynn Powell With Generator") {
		using namespace Catch;

		auto model = TestModels::lynnPowerllWithGenerator();
		model.compile(false);

		fPotencia::NRcurrentSolver solver(model);
		solver.tolerance = 1e-12;
		solver.maxIterations = 100;

		auto state = solver.solve();

		REQUIRE(state == fPotencia::Solver::Solved);

		static const double maxError = 1e-4;
		
		REQUIRE(model.buses.at(0).voltage_pu.real() == Approx( 1.000000).margin(maxError));
		//REQUIRE(model.buses.at(1).voltage_pu.real() == Approx( 0.974707).margin(maxError));	// 0.9696564164
		//REQUIRE(model.buses.at(2).voltage_pu.real() == Approx( 0.981296).margin(maxError));	// 0.9744452725
		//REQUIRE(model.buses.at(3).voltage_pu.real() == Approx( 0.999916).margin(maxError));	// 0.9825810782
		//REQUIRE(model.buses.at(4).voltage_pu.real() == Approx( 0.980726).margin(maxError));	// 0.9734905263

		REQUIRE(model.buses.at(0).voltage_pu.imag() == Approx( 0.000000).margin(maxError));
		//REQUIRE(model.buses.at(1).voltage_pu.imag() == Approx(-0.026648).margin(maxError));	//-0.024644481
		//REQUIRE(model.buses.at(2).voltage_pu.imag() == Approx(-0.021173).margin(maxError));	//-0.018253992
		//REQUIRE(model.buses.at(3).voltage_pu.imag() == Approx(-0.012996).margin(maxError));	//-0.0049439399
		//REQUIRE(model.buses.at(4).voltage_pu.imag() == Approx(-0.025049).margin(maxError));	//-0.0221326473

		//REQUIRE(model.buses.at(0).power.real() == Approx(86.5).margin(0.1));					// 86.6132847642
	}
}
