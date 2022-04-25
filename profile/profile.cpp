#include <iostream>
#include <chrono>
#include <Solver_NRpolar.h>
#include <test-models.h>

int main() {
	const auto start = std::chrono::high_resolution_clock::now();

	auto model = TestModels::ieee300Model();

	model.compile(false);

	fPotencia::NRpolarSolver solver(model);
	solver.maxIterations = 20;
	solver.tolerance = 1e-6;

	auto result = solver.solve(true);

	const auto end = std::chrono::high_resolution_clock::now();

	std::cout << "System did " << (result != fPotencia::Solver::Result::Solved ? "not " : "") << "converge in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

	// Current record in release: 224ms
}
