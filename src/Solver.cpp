#include "Solver.h"

namespace fPotencia {
	Solver::Solver():
			tolerance(DEFAULT_SOLUTION_TOLERANCE),
			maxIterations(DEFAULT_MAX_ITERATIONS)
	{
	}


	Solver::~Solver() noexcept
	{
	}
} // namespace fPotencia
