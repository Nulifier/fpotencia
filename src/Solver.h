#pragma once

namespace fPotencia {
	class Circuit;

	class Solver {
	public:
		virtual ~Solver() = default;

		//! \brief Result indicator flag for the #powerFlow() method
		enum Result
		{
			Solved,
			NotSolved,
			NotSolveable,
			SolverUnfitForGrid
		};

		/*!
		 * \brief Allowable tolerance of the solver instance
		 *
		 * Solving the power flow equations is an interative process. For
		 * every iteration, the parameters are adjusted in order to reach
		 * convergence. The tolerance defines the allowable deviation/error
		 * after which the process halts.
		 */
		double tolerance = 1e-9;

		/*!
		 * \brief Default maximum number of iterations after which the solver
		 *  declares failure
		 */
		unsigned int maxIterations = 100;

		/*!
		 * Calculates the power flow in the given circuit.
		 *
		 * This is the main method for any solver: It calculates the power
		 * flow using the solver's specific method. Each individual solver
		 * implements it.
		 *
		 * The method might throw an UnfittingSolverException if the solver
		 * selected cannot work on the given grid. Details on this can be
		 * found in the individual solver's documentation.
		 *
		 * \return A flag indicating success or failure
		 *
		 * \throw UnfittingSolverError
		 *
		 * \sa Solver::Result
		 */
		virtual Result solve(bool printIterations = false) = 0;

		/**
		 * Checks if this solver can solve this model.
		 * @note The solver can still fail to solve the model, this just runs
		 *       checks that can be done before.
		 */
		virtual bool canSolve() const = 0;
	};
}
