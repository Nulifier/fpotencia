#pragma once

namespace fPotencia {
	class Circuit;

	class Solver {
	public:
		//! \brief Result indicator flag for the #powerFlow() method
		enum Result
		{
			Solved,
			NotSolved,
			NotSolveable,
			SolverUnfitForGrid
		};

		/*!
		 * \brief The default tolerance for the solution
		 *
		 * Solving the power flow equations is an interative process. For
		 * every iteration, the parameters are adjusted in order to reach
		 * convergence. The tolerance defines the allowable deviation/error
		 * after which the process halts.
		 */
		static constexpr const double DEFAULT_SOLUTION_TOLERANCE = 1e-9;

		/*!
		 * \brief Default maximum number of iterations after which the solver
		 *  declares failure
		 */
		static constexpr const unsigned DEFAULT_MAX_ITERATIONS = 100;


		/*!
		 * \brief Allowable tolerance of the solver instance
		 *
		 * \sa DEFAULT_SOLUTION_TOLERANCE
		 */
		double tolerance = DEFAULT_SOLUTION_TOLERANCE;


		/*!
		 * \brief Maximum number of iterations
		 *
		 * \sa DEFAULT_MAX_ITERATIONS
		 */
		unsigned int maxIterations = DEFAULT_MAX_ITERATIONS;

		/*!
		 * \brief Calculates the power flow in the given circuit
		 *
		 * This is the main method for any solver: It calculates the power
		 * flow using the solver's specific method. Each individual solver
		 * implements it.
		 *
		 * The method might throw an UnfittingSolverException if the solver
		 * selected cannot work on the given grid. Details on this can be
		 * found in the individual solver's documentation.
		 *
		 * \param[inout] grid The given circuit
		 *
		 * \return A flag indicating success or failure
		 *
		 * \throw UnfittingSolverError
		 *
		 * \sa Solver::Result
		 */
		virtual Result powerFlow(Circuit& grid) = 0;
	};
}
