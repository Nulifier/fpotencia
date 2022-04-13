#pragma once

#include <stdexcept>

/*!
 * \brief The UnfittingSolverException class indicates a logic error when the
 *  wrong solver was chosen to calculate the power flow for a given grid.
 */
class UnfittingSolverException: public std::logic_error {
public:
	UnfittingSolverException(const char* what);
	virtual ~UnfittingSolverException() noexcept;
};
