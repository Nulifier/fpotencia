#include <stdexcept>
#include "UnfittingSolverException.h"


UnfittingSolverException::UnfittingSolverException(const char* what):
	std::logic_error(what)
{
}


UnfittingSolverException::~UnfittingSolverException() noexcept
{
}
