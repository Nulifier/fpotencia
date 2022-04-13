cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 0.8 LANGUAGES CXX)

include(FindDoxygen)
#include(FindPkgConfig)
include(GNUInstallDirs)

find_package(GTest)
#pkg_check_modules(eigen3 REQUIRED eigen3>=3.2.8)

option(BUILD_UNSUPPORTED_SOLVERS
	"Build currently unsupported (i.e., test-failing) Load Flow solvers"
	OFF)

option(WITH_OPENMP
	"Include support for parallelization through OpenMP"
	OFF)

if(WITH_OPENMP)
	include(FindOpenMP)
	find_package(OpenMP REQUIRED)

	if(OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
endif()


configure_file("config.h.in" "config.h")
configure_file("Doxyfile.in" "Doxyfile")

add_subdirectory(src)


if(GTEST_FOUND)
	enable_testing()
	add_subdirectory(test)
endif()


if(DOXYGEN_FOUND)
	add_custom_target(doc
		COMMAND ${DOXYGEN_EXECUTABLE}
		DEPENDS Doxyfile.in)
endif()
