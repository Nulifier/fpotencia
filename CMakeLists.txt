cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 0.8 LANGUAGES CXX)

include(FindDoxygen)
#include(FindPkgConfig)
include(GNUInstallDirs)

find_package(GTest)
find_package(Eigen3)
#pkg_check_modules(eigen3 REQUIRED eigen3>=3.2.8)

option(BUILD_UNSUPPORTED_SOLVERS
	"Build currently unsupported (i.e., test-failing) Load Flow solvers"
	OFF)

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
