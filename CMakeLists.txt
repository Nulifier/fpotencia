cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 0.8 LANGUAGES CXX)

include(FindDoxygen)
include(GNUInstallDirs)

option(BUILD_PROFILE "Build the profile for fPotencia" OFF)
if(NOT TARGET Eigen3::Eigen)
	find_package(Eigen3)
endif()

option(BUILD_UNSUPPORTED_SOLVERS
	"Build currently unsupported (i.e., test-failing) Load Flow solvers"
	OFF)

option(BUILD_TESTS
	"Build the tests for fPotencia"
	ON)

configure_file("config.h.in" "config.h")
configure_file("Doxyfile.in" "Doxyfile")

add_subdirectory(src)

if(BUILD_TESTS)
	find_package(Catch2 REQUIRED)
	enable_testing()
	add_subdirectory(test)
endif()

if(BUILD_PROFILE)
	add_subdirectory(profile)
endif()
if(DOXYGEN_FOUND)
	add_custom_target(doc
		COMMAND ${DOXYGEN_EXECUTABLE}
		DEPENDS Doxyfile.in)
endif()
