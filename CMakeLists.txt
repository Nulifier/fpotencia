cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 0.8 LANGUAGES CXX)

include(FindDoxygen)
include(GNUInstallDirs)

option(BUILD_PROFILE "Build the profile for fPotencia" OFF)

if(NOT TARGET Eigen3::Eigen)
	find_package(Eigen3)
endif()

add_subdirectory(src)

include(CTest)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
	find_package(Catch2 3 REQUIRED)
	add_subdirectory(test)
endif()

if(BUILD_PROFILE)
	add_subdirectory(profile)
endif()

if(DOXYGEN_FOUND)
	configure_file("Doxyfile.in" "Doxyfile")
	add_custom_target(doc
		COMMAND ${DOXYGEN_EXECUTABLE}
		DEPENDS Doxyfile.in)
endif()
