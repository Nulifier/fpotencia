configure_file("config.h.in" "config.h")

add_library(fpotencia
	Solver.cpp
	UnfittingSolverException.cpp
	LineConstructors.cpp
	Solution.cpp
	Circuit.cpp
	Generator.cpp
	TransformerConstructors.cpp
	Solver_NRcurrent.cpp
	Shunt.cpp
	Load.cpp
	Line.cpp
	Bus.cpp
	Solver_NRpolar.cpp
	CommonFunctions.cpp
	Transformer.cpp
	ExternalGrid.cpp
)

# Set the C++ standard to C++17
target_compile_features(fpotencia PUBLIC cxx_std_17)
set_target_properties(fpotencia PROPERTIES CXX_EXTENSIONS OFF)

target_include_directories(fpotencia
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}
)

if(USE_OPENMP AND NOT OpenMP_CXX_FOUND)
	message(WARNING "Requested OpenMP support but the library was not installed")
endif()

target_link_libraries(fpotencia
	PUBLIC
		Eigen3::Eigen
		# Include OpenMP if requested and found
		$<$<AND:$<BOOL:OpenMP_CXX_FOUND>,$<BOOL:USE_OPENMP>>:OpenMP::OpenMP_CXX>
)

install(TARGETS fpotencia
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${fpotencia_HEADERS}
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fpotencia)
