set(fpotencia_SOURCES
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
    ExternalGrid.cpp)

set(fpotencia_HEADERS
    Solver.h
    UnfittingSolverException.h
    Circuit.h
    TransformerConstructors.h
    Shunt.h
    Transformer.h
    Line.h
    enumaratons.h
    Load.h
    Bus.h
    fpotencia_libs.h
    Solution.h
    fpotencia.h
    ExternalGrid.h
    Solver_NRpolar.h
    LineConstructors.h
    Solver_NRcurrent.h
    CommonFunctions.h
    Generator.h)

set(fpotencia_UNSUPPORTED_SOLVERS_HEADERS
    Solver_Jacobi.h
    Solver_ZBusGS.h
    Solver_NRrect.h)

set(fpotencia_UNSUPPORTED_SOLVERS_SOURCES
    Solver_Jacobi.cpp
    Solver_ZBusGS.cpp
    Solver_NRrect.cpp)

if(BUILD_UNSUPPORTED_SOLVERS)
    list(APPEND ${fpotencia_SOURCES} ${fpotencia_UNSUPPORTED_SOLVERS_SOURCES})
endif()

add_library(fpotencia ${fpotencia_SOURCES})

target_include_directories(fpotencia
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(fpotencia
    PUBLIC
        #Boost::graph
        Eigen3::Eigen
)

set_property(TARGET fpotencia PROPERTY CXX_STANDARD 17)
set_property(TARGET fpotencia PROPERTY SOVERSION ${VERSION})

install(TARGETS fpotencia
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${fpotencia_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fpotencia)
