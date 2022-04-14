# fPotencia --- A Library for Power System Load Flow Analysis


## About

fPotencia is a C++ library that helps to define a power grid and perform power
system load flow analysis on it. It currently offers two Newton-Raphson solver
algorithms:

  - Newton-Raphson using polar coordinates (including µ acceleration
    parameter)
  - Newton-Raphson with current injection


## Building and Installing

### Prerequisites

  - CMake >= 3.0.0
  - Eigen >= 3.2.8

Optionally, if you want testing:

  - GTest

Optionally, if the API documentation is desired:

  - Doxygen (for API documentation)

Of course, a working C++ compiler and the C++ Standard Template Library is
also required. The C++ compiler must understand C++11.

### Building

To create a release build of a shared library, e.g., for packaging, issue:

    mkdir Build && cd Build
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON \        # If OFF, a static lib will be built
        ..
    make
    make test                           # Optionally
    make doc                            # Optionally
    make install

If you would like to contribute to fPotencia, `-DCMAKE_BUILD_TYPE=Debug` is
obviously the preferred choice.


## Getting Started

The best way to get an immediate overview of how fPotencia works is to check
out the test grid definitions in `test/SolverTest.cpp`. They show how a grid
is defined.

Of course, for a deeper dive, the API documentation is preferrable. However,
if you just want to use fPotencia to perform a load flow analysis of a grid
you supply, the example grids in `test/SolverTest.cpp` are the only
instruction you're going to need.


## Roadmap

### 0.9

  - API overhaul
  - Performance Improvements

### 1.0

  - More test cases
  - Reliable PV <-> PQ switching
