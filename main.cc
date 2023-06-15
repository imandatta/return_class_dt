// clang++ -std=c++17 main.cc

// #include <mpi.h>

// #include <iostream>
// #include <memory>
// #include <array>
// #include <string>
// #include <vector>
// #include <limits>
// #include <ctime>
// #include <cmath>

#include "apps/apps.h"
#include "solver/solver.h"
#include "mpi_boilerplate/mpi_boilerplate.h"

// typedef double real;


int main(int argc, char** argv)
{
	mpi_boilerplate::init(argc, argv);
	solver::solver();
	mpi_boilerplate::finalize();

    return 0;
}
