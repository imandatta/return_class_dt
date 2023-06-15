#pragma once

#include <mpi.h>

namespace mpi_boilerplate
{
bool init(int argc, char** argv);
void finalize();
} // namespace mpi_boilerplate
