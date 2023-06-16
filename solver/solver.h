#pragma once

// #include <mpi.h>

// #include <iostream>
#include <memory>
// #include <array>
// #include <string>
#include <vector>
// #include <limits>
// #include <ctime>
// #include <cmath>

#include "definitions.h"
#include "apps/apps.h"



namespace solver
{

// Solver
void solver();

// Mimic RK solve with dg
const std::vector<std::shared_ptr<apps::WriteoutDiagnostic>> patch_time_advance();

} // namespace solver
