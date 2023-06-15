#include "apps.h"

#include <cmath>

namespace apps
{

// app0
const std::shared_ptr<WriteoutDiagnostic> App0::call(const std::array<real, 3>& x) const
{
    // Do the work of calculating physics

    // Calculate a timestep
    // For this test, generate a random number between 0 and 1
    // https://en.cppreference.com/w/cpp/numeric/random/rand
    // https://www.digitalocean.com/community/tutorials/random-number-generator-c-plus-plus
    // https://stackoverflow.com/questions/686353/random-float-number-generation
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // const real lo = 0;
    // const real hi = 1.0;
    // real dt = lo + static_cast <real> (rand()) /( static_cast <real>
    // (RAND_MAX/(hi-lo))); std::cout << ")rank = " << rank <<  " dt = " << dt
    // 		  <<  " x = (" << x[0] << ", " << x[1] << ", " << x[2]
    // 		  << ") physics = " << physics << std::endl;

    const real dt = 1.0;
    // return a timestep
    _wd->setValues(dt, x, _physics);

    return _wd;
}

// app1
const std::shared_ptr<WriteoutDiagnostic> App1::call(const std::array<real, 3>& x) const
{
    // Do the work of calculating physics

    // const real dt = x[0] * x[0] - 6;
    const real dt = std::sin(2 * 3.141592654 * x[0] / 1);
    // return a timestep
    _wd->setValues(dt, x, _physics);

    return _wd;
}

// app2
const std::shared_ptr<WriteoutDiagnostic> App2::call(const std::array<real, 3>& x) const
{
    {
        // Do the work of calculating physics

        const real dt = 0.5 * std::sqrt(x[0]);
        // return a timestep
        _wd->setValues(dt, x, _physics);

        return _wd;
    }
}

} // namespace apps
