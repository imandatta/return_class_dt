// clang++ -std=c++17 main.cc

#include <mpi.h>

#include <iostream>
#include <memory>
#include <array>
#include <string>
#include <vector>
#include <limits>
#include <ctime>
#include <cmath>

typedef double real;

class WriteoutDiagnostic
{
private:
    real _dt;
    std::array<real, 3> _x;
    std::string _physics;

public:
    // Constructor
    WriteoutDiagnostic()
    {
    }
    // Constructor
    WriteoutDiagnostic(real dt, const std::array<real, 3>& x, const std::string& physics)
        : _dt(dt), _x(x), _physics(physics)
    {
    }

    // Minimum dt
    static const WriteoutDiagnostic& minDt(const WriteoutDiagnostic& a,
                                           const WriteoutDiagnostic& b)
    {
        return a._dt < b._dt ? a : b;
    }

    // Setters
    void setDt(real dt)
    {
        _dt = dt;
    }
    void setX(const std::array<real, 3>& x)
    {
        _x = x;
    }
    void setPhysics(const std::string& physics)
    {
        _physics = physics;
    }

    void setValues(real dt, const std::array<real, 3>& x, const std::string& physics)
    {
        setDt(dt);
        setX(x);
        setPhysics(physics);
    }

    // Getters
    real getDt() const
    {
        return _dt;
    }

    const std::array<real, 3>& getX() const
    {
        return _x;
    }

    const std::string& getPhysics() const
    {
        return _physics;
    }
};

class AppBase
{
protected:
    std::shared_ptr<WriteoutDiagnostic> _wd;
    std::string _physics;

public:
    // Constructor
    AppBase()
    {
        _wd = std::make_shared<WriteoutDiagnostic>();
    }

    virtual void setup()
    {
    }

    virtual const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const
    {
        std::cerr << "Unimplemented call()\n";
        exit(-1);
    }
};

class App0 : public AppBase
{
private:
public:
    // Constructor
    App0()
    {
        _physics = "app0";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override

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
};

class App1 : public AppBase
{
private:
public:
    // Constructor
    App1()
    {
        _physics = "app1";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override

    {
        // Do the work of calculating physics

        // const real dt = x[0] * x[0] - 6;
        const real dt = std::sin(2 * 3.141592654 * x[0] / 1);
        // return a timestep
        _wd->setValues(dt, x, _physics);

        return _wd;
    }
};

class App2 : public AppBase
{
private:
public:
    // Constructor
    App2()
    {
        _physics = "app2";
    }

    // Setup
    void setup() override
    {
        AppBase::setup();
    }

    const std::shared_ptr<WriteoutDiagnostic>
    call(const std::array<real, 3>& x) const override

    {
        // Do the work of calculating physics

        const real dt = 0.5 * std::sqrt(x[0]);
        // return a timestep
        _wd->setValues(dt, x, _physics);

        return _wd;
    }
};

// Mimic RK solve with dg
const std::shared_ptr<WriteoutDiagnostic> patch_time_advance()
{
    // setup apps
    const int num_apps = 3;
    std::vector<std::unique_ptr<AppBase>> apps;
    // for (int i = 0; i < num_apps; i++)
    // {
    //     apps.emplace_back(std::make_unique<App0>());
    //     apps.back()->setup();
    // }

    // app0
    apps.emplace_back(std::make_unique<App0>());
    apps.back()->setup();

    // app1
    apps.emplace_back(std::make_unique<App1>());
    apps.back()->setup();

    // app2
    apps.emplace_back(std::make_unique<App2>());
    apps.back()->setup();

    // start with infinite dt
    real sugg_dt = std::numeric_limits<real>::infinity();
    std::array<real, 3> xloc = {0, 0, 0};
    std::string slowest_physics = "nothing";
    std::shared_ptr<WriteoutDiagnostic> sugg_wd =
        std::make_shared<WriteoutDiagnostic>(sugg_dt, xloc, slowest_physics);

    // setup a vector of apps

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // assume in x = [0, 1], y = 0, z = 0
    const real domain_min = 0.0;
    const real domain_max = 1.0;
    const real domain_length = domain_max - domain_min;
    const real patch_length = domain_length / static_cast<real>(size);

    // this is like nodes in a patch
    const int num_nodes = 10;
    const real node_dx = patch_length / static_cast<real>(num_nodes);
    // iterate through nodes
    for (int node_idx = 0; node_idx < num_nodes; node_idx++)
    {
        // generate some arbitrary position for this test
        // const std::array<real, 3> x = {node_idx * 1.0 + 10, node_idx * 10.0 - 30,
        // node_idx / 2.0 + 6};

        // assume in x = [0, 1], y = 0, z = 0
        const std::array<real, 3> x = {
            patch_length * rank + node_dx * node_idx, 0.0, 0.0};
        // iterate through apps
        for (int app_idx = 0; app_idx < num_apps; app_idx++)
        {
            auto& app = apps[app_idx];
            // std::vector<apps>; avoids having a physics string
            // grab the minimum timestep and x for each app
            // maybe use mpi_double
            const real dt = 3 * (app_idx - 1) * 0.1 + 0.5 + num_apps - 4 * app_idx + 1 -
                            3 * rank - rank * rank;

            // std::string physics = "app" + std::to_string(app_idx);

            // App app;
            // app.setup();
            std::shared_ptr<WriteoutDiagnostic> wd = app->call(x);
            *sugg_wd = WriteoutDiagnostic::minDt(*sugg_wd, *wd);

            // another app that does the same thing but prints
        }

        // std::cout << "dt = " << wd->getDt() << std::endl;
        // std::cout << "(x, y, z) = (" << wd->getX()[0] << ", " << wd->getX()[1] << ", "
        // << wd->getX()[2] << ")\n"; std::cout << "physics = " << wd->getPhysics() <<
        // std::endl;
    }

    // std::cout << "\nApp with minimum timestep is:\n";
    // std::cout << "dt = " << sugg_wd->getDt() << std::endl;
    // std::cout << "(x, y, z) = (" << sugg_wd->getX()[0] << ", " << sugg_wd->getX()[1] <<
    // ", " << sugg_wd->getX()[2] << ")\n"; std::cout << "physics = " <<
    // sugg_wd->getPhysics() << std::endl;

    return sugg_wd;
}

// Solver
void solver()
{
    // Do the time advance and get the write diagostic with min_dt out per patch
    const std::shared_ptr<WriteoutDiagnostic> wd = patch_time_advance();

    // Now we need to coalesce the wd
    real my_dt = wd->getDt();

    // I want to allreduce dt but also carry the x and phyics
    // Can use MPI_MINLOC to reduct dt but also carry the rank of the wd which holds it
    // https://www.open-mpi.org/doc/v4.1/man3/MPI_Reduce.3.php
    struct DtInfo
    {
        real dt;
        int rank;
    };

    DtInfo local_info;
    DtInfo global_info;

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Seed for random number generation based on rank
    // std::srand(std::time(nullptr) + my_rank * 1E12); // use current time as seed for
    // random generator std::srand(std::time(nullptr)); // use current time as seed for
    // random generator std::srand(my_rank + 1); std::srand(my_rank); // use current time
    // as seed for random generator

    local_info.dt = my_dt;
    local_info.rank = my_rank;
    // MPI_Allreduce(&my_dt, &comb_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(
        &local_info, &global_info, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    real comb_dt = global_info.dt;
    std::array<real, 3> comb_x = wd->getX();
    MPI_Bcast(comb_x.data(), 3, MPI_DOUBLE, global_info.rank, MPI_COMM_WORLD);

    // https://stackoverflow.com/questions/46331763/what-is-the-most-elegant-way-of-broadcasting-stdstring-in-mpi
    std::string app = wd->getPhysics();
    int app_string_size = app.size();
    // Variable app receives a value on mpiid = global_info.rank
    MPI_Bcast(&app_string_size, 1, MPI_INT, global_info.rank, MPI_COMM_WORLD);
    if (my_rank != global_info.rank)
    {
        app.resize(app_string_size);
    }

    MPI_Bcast(const_cast<char*>(app.data()),
              app_string_size,
              MPI_CHAR,
              global_info.rank,
              MPI_COMM_WORLD);

    if (local_info.rank == 0)
    {
        std::cout << "combined_dt = " << comb_dt << " min rank = " << global_info.rank
                  << " x = (" << comb_x[0] << ", " << comb_x[1] << ", " << comb_x[2]
                  << ") physics = " << app << std::endl;
    }
}

bool init(int argc, char** argv)
{
    // initialize MPI
    int mpi_thread_support_provided;
    int result =
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_thread_support_provided);
    if (result != MPI_SUCCESS)
    {
        std::cerr << "couldn't initialize MPI correctly\n";
        return result;
    }
    if ((mpi_thread_support_provided != MPI_THREAD_FUNNELED) || (result != MPI_SUCCESS))
    {
        std::cerr << "Error: MPI could not be initialized with required support for "
                     "MPI_THREAD_FUNNELED. "
                  << std::endl;
        return false;
    }

    return true;
}

void finalize()
{
    MPI_Finalize();
}

int main(int argc, char** argv)
{
    init(argc, argv);
    solver();
    finalize();

    return 0;
}
