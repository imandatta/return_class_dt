#include "solver.h"

#include <mpi.h>

#include <limits>
#include <vector>
#include <array>
#include <algorithm>

namespace solver
{

struct DtInfo
{
    real dt;
    int rank;
};

void solver()
{
    // Do the time advance and get the write diagostic with min_dt out per patch
    const std::vector<std::shared_ptr<apps::WriteoutDiagnostic>>& wds =
        patch_time_advance();

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    const int num_apps = wds.size();
    std::vector<std::string> min_phys;
    std::vector<real> min_dt;
    std::vector<std::array<real, 3>> min_x;

    // for (auto& wd : wds)
    for (int app_idx = 0; app_idx < num_apps; app_idx++)
    {
        const std::shared_ptr<apps::WriteoutDiagnostic>& wd = wds[app_idx];

        // Now we need to coalesce the wd
        real my_dt = wd->getDt();

        // I want to allreduce dt but also carry the x and phyics
        // Can use MPI_MINLOC to reduct dt but also carry the rank of the wd which holds
        // it https://www.open-mpi.org/doc/v4.1/man3/MPI_Reduce.3.php
        DtInfo local_info;
        DtInfo global_info;

        // Seed for random number generation based on rank
        // std::srand(std::time(nullptr) + my_rank * 1E12); // use current time as seed
        // for random generator std::srand(std::time(nullptr)); // use current time as
        // seed for random generator std::srand(my_rank + 1); std::srand(my_rank); // use
        // current time as seed for random generator

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

        min_dt.push_back(comb_dt);
        min_x.push_back(comb_x);
        min_phys.push_back(app);

        if (local_info.rank == 0)
        {
            std::cout << "physics = " << app << " combined_dt = " << comb_dt
                      << " min rank = " << global_info.rank << " x = (" << comb_x[0]
                      << ", " << comb_x[1] << ", " << comb_x[2] << ")" << std::endl;
        }
    }

    std::vector<real>::iterator tmp = std::min_element(min_dt.begin(), min_dt.end());
    int min_dt_idx = std::distance(min_dt.begin(), tmp);
    if (my_rank == 0)
    {
        std::cout << "Global dt_min = " << min_dt[min_dt_idx] << " x = ("
                  << min_x[min_dt_idx][0] << ", " << min_x[min_dt_idx][1] << ", "
                  << min_x[min_dt_idx][2] << ")"
                  << " physics = " << min_phys[min_dt_idx] << std::endl;
    }
}

const std::vector<std::shared_ptr<apps::WriteoutDiagnostic>> patch_time_advance()
{
    // setup apps
    const int num_apps = 3;
    std::vector<std::unique_ptr<apps::AppBase>> apps;
    // for (int i = 0; i < num_apps; i++)
    // {
    //     apps.emplace_back(std::make_unique<App0>());
    //     apps.back()->setup();
    // }

    // app0
    apps.emplace_back(std::make_unique<apps::App0>());
    apps.back()->setup();

    // app1
    apps.emplace_back(std::make_unique<apps::App1>());
    apps.back()->setup();

    // app2
    apps.emplace_back(std::make_unique<apps::App2>());
    apps.back()->setup();

    // start with infinite dt
    std::vector<std::shared_ptr<apps::WriteoutDiagnostic>> sugg_wd;
    for (int app_idx = 0; app_idx < num_apps; app_idx++)
    {
        real sugg_dt = std::numeric_limits<real>::infinity();
        std::array<real, 3> xloc = {0, 0, 0};
        std::string slowest_physics = "nothing";
        // std::shared_ptr<apps::WriteoutDiagnostic> sugg_wd =
        // std::make_shared<apps::WriteoutDiagnostic>(sugg_dt, xloc, slowest_physics);
        sugg_wd.emplace_back(
            std::make_shared<apps::WriteoutDiagnostic>(sugg_dt, xloc, slowest_physics));
    }

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
        const std::array<real, 3> x = {patch_length * rank + node_dx * node_idx,
                                       static_cast<real>(rank) + 8.25,
                                       static_cast<real>(rank) * -6};
        // iterate through apps
        for (int app_idx = 0; app_idx < num_apps; app_idx++)
        {
            auto& app = apps[app_idx];
            std::shared_ptr<apps::WriteoutDiagnostic> wd = app->call(x);
            // *sugg_wd = apps::WriteoutDiagnostic::minDt(*sugg_wd, *wd);
            *sugg_wd[app_idx] = apps::WriteoutDiagnostic::minDt(*sugg_wd[app_idx], *wd);

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

} // namespace solver
