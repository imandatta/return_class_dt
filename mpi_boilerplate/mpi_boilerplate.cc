#include "mpi_boilerplate.h"

#include <iostream>


namespace mpi_boilerplate
{

	
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

} // namespace mpi_boilerplate
