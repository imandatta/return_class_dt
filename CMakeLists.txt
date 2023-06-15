cmake_minimum_required(VERSION 3.24.1)
project(mytest)


set(CMAKE_CXX_STANDARD 17)


add_subdirectory(apps)
set(APP_SOURCES apps/apps.cc)
set(SOLVER_SOURCES solver/solver.cc)
set(MPI_BOILERPLATE_SOURCES mpi_boilerplate/mpi_boilerplate.cc)
add_executable(mytest main.cc
					  ${APP_SOURCES}
					  ${SOLVER_SOURCES}
					  ${MPI_BOILERPLATE_SOURCES})
target_include_directories(mytest PRIVATE ${PROJECT_SOURCE_DIR})


# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(mytest ${MPI_LIBRARIES} )



