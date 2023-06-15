cmake_minimum_required(VERSION 3.24.1)
project(mytest)


set(CMAKE_CXX_STANDARD 17)


add_subdirectory(apps)
set(SOURCES apps/apps.cc)
add_executable(mytest main.cc ${SOURCES})





# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(mytest ${MPI_LIBRARIES} )



