cmake_minimum_required(VERSION 3.24.1)
project(mytest)


set(CMAKE_CXX_STANDARD 17)


add_executable(mytest main.cc)



# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(mytest ${MPI_LIBRARIES})
