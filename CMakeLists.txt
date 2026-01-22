cmake_minimum_required(VERSION 3.20)

set(CMAKE_Fortran_COMPILER gfortran)
project(REcoM LANGUAGES Fortran)
find_package(MPI REQUIRED)

# Get the list of Fortran sources in the folder
file(GLOB_RECURSE SOURCES "*.?90")

# Sets the fortran flag to allow files with no line length limit
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

# This is currently required because recom_sinking currently relies on this definition
add_compile_definitions(__recom)

# Create the target library for REcoM and add its dependencies
add_library(recom ${SOURCES})
target_link_libraries(recom PUBLIC MPI::MPI_Fortran)
