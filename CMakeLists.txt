cmake_minimum_required(VERSION 3.20)

project(REcoM LANGUAGES Fortran)
find_package(MPI REQUIRED)

# Get the list of Fortran sources in the folder
file(GLOB_RECURSE SOURCES "*.?90")

# Sets the fortran flag to allow files with no line length limit
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

# This is currently required because recom_sinking currently relies on this definition
add_compile_definitions(__recom)

# Create the target library for REcoM and add its dependencies
add_library(recom STATIC ${SOURCES})

if(DEFINED ENV{GITHUB_ACTIONS})
    target_compile_options(recom PRIVATE -O2 -g -fbacktrace -ffloat-store -finit-local-zero  -finline-functions -fimplicit-none  -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp)
endif()

target_link_libraries(recom PUBLIC MPI::MPI_Fortran)
