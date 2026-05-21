cmake_minimum_required(VERSION 3.20)

project(REcoM LANGUAGES Fortran)
find_package(MPI REQUIRED)

# Get the list of Fortran sources in the folder
file(GLOB_RECURSE SOURCES "*.?90")

# This is currently required because recom_sinking relies on this definition
add_compile_definitions(__recom)

# Create the target library for REcoM and add its dependencies
add_library(recom STATIC ${SOURCES})

# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL  IntelLLVM )

    # Base compiler flags
    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -r8 -i4 -fp-model precise -no-prec-div -fimf-use-svml -init=zero -no-wrap-margin -fpe0 -fpp)
    target_compile_options(${PROJECT_NAME} PRIVATE -march=core-avx2 -fPIC -qopt-malloc-options=2 -qopt-prefetch=5 -unroll-aggressive)

    # compiler flags not supported by IntelLLVM
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
        target_compile_options(${PROJECT_NAME} PRIVATE -no-prec-sqrt -ip )
    endif()

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )

    # Sets the fortran flag to allow files with no line length limit
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
    target_compile_options(recom PRIVATE -ffree-line-length-none)

endif()

if(DEFINED ENV{GITHUB_ACTIONS})

    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
        target_compile_options(recom PRIVATE -O2 -g -fbacktrace -ffloat-store -finit-local-zero  -finline-functions -fimplicit-none -fdefault-real-8 -fdefault-double-8 -cpp)
    endif()

endif()

target_link_libraries(recom PUBLIC MPI::MPI_Fortran)
