# To generate build files:
#   mkdir build
#   cd build
#   cmake ..
#   make
# Regenerate after changes:
#   cmake .
# Show each build command line as it is launched
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..
# Thank you Boost, Sundials

include(CMakePrintHelpers)
cmake_minimum_required(VERSION 3.20)

project(RBC3D LANGUAGES Fortran C CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# use mkl implementation of lapack
set(BLA_VENDOR Intel10_64lp)

if (NOT DEFINED ENV{PETSC_DIR} OR NOT DEFINED ENV{PETSC_ARCH})
    message(FATAL_ERROR "PETSc environment variables are not set. Please export PETSC_DIR and PETSC_ARCH.")
endif()

# PETSC_DIR and PETSC_ARCH env variable must be set
set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

cmake_print_variables(PETSC)

# use the same fortran compiler petsc was configured with
execute_process (COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
cmake_print_variables(FORTRAN_COMPILER)

if (FORTRAN_COMPILER)
    message(${FORTRAN_COMPILER})
    set(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    add_compile_options(-fallow-argument-mismatch -freal-4-real-8 -O3)
    add_link_options(-ldl -lstdc++ -Wl,--copy-dt-needed-entries)

    # set(CMAKE_Fortran_LINK_EXECUTABLE "${CMAKE_Fortran_LINK_EXECUTABLE} -ldl")
else ()
    message(FATAL_ERROR "PETSc Fortran compiler not found. Please export PETSC_DIR and PETSC_ARCH.")
endif (FORTRAN_COMPILER)

find_package(PkgConfig REQUIRED)
# find_package(LAPACK REQUIRED)

if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
else()
    message(FATAL_ERROR "PkgConfig was not found.")
endif()

if (NOT PETSC_FOUND)
    message(FATAL_ERROR "PETSc was not found through PkgConfig.")
else()
    message("PETSc was not found")
    cmake_print_variables(PETSC_INCLUDE_DIRS)
    cmake_print_variables(PETSC_LINK_LIBRARIES)
endif()

if (LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
endif(LAPACK_FOUND)

# look for FFTWROOT or FFTW_ROOT
if (DEFINED ENV{FFTWROOT})
    # set(FFTW_INCLUDE "${FFTW_INCLUDE};$ENV{FFTWROOT}/include")
    set(FFTW_LIB "-L$ENV{FFTWROOT}/lib -lfftw3")
    cmake_print_variables(FFTW_LIB)
elseif(DEFINED ENV{FFTWROOT})
    set(FFTW_LIB "-L$ENV{FFTW_ROOT}/lib -lfftw3")
    cmake_print_variables(FFTW_LIB)
else ()
    message(FATAL_ERROR "FFTWROOT and FFTW_ROOT environment variables not found. Please module load fftw.")
endif()

# look for NETCDF_FORTRANROOT or NETCDF_FORTRAN_ROOT 
if (DEFINED ENV{NETCDF_FORTRANROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRANROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRANROOT}/lib -lnetcdff")
    cmake_print_variables(NETCDF_INCLUDE)
    cmake_print_variables(NETCDF_LIB)
elseif (DEFINED ENV{NETCDF_FORTRAN_ROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRAN_ROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRAN_ROOT}/lib -lnetcdff")
else ()
    message(FATAL_ERROR "NETCDF_FORTRANROOT or NETCDF_FORTRAN_ROOT environment variables not found. Please module load netcdf-fortran.")
endif()

# spherepack library
set(SPHPK_LIB "-L${CMAKE_CURRENT_SOURCE_DIR}/packages/spherepack3.2/lib -lspherepack")

# build common
add_subdirectory(common)
install(TARGETS common LIBRARY DESTINATION "bin")

include_directories("${CMAKE_CURRENT_BINARY_DIR}/common")

set(examples_dir "examples")

# glob for the directories in /examples
file(GLOB examples_list RELATIVE ${CMAKE_SOURCE_DIR}/${examples_dir} ${CMAKE_SOURCE_DIR}/${examples_dir}/*)
set(dirlist "")
foreach (example ${examples_list})
    if (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/${examples_dir}/${example})
        list(APPEND dirlist ${example})
    endif()
endforeach()

cmake_print_variables(dirlist)

set(ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

foreach (casedir ${dirlist})
    cmake_print_variables(${casedir})

    # make binary dir inside /build for case files
    set(case_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${casedir}")
    file(MAKE_DIRECTORY ${case_binary_dir})
    cmake_print_variables(case_binary_dir)

    # find all F90 files in the examples/case directory
    file(GLOB case_F90s "${CMAKE_CURRENT_SOURCE_DIR}/${examples_dir}/${casedir}/*.F90")
    cmake_print_variables(case_F90s)

    # create a custom target with the name of the current case directory
    add_custom_target(${casedir})

    # loop through each F90 file
    foreach(F90_filepath ${case_F90s})
        # get the file name without extension
        get_filename_component(file_target ${F90_filepath} NAME_WE)
        add_executable(${file_target} ${F90_filepath})

        # add target to its own case build directory
        set_target_properties(${file_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${case_binary_dir})

        # set include directories and link libraries
        # petsc might not be necessary for all files (i.e. initcond) but it doesn't hurt to add it
        target_include_directories(${file_target} PRIVATE ${PETSC_INCLUDE_DIRS})
        target_link_libraries(${file_target} PRIVATE common ${PETSC_LINK_LIBRARIES})

        # add dependency to the custom target
        add_dependencies(${casedir} ${file_target})
    endforeach()
endforeach()

# # make binary dir inside /build for case files
# set(CASE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/case")
# file(MAKE_DIRECTORY ${CASE_BINARY_DIR})

# add_subdirectory(examples/case case)
# add_subdirectory(examples/case_sickles case_sickles)
# # Find all Fortran files in the examples/case directory
# file(GLOB case_F90s "${CMAKE_CURRENT_SOURCE_DIR}/examples/case/*.F90")
# cmake_print_variables(case_F90s)

# # Create a custom target with the name of the current directory
# add_custom_target(case)

# # Loop through each Fortran file
# foreach(F90_filepath ${case_F90s})
#     # Get the file name without extension
#     get_filename_component(file_target ${F90_filepath} NAME_WE)
#     add_executable(${file_target} ${F90_filepath})
#     cmake_print_variables(file_target)

#     set_target_properties(${file_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CASE_BINARY_DIR})

#     # Set include directories and link libraries
#     target_include_directories(${file_target} PRIVATE ${PETSC_INCLUDE_DIRS})
#     target_link_libraries(${file_target} PRIVATE common ${PETSC_LINK_LIBRARIES})

#     # Add dependency to the custom target
#     add_dependencies(case ${file_target})
# endforeach()
