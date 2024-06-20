# To generate build files for common and all example cases:
#   mkdir build
#   cd build
#   cmake ..
#   make
# Regenerate case files after changes:
#   cmake .
#   make case
# Show each build command line as it is launched
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..

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
    set(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    # same compile options as Makefile.in
    add_compile_options(-fallow-argument-mismatch -freal-4-real-8 -O3)
    add_link_options(-ldl -lstdc++ -Wl,--copy-dt-needed-entries)
else ()
    message(FATAL_ERROR "PETSc Fortran compiler not found. Please ensure PETSc was configured with a Fortran compiler.")
endif (FORTRAN_COMPILER)

find_package(PkgConfig REQUIRED)

if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
else()
    message(FATAL_ERROR "PkgConfig was not found.")
endif()

if (NOT PETSC_FOUND)
    message(FATAL_ERROR "PETSc was not found through PkgConfig.")
else()
    # message("PETSc was found")
    # cmake_print_variables(PETSC_INCLUDE_DIRS)
    # cmake_print_variables(PETSC_LINK_LIBRARIES)
endif()

if (LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
endif(LAPACK_FOUND)

# look for FFTWROOT or FFTW_ROOT
if (DEFINED ENV{FFTWROOT})
    set(FFTW_LIB "-L$ENV{FFTWROOT}/lib -lfftw3")
    # cmake_print_variables(FFTW_LIB)
elseif(DEFINED ENV{FFTWROOT})
    set(FFTW_LIB "-L$ENV{FFTW_ROOT}/lib -lfftw3")
    # cmake_print_variables(FFTW_LIB)
else ()
    message(FATAL_ERROR "FFTWROOT and FFTW_ROOT environment variables not found. Please module load fftw.")
endif()

# look for NETCDF_FORTRANROOT or NETCDF_FORTRAN_ROOT 
if (DEFINED ENV{NETCDF_FORTRANROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRANROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRANROOT}/lib -lnetcdff")
    # cmake_print_variables(NETCDF_INCLUDE)
    # cmake_print_variables(NETCDF_LIB)
elseif (DEFINED ENV{NETCDF_FORTRAN_ROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRAN_ROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRAN_ROOT}/lib -lnetcdff")
else ()
    message(FATAL_ERROR "NETCDF_FORTRANROOT or NETCDF_FORTRAN_ROOT environment variables not found. Please module load netcdf-fortran.")
endif()

# spherepack library
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/packages/spherepack3.2/lib")
    set(SPHPK_LIB "-L${CMAKE_CURRENT_SOURCE_DIR}/packages/spherepack3.2/lib -lspherepack")
else()
    message(FATAL_ERROR "Spherepack is not installed in ${CMAKE_CURRENT_SOURCE_DIR}/packages/spherepack3.2")
endif()

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

# cmake_print_variables(dirlist)

foreach (casedir ${dirlist})
    # cmake_print_variables(casedir)

    # make binary dir inside /build for case executables
    set(case_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${casedir}")
    file(MAKE_DIRECTORY ${case_binary_dir})
    make_directory("${case_binary_dir}/D")
    file(GLOB inputs "${CMAKE_CURRENT_SOURCE_DIR}/${examples_dir}/${casedir}/Input/*")
    # cmake_print_variables(inputs)
    file(COPY ${inputs} DESTINATION "${case_binary_dir}/Input")
    # cmake_print_variables(case_binary_dir)

    # find all F90 files in the examples/case directory
    file(GLOB case_F90s "${CMAKE_CURRENT_SOURCE_DIR}/${examples_dir}/${casedir}/*.F90")
    # cmake_print_variables(case_F90s)

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
