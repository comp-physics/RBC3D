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
# See compile commands in build directory
#   cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..

include(CMakePrintHelpers)
cmake_minimum_required(VERSION 3.20)

project(RBC3D LANGUAGES Fortran C CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(USE_MKL OFF)

if (DEFINED ENV{MKLROOT})
    set(USE_MKL ON)
    message(STATUS "MKLROOT: $ENV{MKLROOT}")   
    # use mkl implementation of lapack
    set(BLA_VENDOR Intel10_64lp)
endif()

if (NOT DEFINED ENV{PETSC_DIR})
    message(FATAL_ERROR "PETSc environment variables are not set. Please export PETSC_DIR.")
endif()

if (NOT DEFINED ENV{PETSC_ARCH})
    # PETSC_ARCH is not part of PETSc module
    message(WARNING "PETSC_ARCH is not set. Please export PETSC_ARCH if using manually installed PETSc package.")
    set(PETSC $ENV{PETSC_DIR})
else ()
    # full path if PETSC_DIR and PETSC_ARCH env variables are set
    set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
    if ("$ENV{PETSC_ARCH}" STREQUAL "arch-darwin-c-opt")
        set(MAC on)
    elseif("$ENV{PETSC_ARCH}" STREQUAL "arch-linux-c-opt")
        set(LINUX on)
    else ()
        message(FATAL_ERROR "PETSC_ARCH $ENV{PETSC_ARCH} not recognized. Please add support in CMakeLists.txt")
    endif()
    
endif()

set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

cmake_print_variables(PETSC)

# use the same fortran compiler petsc was configured with
execute_process (COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
cmake_print_variables(FORTRAN_COMPILER)

if (FORTRAN_COMPILER)
    set(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    # same compile options as Makefile.in
    add_compile_options(-fallow-argument-mismatch -freal-4-real-8 -O3)
    add_link_options(-ldl -lstdc++)
    if (LINUX)
        # ice cluster install needs this extra link option
        add_link_options(-Wl,--copy-dt-needed-entries)
    endif()
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
    if (NOT DEFINED ENV{PETSC_ARCH})
        # not having PETSC_ARCH set makes include dirs null for some reason
        set(PETSC_INCLUDE_DIRS "$ENV{PETSC_DIR}/include")
    endif()
endif()

if (LAPACK_FOUND)
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
elseif(NOT USE_MKL)
    # LAPACK_LIBRARIES doesn't need to be set if MKL is module loaded
    set(LAPACK_LIBRARIES "-L${CMAKE_CURRENT_SOURCE_DIR}/packages/lapack-3.11 -llapack -lrefblas")
    cmake_print_variables(LAPACK_LIBRARIES)
endif(LAPACK_FOUND)

# look for FFTWROOT or FFTW_ROOT
if (DEFINED ENV{FFTWROOT})
    set(FFTW_INCLUDE "$ENV{FFTWROOT}/include")
    set(FFTW_LIB "-L$ENV{FFTWROOT}/lib -lfftw3")
elseif(DEFINED ENV{FFTW_ROOT})
    set(FFTW_INCLUDE "$ENV{FFTW_ROOT}/include")
    set(FFTW_LIB "-L$ENV{FFTW_ROOT}/lib -lfftw3")
else ()
    message(WARNING "FFTWROOT and FFTW_ROOT environment variables from module load fftw not found. Assuming you installed FFTW in packages.")
    set(FFTW_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/packages/fftw-3.3.10/build/include")
    set(FFTW_LIB "-L${CMAKE_CURRENT_SOURCE_DIR}/packages/fftw-3.3.10/build/lib -lfftw3")
endif()

# look for NETCDF_FORTRANROOT or NETCDF_FORTRAN_ROOT 
if (DEFINED ENV{NETCDF_FORTRANROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRANROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRANROOT}/lib -lnetcdff")
elseif (DEFINED ENV{NETCDF_FORTRAN_ROOT})
    set(NETCDF_INCLUDE "$ENV{NETCDF_FORTRAN_ROOT}/include")
    set(NETCDF_LIB "-L$ENV{NETCDF_FORTRAN_ROOT}/lib -lnetcdff")
else ()
    message(WARNING "NETCDF_FORTRANROOT and NETCDF_FORTRAN_ROOT environment variables from module load fftw not found. Assuming you installed NETCDF-FORTRAN in packages.")
    set(NETCDF_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/packages/NETCDF_INST/include")
    # -rpath linker option to avoid having to add netcdf library to LD_LIBRARY_PATH
    set(NETCDF_LIB "-Wl,-rpath -Wl,${CMAKE_CURRENT_SOURCE_DIR}/packages/NETCDF_INST/lib -L${CMAKE_CURRENT_SOURCE_DIR}/packages/NETCDF_INST/lib -lnetcdff")
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

foreach (casedir ${dirlist})

    # make binary dir inside /build for case executables
    set(case_binary_dir "${CMAKE_CURRENT_BINARY_DIR}/${casedir}")
    file(MAKE_DIRECTORY ${case_binary_dir})
    # copy input directory to case binary dir and create D directory
    make_directory("${case_binary_dir}/D")
    file(GLOB inputs "${CMAKE_CURRENT_SOURCE_DIR}/${examples_dir}/${casedir}/Input/*")
    file(COPY ${inputs} DESTINATION "${case_binary_dir}/Input")

    # find all F90 files in the examples/case directory
    file(GLOB case_F90s "${CMAKE_CURRENT_SOURCE_DIR}/${examples_dir}/${casedir}/*.F90")

    # create a custom target with the name of the current case directory
    add_custom_target(${casedir})

    # loop through each F90 file
    foreach(F90_filepath ${case_F90s})
        # get the file name without extension
        get_filename_component(file_target ${F90_filepath} NAME_WE)
        add_executable(${file_target} ${F90_filepath})

        # add each target to its own case build directory
        set_target_properties(${file_target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${case_binary_dir})

        # set include directories and link libraries
        # petsc might not be necessary for all files (i.e. initcond) but it doesn't hurt to add it
        target_include_directories(${file_target} PRIVATE ${PETSC_INCLUDE_DIRS})
        target_link_libraries(${file_target} PRIVATE ${LAPACK_LIBRARIES} common ${PETSC_LINK_LIBRARIES})

        # add dependency to the custom target
        add_dependencies(${casedir} ${file_target})
    endforeach()
endforeach()
