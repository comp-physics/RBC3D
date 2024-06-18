# To build /common (generate libcommon in build dir):
#   mkdir build
#   cd build
#   cmake -DBUILD_COMMON+ON ..
#   make
# Regenerate after common changes:
#   cmake .
# Show each build command line as it is launched
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON ..
# Thank you Boost

include(CMakePrintHelpers)
cmake_minimum_required(VERSION 3.20)

project(RBC3D LANGUAGES Fortran C CXX)

option(BUILD_COMMON  "Build common"  OFF)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# use mkl implementation of lapack
set(BLA_VENDOR Intel10_64lp)

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
endif()

if (PETSC_FOUND)
    message("PETSc found")
    cmake_print_variables(PETSC_INCLUDE_DIRS)
    cmake_print_variables(PETSC_LINK_LIBRARIES)
  # set(PETSC_INCLUDE "$ENV{PETSC_DIR}/include")
  # cmake_print_variables(PETSC_INCLUDE)

  # set(PETSC_ARCH_INCLUDE "${PETSC}/include")
  # cmake_print_variables(PETSC_ARCH_INCLUDE)

#   list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})
#   # include_directories(${PETSC_INCLUDE_DIRS})
#   set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")   
#   list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})   
#   set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})   
#   set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
else ()
    message(FATAL_ERROR "PETSc not found. Please export PETSC_DIR and PETSC_ARCH.")
endif(PETSC_FOUND)

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
if (BUILD_COMMON)
    file(GLOB common_F90s "${CMAKE_CURRENT_SOURCE_DIR}/common/*.F90")
    add_library(common STATIC "${common_F90s}")
    # target_link_libraries(common LAPACK::LAPACK PkgConfig::PETSC)
    target_include_directories(common PRIVATE ${PETSC_INCLUDE_DIRS} ${NETCDF_INCLUDE})
    target_link_libraries(common PRIVATE ${PETSC_LINK_LIBRARIES} ${SPHPK_LIB} ${FFTW_LIB} ${NETCDF_LIB})
    install(TARGETS common RUNTIME DESTINATION bin)
endif()

# set(TARGET_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/examples/case/*.F90 ${CMAKE_CURRENT_SOURCE_DIR}/examples/case/*.F90)
# cmake_print_variables(TARGET_SOURCES)
file(GLOB inputs "${CMAKE_CURRENT_SOURCE_DIR}/examples/case/Input/*")
cmake_print_variables(inputs)
file(COPY ${inputs}
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/Input")
cmake_print_variables(CMAKE_CURRENT_BINARY_DIR)
make_directory("${CMAKE_CURRENT_BINARY_DIR}/D")
add_executable(initcond "${CMAKE_CURRENT_SOURCE_DIR}/examples/case/initcond.F90")
add_executable(tube "${CMAKE_CURRENT_SOURCE_DIR}/examples/case/tube.F90")
target_link_libraries(initcond common)
target_include_directories(tube PRIVATE ${PETSC_INCLUDE_DIRS})
target_link_libraries(tube common ${PETSC_LINK_LIBRARIES})
