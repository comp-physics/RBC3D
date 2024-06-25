file(GLOB common_F90s "*.F90")
add_library(common STATIC "${common_F90s}")

target_include_directories(common PRIVATE ${PETSC_INCLUDE_DIRS} ${NETCDF_INCLUDE} ${FFTW_INCLUDE})
target_link_libraries(common ${LAPACK_LIBRARIES} ${PETSC_LINK_LIBRARIES} ${SPHPK_LIB} ${FFTW_LIB} ${NETCDF_LIB})

