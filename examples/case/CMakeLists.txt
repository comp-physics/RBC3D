file(GLOB case_F90s "*.F90")
cmake_print_variables(case_F90s)

foreach(F90_filepath ${case_F90s})
    get_filename_component(file_target ${F90_filepath} NAME_WE)
    add_executable(${file_target} ${F90_filepath})
    cmake_print_variables(file_target)

    target_include_directories(${file_target} PRIVATE ${PETSC_INCLUDE_DIRS})
    target_link_libraries(${file_target} PRIVATE common ${PETSC_LINK_LIBRARIES})
    message("HERE")
endforeach()

# add_executable(initcond "initcond.F90")
# add_executable(tube "tube.F90")

# target_link_libraries(initcond PUBLIC common)
# target_include_directories(tube PRIVATE ${PETSC_INCLUDE_DIRS})
# target_link_libraries(tube common ${PETSC_LINK_LIBRARIES})