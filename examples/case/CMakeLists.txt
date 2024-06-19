# get current subdirectory name
get_filename_component(current_dir ${CMAKE_CURRENT_LIST_DIR} NAME)
message(STATUS "Current directory: ${current_dir}")

file(GLOB case_F90s "*.F90")
cmake_print_variables(case_F90s)

# make subdirectory a custom target with F90 files as dependencies
# so `make case` will make all F90 files inside the case directory
add_custom_target(${current_dir})

foreach(F90_filepath ${case_F90s})
    # get just the name of the file from the filepath for the targets
    get_filename_component(file_target ${F90_filepath} NAME_WE)
    add_executable(${file_target} ${F90_filepath})
    cmake_print_variables(file_target)

    # give all F90 files the same include and link options
    # initcond doesn't need this, but it can't hurt
    target_include_directories(${file_target} PRIVATE ${PETSC_INCLUDE_DIRS})
    target_link_libraries(${file_target} PRIVATE common ${PETSC_LINK_LIBRARIES})

    add_dependencies(${current_dir} ${file_target})
endforeach()
