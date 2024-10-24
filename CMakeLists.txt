cmake_minimum_required(VERSION 3.10)

# Project name
project(LatticeBoltzmann_PACS)

# Function to add all subdirectories.
function(add_all_subdirectories)
    file(GLOB SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
    foreach(subdir ${SUBDIRS})
        if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
            add_subdirectory(${subdir})
        endif()
    endforeach()
endfunction()

# Generate documentation.
# add_subdirectory(${CMAKE_SOURCE_DIR}/doc/)

# Generate source code.
add_subdirectory(${CMAKE_SOURCE_DIR}/src/)

# Generate examples.
add_subdirectory(${CMAKE_SOURCE_DIR}/examples/)

# Add python_scripts subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/src/python_scripts)

