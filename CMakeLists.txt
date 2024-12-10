cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

# add_compile_options(-Wall -Wextra)

# add_compile_options(-Wno-dev -Wno-deprecated-declarations)

# Project name
project(LatticeBoltzmann_PACS LANGUAGES CXX)

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

# Generate images.
add_subdirectory(${CMAKE_SOURCE_DIR}/images/)

# Add python_scripts subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/src/python_scripts)