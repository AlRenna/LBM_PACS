cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

# add_compile_options(-Wall -Wextra)

# Project name
project(LatticeBoltzmann_PACS LANGUAGES CXX) # CUDA)

# Function to add all subdirectories.
function(add_all_subdirectories)
    file(GLOB SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
    foreach(subdir ${SUBDIRS})
        if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
            add_subdirectory(${subdir})
        endif()
    endforeach()
endfunction()

# Generate source code.
add_subdirectory(${CMAKE_SOURCE_DIR}/src/)

# Generate examples.
add_subdirectory(${CMAKE_SOURCE_DIR}/examples/)

# Generate images.
add_subdirectory(${CMAKE_SOURCE_DIR}/images/)

# Add python_scripts subdirectory
add_subdirectory(${CMAKE_SOURCE_DIR}/src/python_scripts)

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
endif ()