cmake_minimum_required(VERSION 3.10)

# Collect all source and header files in the numerics directory
file(GLOB_RECURSE UTILS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

# Create a library from the source files
add_library(UtilsLib ${Utils_SOURCES})

# Link the library with the main library
target_link_libraries(UtilsLib LatticeLib)