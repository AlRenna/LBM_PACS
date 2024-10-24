cmake_minimum_required(VERSION 3.10)

# Collect all source and header files in the src directory
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

# Create a library from the source files
add_library(LatticeLib ${SOURCES})

# Optionally, specify include directories if needed
# target_include_directories(LatticeLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add subdirectories for further libraries or executables
# add_subdirectory(io)
# add_subdirectory(numerics)
# add_subdirectory(python_scripts)