cmake_minimum_required(VERSION 3.10)

# Collect all source files in the test1 directory
file(GLOB_RECURSE TEST1_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Create an executable from the source files
add_executable(test1 ${TEST1_SOURCES})

# Link the executable with the main library
target_link_libraries(test1 LatticeLib)