cmake_minimum_required(VERSION 3.10)

# Collect all source files in the test_lid_driven directory
file(GLOB_RECURSE TEST_INFLOW_OUTFLOW_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Create an executable from the source files
add_executable(test_inflow_outflow ${TEST_INFLOW_OUTFLOW_SOURCES})

# Link the executable with the main library
target_link_libraries(test_inflow_outflow PRIVATE LatticeLib)

find_package(CUDA)

if (CUDA_FOUND)
  target_link_libraries(test_inflow_outflow PRIVATE LBMCudaLib)
  target_compile_definitions(test_inflow_outflow PRIVATE -DUSE_CUDA)
  message(STATUS "setting USE_CUDA")
endif()

# Copy params.json to the build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/params.json DESTINATION ${CMAKE_BINARY_DIR}/examples/test_inflow_outflow)
