cmake_minimum_required(VERSION 3.20)

# Find OpenMP package
find_package(OpenMP REQUIRED)

# Include FetchContent module
include(FetchContent)

# Download nlohmann-json library
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.10.5 # You can specify the version you need
)
FetchContent_MakeAvailable(nlohmann_json)

# Collect all source and header files in the src directory
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

# Create a library from the source files
add_library(LatticeLib ${SOURCES})

# Link nlohmann-json and OpenMP
target_link_libraries(LatticeLib PUBLIC nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

# Check if CUDA is available and create custom CUDA library for LBM
if (CMAKE_CUDA_COMPILER)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  include_directories(${CUDA_INCLUDE_DIRS})
  add_library(LBMCudaLib lattice_gpu.cu)
  target_link_libraries(LBMCudaLib PRIVATE ${CUDA_LIBRARIES} OpenMP::OpenMP_CXX)
  target_compile_definitions(LBMCudaLib PRIVATE -DUSE_CUDA)
endif()

# Link CUDA if available
if (CMAKE_CUDA_COMPILER)
  target_link_libraries(LatticeLib PRIVATE LBMCudaLib)
  target_compile_definitions(LatticeLib PRIVATE -DUSE_CUDA)
endif()

# Add the root directory to the include directories
target_include_directories(LatticeLib PUBLIC ${CMAKE_SOURCE_DIR})
