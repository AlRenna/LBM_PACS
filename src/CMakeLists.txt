cmake_minimum_required(VERSION 3.10)

# Include FetchContent module
include(FetchContent)

# find openmp package
find_package(OpenMP REQUIRED)

# Specify the CUDA toolkit path
set(CUDAToolkit_ROOT /usr/local/cuda-12.6)

find_package(CUDAToolkit)
if (CUDAToolkit_FOUND)
  message(STATUS "CUDA Toolkit found at ${CUDAToolkit_ROOT}")
else()
  message(FATAL_ERROR "CUDA Toolkit not found")
endif()

if (CUDAToolkit_FOUND)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_COMPILER "${CUDAToolkit_ROOT}/bin/nvcc")
  enable_language(CUDA)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  include_directories(${CUDA_INCLUDE_DIRS})
  file(GLOB_RECURSE CUDA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cu" "${CMAKE_CURRENT_SOURCE_DIR}/*.cuh")
  add_library(LBMCudaLib ${CUDA_SOURCES})
  target_link_libraries(LBMCudaLib PUBLIC ${CUDAToolkit_LIBRARIES} OpenMP::OpenMP_CXX)
  
  set_target_properties(LBMCudaLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_compile_definitions(LBMCudaLib PRIVATE -DUSE_CUDA)
endif()

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

if (CUDAToolkit_FOUND)
  target_link_libraries(LatticeLib PUBLIC LBMCudaLib)
endif()

# Link nlohmann-json and OpenMP
target_link_libraries(LatticeLib PUBLIC nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

# Add the root directory to the include directories
target_include_directories(LatticeLib PUBLIC ${CMAKE_SOURCE_DIR})
