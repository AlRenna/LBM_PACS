cmake_minimum_required(VERSION 3.10)

# Include FetchContent module
include(FetchContent)

# find openmp package
find_package(OpenMP REQUIRED)

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

# Link nlohmann-json and OpenMP to your target
target_link_libraries(LatticeLib PUBLIC nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

# Add the root directory to the include directories
target_include_directories(LatticeLib PUBLIC ${CMAKE_SOURCE_DIR})
