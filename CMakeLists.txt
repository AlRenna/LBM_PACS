cmake_minimum_required(VERSION 3.10)

# Set the project name
project(LatticeSimulation)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the src directory to the include path
include_directories(src)

# Find all .cpp and .hpp files in the src directory
file(GLOB_RECURSE SRC_FILES src/*.cpp src/*.hpp)

# Find all .cpp files in the examples directory
file(GLOB EXAMPLES_FILES examples/*.cpp)

# Add the executable
add_executable(LatticeSimulation ${EXAMPLES_FILES} ${SRC_FILES})

# Link libraries if needed (e.g., Python libraries)
# target_link_libraries(LatticeSimulation <library_name>)