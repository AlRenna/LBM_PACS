cmake_minimum_required(VERSION 3.10)

# Collect all Python files in the python_scripts directory
file(GLOB PYTHON_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

# Specify the destination directory in the build folder
set(PYTHON_SCRIPTS_DESTINATION "${CMAKE_BINARY_DIR}/src/python_scripts")

# Copy Python scripts to the build directory
file(COPY ${PYTHON_SCRIPTS} DESTINATION ${PYTHON_SCRIPTS_DESTINATION})

# Add a custom target to ensure the scripts are copied during the build process
add_custom_target(copy_python_scripts ALL DEPENDS ${PYTHON_SCRIPTS})