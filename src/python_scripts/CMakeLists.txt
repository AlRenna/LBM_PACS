cmake_minimum_required(VERSION 3.10)

# Collect all Python files in the python_scripts directory
file(GLOB PYTHON_SCRIPTS "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

# Specify the destination directory in the build folder
set(PYTHON_SCRIPTS_DESTINATION "${CMAKE_BINARY_DIR}/src/python_scripts")

# Create symbolic links for Python scripts in the build directory
foreach(PYTHON_SCRIPT ${PYTHON_SCRIPTS})
    get_filename_component(SCRIPT_NAME ${PYTHON_SCRIPT} NAME)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${PYTHON_SCRIPT} ${PYTHON_SCRIPTS_DESTINATION}/${SCRIPT_NAME})
endforeach()
