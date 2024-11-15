cmake_minimum_required(VERSION 3.10)

# Collect all .png files in the images directory
file(GLOB IMAGES "${CMAKE_CURRENT_SOURCE_DIR}/*.png")

# Specify the destination directory in the build folder
set(IMAGES_DESTINATION "${CMAKE_BINARY_DIR}/images")

# Copy .png files to the build directory
file(COPY ${IMAGES} DESTINATION ${IMAGES_DESTINATION})

# Add a custom target to ensure the scripts are copied during the build process
add_custom_target(copy_images ALL DEPENDS ${IMAGES})