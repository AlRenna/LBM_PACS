cmake_minimum_required(VERSION 3.20)

# Collect all .png files in the images directory
file(GLOB IMAGES "${CMAKE_CURRENT_SOURCE_DIR}/*.png")

# Specify the destination directory in the build folder
set(IMAGES_DESTINATION "${CMAKE_BINARY_DIR}/images")

# Copy .png files to the build directory
file(COPY ${IMAGES} DESTINATION ${IMAGES_DESTINATION})

# Add a custom command to update the images folder if new files are added
add_custom_target(update_images ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${IMAGES_DESTINATION}
    DEPENDS ${IMAGES}
    COMMENT "Updating images folder"
)
