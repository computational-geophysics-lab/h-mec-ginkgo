# Find all .cpp files in the current directory
file(GLOB PROGRAM_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Extract the base name of the source file(s) without the extension
get_filename_component(PROGRAM_NAME ${PROGRAM_SOURCE} NAME_WE)

# Add the executable for the program
add_executable(${PROGRAM_NAME} ${PROGRAM_SOURCE})
target_link_libraries(${PROGRAM_NAME} Ginkgo::ginkgo)

# Add tests if test/ folder exist
if(BUILD_TESTS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_tests_for_subdirectory(${PROGRAM_NAME})
endif()

# Check if the data folder exists and copy the contents if it does
set(DATA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(DATA_DESTINATION_DIR "${CMAKE_CURRENT_BINARY_DIR}/data")

# Copy data over if data/ folder exist
if(EXISTS ${DATA_SOURCE_DIR})
    # Create the destination directory in case it doesn't exist
    file(MAKE_DIRECTORY ${DATA_DESTINATION_DIR})
    # Find all files in the data directory
    file(GLOB DATA_FILES "${DATA_SOURCE_DIR}/*")
    # Copy each file found to the build directory
    foreach(DATA_FILE ${DATA_FILES})
        get_filename_component(FILE_NAME ${DATA_FILE} NAME)
        configure_file(${DATA_FILE} ${DATA_DESTINATION_DIR}/${FILE_NAME} COPYONLY)
    endforeach()
endif()
