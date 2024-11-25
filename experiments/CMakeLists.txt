cmake_minimum_required(VERSION 3.15)

project(experiments)

# Initialize the list
set(EXPERIMENTS_LIST "")

# Conditionally add Ginkgo experiments
if(BUILD_GINKGO)
    list(APPEND EXPERIMENTS_LIST h-mec-rsf-ginkgo-v1)
endif()

# Conditionally add Eigen experiments
if(BUILD_EIGEN)
    list(APPEND EXPERIMENTS_LIST h-mec-rsf-eigen-v2)
endif()

if(BUILD_MKL)
    list(APPEND EXPERIMENTS_LIST
    h-mec-rsf-eigen-v8
    h-mec-rsf-eigen-v10)
endif()

# Print the final list for debugging
message(STATUS "Experiments to build: ${EXPERIMENTS_LIST}")

# Build for all examples
foreach(experiment ${EXPERIMENTS_LIST})
    add_subdirectory(${experiment})
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}")
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}/output_data")
	#Write a file called StartingTimestep.txt in which to specify if the simulation starts new or continues from a previous simulation
endforeach()
