cmake_minimum_required(VERSION 3.15)

project(experiments)

# List of example folder (add more if needed)
set(EXPERIMENTS_LIST
    h-mec-rsf-ginkgo-v1
    h-mec-rsf-eigen-v2
    # h-mec-rsf-eigen-v8	# pardiso solver requires intelMKL
    # h-mec-rsf-eigen-v10	# pardiso solver requires intelMKL
)

# Build for all examples
foreach(experiment ${EXPERIMENTS_LIST})
    add_subdirectory(${experiment})
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}")
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}/output_data")
	#Write a file called StartingTimestep.txt in which to specify if the simulation starts new or continues from a previous simulation
endforeach()
