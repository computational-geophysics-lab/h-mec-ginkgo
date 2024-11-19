cmake_minimum_required(VERSION 3.15)

project(experiments)

# List of example folder (add more if needed)
set(EXPERIMENTS_LIST
    h-mec-rsf-ginkgo-v1
    test_program
)

# Build for all examples
foreach(experiment ${EXPERIMENTS_LIST})
    add_subdirectory(${experiment})
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}")
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}/input_data")
	file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/experiments/${experiment}/output_data")
	#Write a file called StartingTimestep.txt in which to specify if the simulation starts new or continues from a previous simulation
	file(WRITE "${CMAKE_BINARY_DIR}/experiments/${experiment}/input_data/StartingTimestep.txt" "0")
endforeach()
