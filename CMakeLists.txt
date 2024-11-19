# Required minimum version of CMake
cmake_minimum_required(VERSION 3.16)

# Project name
project(H-MECs-Ginkgo)

#Set C++17 Standard
set(CMAKE_CXX_STANDARD 17)

# Include the custom SetupTests.cmake file
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SetupDependency)

# Configuration options for setting up tests for each subdirectory
option(BUILD_BENCHMARKS "Enable Google Benchmarking" OFF)
option(BUILD_TESTS "Generate build files for unit tests" OFF)
option(VISU_OPENCV "Visualization using OpenCV" ON)

# Configure Google Test v1.12.0
if(BUILD_TESTS)
    enable_testing()
    # Find and include Google Test
    find_or_fetch_library("GTest" "f8d7d77c06936315286eb55f8de22cd23c188571" "1.12.0" "https://github.com/google/googletest.git")

    # Find and include Google Test
    include(GoogleTest)  
    include(SetupTests)
endif()

#Find Ginkgo
#list(APPEND CMAKE_MODULE_PATH "/home/niklmeier/ginkgo/build")
#list(APPEND CMAKE_PREFIX_PATH "/home/niklmeier/ginkgo/build")
find_package(Ginkgo)
if(Ginkgo_FOUND)
    message(STATUS "Ginkgo version: ${Ginkgo_VERSION}")
endif()

# Configure Benchmark v1.8.3
if(BUILD_BENCHMARKS)
    find_or_fetch_library("Benchmark" "344117638c8ff7e239044fd0fa7085839fc03021" "1.8.3" "https://github.com/google/benchmark.git")
endif()

# Configure Visualization with OpenCV v4.6.0

if(VISU_OPENCV)
    list(APPEND CMAKE_MODULE_PATH "/cluster/home/   niklmeier/opencv/build")
    list(APPEND CMAKE_PREFIX_PATH "/cluster/home/niklmeier/opencv/build")
    find_package(OpenCV)
    #find_or_fetch_library("OpenCV" "b0dc474160e389b9c9045da5db49d03ae17c6a6b" "4.6.0" "https://github.com/opencv/opencv.git")
endif()


# Find and include HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(HDF5)

# Directory containing the codes
add_subdirectory(experiments)
