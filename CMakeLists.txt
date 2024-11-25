# Required minimum version of CMake
cmake_minimum_required(VERSION 3.16)

# Project name
project(H-MECs-Ginkgo)

#Set C++17 Standard
set(CMAKE_CXX_STANDARD 17)

# Options
option(BUILD_GINKGO "Build ginkgo version of the code" ON)
option(BUILD_EIGEN "Build eigen versions of the code, excluding ones that need MKL" ON)
option(BUILD_MKL "Build eigen version v10, depends on Intel oneMKL" OFF)


# Include the custom SetupTests.cmake file
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SetupDependency)

if(BUILD_EIGEN)
    # Find Eigen (required for eigen versions of the code)
    find_package(Eigen NO_MODULE)
endif()

if(BUILD_MKL)
    # Find Intel Math Kernel Library (required for h-mec-rsf-eigen-v10)
    find_package(MKL)
endif()


if(BUILD_GINKGO)
    #Find Ginkgo
    # Uncomment if using a local version of ginkgo
    #list(APPEND CMAKE_MODULE_PATH "/path/to/your/ginkgo/build")
    #list(APPEND CMAKE_PREFIX_PATH "/path/to/your/ginkgo/build")
    find_package(Ginkgo REQUIRED)
    if(Ginkgo_FOUND)
        message(STATUS "Ginkgo version: ${Ginkgo_VERSION}")
    endif()
endif()


# Find and include HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX HL)
find_package(HDF5)

# Directory containing the codes
add_subdirectory(experiments)
