#Find all .cpp files in current directory
file(GLOB PROGRAM_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Extract the name without the extension for the executable name
foreach(source_file ${PROGRAM_SOURCE})
	get_filename_component(PROGRAM_NAME ${source_file} NAME_WE)
	get_filename_component(PROGRAM_DIRECTORY ${source_file} DIRECTORY)
	file(GLOB INCLUDE_FILES "${PROGRAM_DIRECTORY}/include/*.cpp")
	add_executable(${PROGRAM_NAME} ${source_file} ${INCLUDE_FILES}
			include/constants.cpp)
endforeach()

target_include_directories(${PROGRAM_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})

#Choose accelerator:
target_compile_options(${PROGRAM_NAME}
        PUBLIC	-DGINKGO_BUILD_REFERENCE=ON
        PUBLIC 	-DGINKGO_BUILD_CUDA=ON
#       PUBLIC	-DGINKGO_BUILD_HIP=OFF
#       PUBLIC  -DGINKGO_BUILD_SYCL=OFF
)

#Link the libraries:
target_link_libraries(${PROGRAM_NAME}
		Ginkgo::ginkgo
		hdf5::hdf5
		hdf5::hdf5_cpp
		${HDF5_CXX_HL_LIBRARIES}
		${HDF5_C_LIBRARIES}
		${HDF5_CPP_LIBRARIES}
		${HDF5_HL_CPP_LIBRARIES}
		${HDF5_LIBRARIES} 
		${HDF5_CXX_LIBRARIES}
		${HDF5_HL_LIBRARIES}
)

