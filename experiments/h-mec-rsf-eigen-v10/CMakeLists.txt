#Find all .cpp files in current directory
file(GLOB PROGRAM_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# Extract the name without the extension for the executable name
foreach(source_file ${PROGRAM_SOURCE})
	get_filename_component(PROGRAM_NAME ${source_file} NAME_WE)
	get_filename_component(PROGRAM_DIRECTORY ${source_file} DIRECTORY)
	file(GLOB INCLUDE_FILES "${PROGRAM_DIRECTORY}/include/*.cpp")
	add_executable(${PROGRAM_NAME} ${source_file} ${INCLUDE_FILES})
	
endforeach()


# Add compiler options
target_compile_definitions(${PROGRAM_NAME}
	PUBLIC -DEIGEN_USE_MKL_ALL
	PUBLIC _LARGEFILE_SOURCE
	PUBLIC _LARGEFILE64_SOURCE
	PUBLIC _BSD_SOURCE
	PUBLIC -DMKL_LP64
)

#Link the libraries:
target_link_libraries(${PROGRAM_NAME}
		MPI::MPI_C
# Uncomment when switching to Ginkgo
#		Ginkgo::ginkgo
		MKL::MKL
		${MKL_CORE_LIBRARY}
		${MKL_INTEL_THREAD_LIBRARY}
		${MKL_INTEL_LP64_LIBRARY}
		${HDF5_CXX_HL_LIBRARIES}
		${HDF5_C_LIBRARIES}
		${HDF5_CPP_LIBRARIES}
		${HDF5_HL_CPP_LIBRARIES}
		${HDF5_LIBRARIES} 
		${HDF5_CXX_LIBRARIES}
		${HDF5_HL_LIBRARIES}

		${MPI_CXX_LIBRARIES}
		-lsz
		-lz
		-lrt
		-lm
		-ldl
		-lpthread
)

#set_target_properties(${PROGRAM_NAME} PROPERTIES INSTALL_RPATH "/cluster/apps/hdf5/1.8.13/x86_64/gcc_4.8.2/serial/lib" RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
