# Function to simplify adding tests for a subdirectory
function(add_tests_for_subdirectory PROGRAM_NAME)
    file(GLOB TEST_SOURCES "test/*-test.cpp")

    # Create a test executable for the specified subdirectory
    add_executable(${PROGRAM_NAME}-tests ${TEST_SOURCES})

    # Link Ginkgo and Google Test to the test executable
    target_link_libraries(${PROGRAM_NAME}-tests
        Ginkgo::ginkgo
        GTest::GTest
        GTest::Main
    )

    # Register the tests with CTest for automatic test discovery
    gtest_discover_tests(${PROGRAM_NAME}-tests)
    
endfunction()