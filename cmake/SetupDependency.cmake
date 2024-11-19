include(FetchContent)

# Define helper
function(find_or_fetch_library lib_name lib_tag version url)
    find_package(${lib_name} ${version} QUIET)
    if (NOT ${lib_name}_FOUND)
        message(STATUS "Fetching ${lib_name} ${lib_tag} ${version}")
        
        set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/third_party/${lib_name}")
        
        FetchContent_Declare(
            ${lib_name}
            GIT_REPOSITORY ${url}
            GIT_TAG ${lib_tag}
        )
        
        FetchContent_MakeAvailable(${lib_name})
    else()
        message(STATUS "${lib_name} ${version} found in the system.")
    endif()
endfunction()