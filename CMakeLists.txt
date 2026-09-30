cmake_minimum_required(VERSION 3.15)

project(Numeric VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")
    message(FATAL_ERROR
        "Git submodules are not initialized. "
        "Run 'git submodule update --init --recursive' before configuring the project.")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/submodules/configuration/cmake/ProjectOptions.cmake")

require_initialized_submodules()
setup_project()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        REFERENCES
            "${PROJECT_SOURCE_DIR}|submodules/configuration"
            "${PROJECT_SOURCE_DIR}/submodules/foundation|submodules/configuration"
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include)

if(NOT TARGET foundation)
    add_subdirectory(submodules)
endif()

add_library(numeric INTERFACE)
target_include_directories(numeric INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-numeric>
)
target_link_libraries(numeric INTERFACE
    foundation
    GMP::GMP
    MPFR::MPFR
)

if(MSVC)
    target_compile_options(numeric INTERFACE
        /wd4099
        /wd4146
        /wd4244
        /wd4250
        /wd4267
        /wd4996
    )
endif()

add_subdirectory(src)
target_link_libraries(NUMERIC_SRC PRIVATE numeric)

if(MSVC)
    target_compile_options(NUMERIC_SRC PRIVATE /fp:strict)
else()
    target_compile_options(NUMERIC_SRC PRIVATE -frounding-math)
endif()

setup_project_library(
    TARGET ariadne-numeric
    OBJECTS
        NUMERIC_SRC
        FOUNDATION_SRC
        UTILITY_SRC
    LINK_LIBRARIES
        numeric
)

setup_standalone_project_tests(
    TARGET ariadne-numeric
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    get_target_property(FOUNDATION_SOURCE_DIR foundation SOURCE_DIR)
    get_target_property(UTILITY_SOURCE_DIR utility SOURCE_DIR)

    install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/numeric"
        DESTINATION include/ariadne-numeric
    )
    install(DIRECTORY "${FOUNDATION_SOURCE_DIR}/include/foundation"
        DESTINATION include/ariadne-numeric
    )
    install(DIRECTORY "${UTILITY_SOURCE_DIR}/include/utility"
        DESTINATION include/ariadne-numeric
    )
endif()
