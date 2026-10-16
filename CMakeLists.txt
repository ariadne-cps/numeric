cmake_minimum_required(VERSION 3.15)

project(Numeric VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
if(NOT TARGET ariadne-configuration)
    find_package(Git REQUIRED)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                    submodule update --init -- ss/configuration
                    RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize configuration submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/configuration/cmake/ProjectOptions.cmake")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    require_same_dependency_commit(Configuration
        CHAINS
            "configuration"
            "foundation|configuration"
    )
endif()

setup_project()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(ss)

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
if(WIN32)
    target_compile_definitions(numeric INTERFACE ARIADNE_NUMERIC_STATIC)
endif()

add_subdirectory(src)
target_link_libraries(NUMERIC_SRC PRIVATE numeric)

if(MSVC)
    target_compile_options(NUMERIC_SRC PRIVATE /fp:strict)
else()
    target_compile_options(NUMERIC_SRC PRIVATE -frounding-math)
endif()

ariadne_register_public_headers(
    TARGET numeric
    DIRECTORIES "${PROJECT_SOURCE_DIR}/include/numeric"
    FILES "${PROJECT_SOURCE_DIR}/include/ariadne-numeric.hpp"
)

setup_project_library(
    TARGET ariadne-numeric
    OBJECTS NUMERIC_SRC
    LINK_LIBRARIES numeric
)

# Benchmarks are registered by setup_standalone_project_tests().

setup_standalone_project_tests(
    TARGET ariadne-numeric
    TEST_DIRECTORY tests
    TEST_TARGET tests
    SOURCES
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    ariadne_install_dependency_bundle(
        TARGET numeric
        DESTINATION include/ariadne-numeric
    )
endif()
