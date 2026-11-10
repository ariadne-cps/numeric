cmake_minimum_required(VERSION 3.15)

project(Numeric VERSION 3.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
if(NOT TARGET ariadne-build)
    find_package(Git REQUIRED)
    # Keep initial Build checkout messages in English.
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env "LC_ALL=C"
                "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}"
                submodule update --init -- ss/bui
        RESULT_VARIABLE _result)
    if(_result)
        message(FATAL_ERROR "Cannot initialize build submodule")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/ss/bui/cmake/ProjectOptions.cmake")
endif()

require_same_dependency_commit(Build CHAINS
    "build" "foundation|build")

# Numeric reaches Utility via Foundation, while Benchmarking reaches it via Logging.
require_same_dependency_commit(Utility CHAINS
    "foundation|utility" "benchmarking|logging|utility")

ariadne_setup_module(numeric)

# Keep Numeric's ABI and strict floating-point compilation policy.
if(WIN32)
    target_compile_definitions(numeric INTERFACE ARIADNE_NUMERIC_STATIC)
endif()

if(MSVC)
    target_compile_options(NUMERIC_SRC PRIVATE /fp:strict)
else()
    target_compile_options(NUMERIC_SRC PRIVATE -frounding-math)
endif()
