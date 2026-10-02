cmake_minimum_required(VERSION 3.15)

project(Numeric VERSION 1.0)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)

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
            "${PROJECT_SOURCE_DIR}/submodules/paradigm|submodules/configuration"
    )
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/include)

if(NOT TARGET paradigm)
    add_subdirectory(submodules)
endif()

add_library(numeric INTERFACE)
target_include_directories(numeric INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/ariadne-numeric>
)
target_link_libraries(numeric INTERFACE
    paradigm
    GMP::GMP
    MPFR::MPFR
)

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
        PARADIGM_SRC
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

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python)
endif()

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    get_target_property(PARADIGM_SOURCE_DIR paradigm SOURCE_DIR)
    get_target_property(UTILITY_SOURCE_DIR utility SOURCE_DIR)

    install(FILES "${PROJECT_SOURCE_DIR}/include/ariadne-numeric.hpp"
        DESTINATION include/ariadne-numeric
    )
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/numeric"
        DESTINATION include/ariadne-numeric
    )
    install(DIRECTORY "${PARADIGM_SOURCE_DIR}/include/paradigm"
        DESTINATION include/ariadne-numeric
    )
    install(DIRECTORY "${UTILITY_SOURCE_DIR}/include/utility"
        DESTINATION include/ariadne-numeric
    )
endif()
