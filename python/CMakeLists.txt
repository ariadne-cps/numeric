set(Python3_FIND_VIRTUALENV FIRST)
set(Python3_FIND_STRATEGY LOCATION)
set(PYBIND11_FINDPYTHON ON)

find_package(Python3 3.11 COMPONENTS Interpreter Development.Module)

option(PYTHON_ENFORCE_TESTS "Enforce Python tests" OFF)

if(Python3_FOUND)
    message(STATUS "Python executable: ${Python3_EXECUTABLE}")
    message(STATUS "Python version:    ${Python3_VERSION}")
    message(STATUS "Python includes:   ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python SOABI:      ${Python3_SOABI}")

    add_subdirectory(source)

    if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
        pybind11_add_module(pyariadne-numeric
            $<TARGET_OBJECTS:pyariadne-numeric-module-obj>
            $<TARGET_OBJECTS:UTILITY_SRC>
            $<TARGET_OBJECTS:PARADIGM_SRC>
            $<TARGET_OBJECTS:NUMERIC_SRC>
        )
        set_target_properties(pyariadne-numeric PROPERTIES OUTPUT_NAME "pyariadne_numeric")
        target_link_libraries(pyariadne-numeric PRIVATE pyariadne-numeric-python numeric)

        add_subdirectory(tests)
        add_dependencies(everything pyariadne-numeric)
    endif()
else()
    if(PYTHON_ENFORCE_TESTS)
        message(FATAL_ERROR "Python >= 3.11 with Interpreter and Development.Module components is required when PYTHON_ENFORCE_TESTS=ON.")
    else()
        message(STATUS "Compatible Python development environment not found: Python bindings will not be built.")
    endif()
endif()
