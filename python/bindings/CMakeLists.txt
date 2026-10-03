include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)
get_target_property(FOUNDATION_SOURCE_DIR foundation SOURCE_DIR)

add_library(pyariadne-numeric-submodule-obj OBJECT
    numeric_submodule.cpp
)
target_include_directories(pyariadne-numeric-submodule-obj PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    "${FOUNDATION_SOURCE_DIR}/python/bindings"
)
target_link_libraries(pyariadne-numeric-submodule-obj PRIVATE
    pyariadne-python-common
    pybind11::module
    numeric
)
set_target_properties(pyariadne-numeric-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(pyariadne-numeric-submodule-obj PRIVATE /bigobj)
endif()

add_library(pyariadne-numeric INTERFACE)
target_sources(pyariadne-numeric INTERFACE
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
)
target_include_directories(pyariadne-numeric INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(pyariadne-numeric INTERFACE
    pyariadne-foundation
)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    add_library(pyariadne-numeric-module-obj OBJECT
        numeric_module.cpp
    )
    target_include_directories(pyariadne-numeric-module-obj PRIVATE
        "${PYTHON_COMMON_SOURCE_DIR}/include"
    )
    target_link_libraries(pyariadne-numeric-module-obj PRIVATE
        pyariadne-python-common
        pybind11::module
    )
    set_target_properties(pyariadne-numeric-module-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
