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
target_compile_definitions(pyariadne-numeric-submodule-obj PRIVATE
    ARIADNE_NUMERIC_EMBEDDED
)
set_target_properties(pyariadne-numeric-submodule-obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(pyariadne-numeric INTERFACE)
target_sources(pyariadne-numeric INTERFACE
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
)
target_include_directories(pyariadne-numeric INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<TARGET_PROPERTY:pyariadne-foundation,INTERFACE_INCLUDE_DIRECTORIES>
)
target_link_libraries(pyariadne-numeric INTERFACE
    pyariadne-foundation
)

ariadne_add_python_module(
    NAME numeric
    MODULE_SOURCE numeric_module.cpp
)
