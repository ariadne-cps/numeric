if(NOT TARGET foundation)
    ensure_submodule_initialized(foundation)
    add_subdirectory(foundation)
endif()
