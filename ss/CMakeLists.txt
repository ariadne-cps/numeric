# Foundation supplies Utility transitively. Logging is not a Numeric
# dependency: optional Benchmarking supplies it directly.
if(NOT TARGET foundation)
    ariadne_add_submodule(foundation)
endif()
