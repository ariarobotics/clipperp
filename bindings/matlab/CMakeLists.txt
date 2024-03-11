find_package(Matlab COMPONENTS MX_LIBRARY)

if(Matlab_FOUND)
    message(STATUS "MATLAB root directory found: ${Matlab_ROOT_DIR}")

    matlab_add_mex(NAME clipperplus_clique_mex SRC clipperplus_clique_mex.cpp LINK_TO Eigen3::Eigen clipperplus)
    set_target_properties(clipperplus_clique_mex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
    set_target_properties(clipperplus_clique_mex PROPERTIES OUTPUT_NAME clipperplus_clique_mex)

    matlab_add_mex(NAME clique_optimization_mex SRC clique_optimization_mex.cpp LINK_TO Eigen3::Eigen clipperplus)
    set_target_properties(clique_optimization_mex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
    set_target_properties(clique_optimization_mex PROPERTIES OUTPUT_NAME clique_optimization_mex)

    matlab_add_mex(NAME clique_corenumber_mex SRC clique_corenumber_mex.cpp LINK_TO Eigen3::Eigen clipperplus)
    set_target_properties(clique_corenumber_mex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
    set_target_properties(clique_corenumber_mex PROPERTIES OUTPUT_NAME clique_corenumber_mex)


    if(OPTION_BUILD_PMC_HEU)
        matlab_add_mex(NAME pmc_heuristic_mex SRC pmc_heuristic_mex.cpp LINK_TO Eigen3::Eigen pmcheu)
        set_target_properties(pmc_heuristic_mex PROPERTIES COMPILE_FLAGS "-fvisibility=default")
        set_target_properties(pmc_heuristic_mex PROPERTIES OUTPUT_NAME pmc_heuristic_mex)
    endif()

else()
    message(WARNING "MATLAB root directory not found. Will not build MATLAB bindings.")
    set(OPTION_BUILD_BINDINGS_MATLAB OFF)
endif()