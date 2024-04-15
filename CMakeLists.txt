cmake_minimum_required(VERSION 3.10)
project(clipperplus VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type." FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

###############################################################################
# Options
###############################################################################

option(DEBUG_FLAG "Enable debugging output" OFF)
option(DEBUG_TIMING_FLAG "Enable timing report" OFF)
option(DEBUG_OPTIM_FLAG "Enable debugging clipper optimization" OFF)

option(BUILD_PMC_HEU "Build PMC heuristic" OFF)
option(BUILD_BINDINGS_MATLAB "Build MATLAB bindings" OFF)
option(BUILD_BINDINGS_PYTHON "Build Python bindings" OFF)
option(BUILD_TESTS "Build testsuite" ON)


###############################################################################
# Dependencies
###############################################################################

include(FetchContent)

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen Version: ${EIGEN3_VERSION_STRING} ${EIGEN3_VERSION}")

find_package(OpenMP REQUIRED)

if(BUILD_TESTS)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/e4fdb87e76b9fc4b01c54ad81aea19d6e994b994.zip
  )
endif()

if(BUILD_BINDINGS_PYTHON)
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.9.2
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


###############################################################################
# Targets
###############################################################################

set(PMC_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_heu.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_cores.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_graph.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_clique_utils.cpp )

set(CLIPPER_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/clipper.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/utils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/dsd.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/invariants/euclidean_distance.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external/clipper/invariants/pointnormal_distance.cpp )
        
set(CLLIPPER_PLUS_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clipperplus_clique.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clipperplus_graph.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clique_corenumber.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clique_optimization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/clipperplus_heuristic.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp
)

# source files
add_library(clipperplus SHARED)
target_sources(clipperplus PRIVATE 
        ${PMC_SOURCE_FILES}
        ${CLIPPER_SOURCE_FILES}        
        ${CLLIPPER_PLUS_SOURCE_FILES} )

# headers:
target_include_directories(clipperplus PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

# link libraries:
target_link_libraries(clipperplus PUBLIC Eigen3::Eigen)
target_link_libraries(clipperplus PUBLIC OpenMP::OpenMP_CXX)



# debug flag  
if(DEBUG_FLAG)
        message(STATUS "Enabling debug flag")
        target_compile_definitions(clipperplus PRIVATE DEBUG)
endif()
if(DEBUG_OPTIM_FLAG)
        message(STATUS "Enabling debug optimization flag")
        target_compile_definitions(clipperplus PRIVATE DEBUG_OPTIM)
endif()
if(DEBUG_TIMING_FLAG)
        message(STATUS "Enabling debug timing flag")
        target_compile_definitions(clipperplus PRIVATE DEBUG_TIMING)
endif()

# version
set_target_properties(clipperplus PROPERTIES
  VERSION ${PROJECT_VERSION})
target_compile_definitions(clipperplus PUBLIC CLIPPERPLUS_VERSION="${PROJECT_VERSION}")


###############################################################################
# Extras
###############################################################################

if(BUILD_BINDINGS_MATLAB)
    message(STATUS "Attempting to build MATLAB bindings.")
    add_subdirectory(bindings/matlab)
endif()

if(BUILD_BINDINGS_PYTHON)
    message(STATUS "Building Python bindings.")
    add_subdirectory(bindings/python)
endif()


if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()


# building PMC heuristic algorithm
if(BUILD_PMC_HEU)
        message(STATUS "Building PMC heuristic.")
        set(PMC_HUE_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/external/pmc/pmc_heuristic.cpp )

        add_library(pmcheu SHARED)
        target_sources(pmcheu PRIVATE 
                ${PMC_SOURCE_FILES}
                ${PMC_HUE_SOURCE_FILES} )

        # headers:
        target_include_directories(pmcheu PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>)

        # link libraries:
        target_link_libraries(pmcheu PUBLIC
                OpenMP::OpenMP_CXX
                Eigen3::Eigen)

        # debug flag  
        if(DEBUG_FLAG)
                target_compile_definitions(pmcheu PRIVATE DEBUG)
        endif()
endif()
