cmake_minimum_required(VERSION 3.20)
cmake_policy(SET CMP0135 NEW)
project(NurbsLab VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile_commands.json for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(NURBSLAB_BUILD_TESTS "Build unit tests" ON)
option(NURBSLAB_BUILD_EXAMPLES "Build examples" ON)

# Library sources
set(NURBSLAB_SOURCES
    src/knot_vector.cpp
    src/nurbs_curve.cpp
    src/algorithms/projection.cpp
    src/algorithms/intersection.cpp
)

# Create library
add_library(nurbslab ${NURBSLAB_SOURCES})

target_include_directories(nurbslab
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(nurbslab PUBLIC cxx_std_20)

# Compiler warnings
if(MSVC)
    target_compile_options(nurbslab PRIVATE /W4)
else()
    target_compile_options(nurbslab PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Tests
if(NURBSLAB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Examples
if(NURBSLAB_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()