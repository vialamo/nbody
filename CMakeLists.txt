cmake_minimum_required(VERSION 3.16)
project(nbody LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# ---------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Eigen3 3.3 REQUIRED)
find_package(OpenMP REQUIRED)

# ---------------------------------------------------------------------
# Core Simulation Library
# ---------------------------------------------------------------------
set(CORE_SOURCES
    src/engine.cpp
    src/gas.cpp
    src/ics.cpp
    src/integrator.cpp
    src/particles.cpp
    src/types.cpp
    src/utils.cpp
)

add_library(SimulationCore ${CORE_SOURCES})

target_link_libraries(SimulationCore
    PUBLIC
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# Tell the compiler where to find your .h files!
target_include_directories(SimulationCore 
    PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${HDF5_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

# ---------------------------------------------------------------------
# Main Executable
# ---------------------------------------------------------------------
add_executable(nbody src/main.cpp)

target_link_libraries(nbody PRIVATE SimulationCore)

if(MSVC)
    target_compile_options(nbody PRIVATE /wd4251)
endif()

# ---------------------------------------------------------------------
# Testing Infrastructure (Catch2)
# ---------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.3 
)
FetchContent_MakeAvailable(Catch2)

enable_testing()

# The test runner looks in the tests/ folder
add_executable(run_tests tests/test_physics.cpp)

target_link_libraries(run_tests PRIVATE SimulationCore Catch2::Catch2WithMain)
# Tests need to include headers from src/
target_include_directories(run_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

include(CTest)
include(Catch)
catch_discover_tests(run_tests)

# ---------------------------------------------------------------------
# Post-build: copy resource files
# ---------------------------------------------------------------------
set(RESOURCE_FILES
    ${CMAKE_SOURCE_DIR}/simulation.ini
)

foreach(RESOURCE_FILE ${RESOURCE_FILES})
    if(EXISTS ${RESOURCE_FILE})
        add_custom_command(
            TARGET nbody POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${RESOURCE_FILE}
                $<TARGET_FILE_DIR:nbody>/
            COMMENT "Copying ${RESOURCE_FILE} to executable directory"
        )
    endif()
endforeach()

# ---------------------------------------------------------------------
# Clang-Tidy
# ---------------------------------------------------------------------
find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
if(CLANG_TIDY_EXE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set_target_properties(SimulationCore PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    set_target_properties(nbody PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
else()
    message(WARNING "clang-tidy not found!")
endif()