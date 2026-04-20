cmake_minimum_required(VERSION 3.16)
project(nbody LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# ---------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------

# HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Eigen (header-only, so no linking required)
find_package(Eigen3 3.3 REQUIRED)

find_package(OpenMP REQUIRED)

# Include directories
include_directories(
    ${HDF5_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}
)

# ---------------------------------------------------------------------
# Target & Sources
# ---------------------------------------------------------------------
set(SOURCES
    main.cpp
    engine.cpp
    gas.cpp
    ics.cpp
    integrator.cpp
    particles.cpp
    utils.cpp
)

add_executable(nbody ${SOURCES})

target_link_libraries(nbody
    PRIVATE
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# Suppress the harmless HDF5 DLL warnings if compiling on Windows (MSVC)
if(MSVC)
    target_compile_options(nbody PRIVATE /wd4251)
endif()

# ---------------------------------------------------------------------
# Post-build: copy resource/config files next to the executable
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
# Status messages (for debug clarity)
# ---------------------------------------------------------------------
message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Eigen include dir: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Executable will be in: ${CMAKE_BINARY_DIR}")

# ---------------------------------------------------------------------
# Code suggestions
# ---------------------------------------------------------------------
find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
if(CLANG_TIDY_EXE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set_target_properties(nbody PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
else()
    message(WARNING "clang-tidy not found!")
endif()
