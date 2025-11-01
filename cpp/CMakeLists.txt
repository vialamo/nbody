cmake_minimum_required(VERSION 3.16)
project(nbody_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

# ---------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------

# SFML
find_package(SFML 2.6.1 REQUIRED COMPONENTS graphics window system)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Eigen (header-only, so no linking required)
find_package(Eigen3 3.3 REQUIRED)

# Include directories
include_directories(
    ${HDF5_INCLUDE_DIRS}
    ${SFML_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}
)

# ---------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------
add_executable(nbody nbody.cpp)

target_link_libraries(nbody
    PRIVATE
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    sfml-graphics
    sfml-window
    sfml-system
)

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
message(STATUS "SFML include dir: ${SFML_INCLUDE_DIRS}")
message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Eigen include dir: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Executable will be in: ${CMAKE_BINARY_DIR}")

