cmake_minimum_required(VERSION 3.16)
project(asimov LANGUAGES CXX)

# ==============================================================================
# GLOBAL BUILD OPTIONS
# ==============================================================================
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(USE_GPU_OFFLOAD "Enable OpenMP GPU Offloading" ON)
option(USE_CLANG_TIDY  "Run clang-tidy during compilation" OFF)
option(USE_FFTW_OMP "Use FFTW3 OpenMP backend (disable to avoid libgomp clash)" OFF)
option(ENABLE_IPO "Enable Interprocedural Optimization (LTO)" ON)

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Eigen3 3.3 REQUIRED)
find_package(OpenMP REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)
find_library(FFTW3_LIBRARY NAMES fftw3 REQUIRED)
if(USE_FFTW_OMP)
    find_library(FFTW3_PARALLEL_LIBRARY NAMES fftw3_omp REQUIRED)
    message(STATUS "FFTW3 Parallel Backend: OpenMP (fftw3_omp)")
else()
    find_library(FFTW3_PARALLEL_LIBRARY NAMES fftw3_threads REQUIRED)
    message(STATUS "FFTW3 Parallel Backend: POSIX Threads (fftw3_threads)")
endif()

# ==============================================================================
# CORE LIBRARY (OBJECT TARGET)
# ==============================================================================
set(CORE_SOURCES
    src/config.cpp
    src/cooling.cpp
    src/cosmology.cpp
    src/debug_utils.cpp
    src/density.cpp
    src/diagnostics.cpp
    src/engine.cpp
    src/gas.cpp
    src/gravity.cpp
    src/hdf5_reader.cpp
    src/hdf5_writer.cpp
    src/hydro.cpp
    src/ics.cpp
    src/integrator.cpp
    src/lbvh.cpp
    src/logger.cpp
    src/math_utils.cpp
    src/mfm.cpp
    src/particles.cpp
    src/pm.cpp
    src/reconstruction.cpp
    src/riemann.cpp
    src/thermodynamics.cpp
    src/types.cpp
    src/utils.cpp
)

add_library(SimulationCore OBJECT ${CORE_SOURCES})

target_link_libraries(SimulationCore PUBLIC
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    OpenMP::OpenMP_CXX
    ${FFTW3_LIBRARY}
    ${FFTW3_PARALLEL_LIBRARY}
)

target_include_directories(SimulationCore PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

target_include_directories(SimulationCore SYSTEM PUBLIC
    ${HDF5_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
)
# ==============================================================================
# HARDWARE OPTIMIZATION
# ==============================================================================
# Base CPU Vectorization
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(SimulationCore PUBLIC 
        "-march=native"
        "-fno-math-errno"
    )
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Enabling Clang SIMD vectorization diagnostics")
        target_compile_options(SimulationCore PUBLIC 
            -Rpass-missed=loop-vectorize 
            -Rpass-analysis=loop-vectorize
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        message(STATUS "Enabling GCC SIMD vectorization diagnostics")
        target_compile_options(SimulationCore PUBLIC 
            "-fopt-info-vec-missed"
        )
    endif()
elseif(MSVC)
    target_compile_options(SimulationCore PUBLIC "/arch:AVX2")
endif()

# ==============================================================================
# GPU OFFLOADING PIPELINE
# ==============================================================================
if(USE_GPU_OFFLOAD)
    message(STATUS "GPU Offloading ENABLED")
    
    target_compile_definitions(SimulationCore PUBLIC USE_GPU)

    execute_process(
        COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
        OUTPUT_VARIABLE GPU_COMPUTE_CAP
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )

    if(GPU_COMPUTE_CAP)
        string(REGEX MATCH "^([0-9]+)\\.([0-9]+)" _ ${GPU_COMPUTE_CAP})
        set(COMPUTE_VERSION "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
    else()
        set(COMPUTE_VERSION "61")
        message(WARNING "Could not auto-detect GPU via nvidia-smi. Defaulting to architecture 61.")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" OR CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        set(NVC_ARCH "cc${COMPUTE_VERSION}")
        message(STATUS "Detected nvc++. Targeting NVIDIA GPU (OpenMP Offload): ${NVC_ARCH}")
        
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -mp=gpu -gpu=${NVC_ARCH} -Minfo=mp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mp=gpu -gpu=${NVC_ARCH}")

    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(GPU_ARCH "sm_${COMPUTE_VERSION}")

        message(STATUS "Applying Clang GPU Offload Flags (${GPU_ARCH})")
        target_compile_options(SimulationCore PUBLIC
            -fopenmp
            -fopenmp-targets=nvptx64-nvidia-cuda 
            --offload-arch=${GPU_ARCH}
        )
        add_link_options(
            -fopenmp
            -fopenmp-targets=nvptx64-nvidia-cuda
            -fuse-ld=lld
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(GPU_ARCH "sm_${COMPUTE_VERSION}")

        message(STATUS "Applying GCC GPU Offload Flags (${GPU_ARCH})")
        target_compile_options(SimulationCore PUBLIC
            -fopenmp
            -foffload=nvptx-none
            "-foffload-options=nvptx-none=-march=${GPU_ARCH}"
            -fcf-protection=none
            -fno-stack-protector # Bypass host stack canaries leaking to the GPU
        )
        add_link_options(
            -fopenmp
            -foffload=nvptx-none
            "-foffload-options=nvptx-none=-march=${GPU_ARCH}"
            -fcf-protection=none
            -fno-stack-protector
        )
    else()
        message(FATAL_ERROR "GPU Offloading is currently only supported via nvc++ or Clang.")
    endif()
else()
    message(STATUS "GPU Offloading DISABLED (Using CPU backend)")
    
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(SimulationCore PUBLIC -foffload=disable)
        add_link_options(-foffload=disable)
    endif()
endif()
# ==============================================================================
# MAIN EXECUTABLE
# ==============================================================================
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE SimulationCore)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /wd4251)
endif()

# Copy resource files post-build
if(EXISTS ${CMAKE_SOURCE_DIR}/config.ini)
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/config.ini
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/config.ini
            ${CMAKE_CURRENT_BINARY_DIR}/config.ini
        DEPENDS ${CMAKE_SOURCE_DIR}/config.ini
        COMMENT "Copying config.ini to build directory (Source was updated)"
    )

    add_custom_target(CopyIni ALL 
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/config.ini
    )

    add_dependencies(${PROJECT_NAME} CopyIni)
endif()

# ==============================================================================
# TESTING INFRASTRUCTURE
# ==============================================================================
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.3 
)
FetchContent_MakeAvailable(Catch2)

enable_testing()

add_executable(run_tests 
    tests/test_cooling.cpp    
    tests/test_math.cpp
    tests/test_gravity.cpp
    tests/test_hydro.cpp
    tests/test_ics.cpp
    tests/test_mfm.cpp
)

target_link_libraries(run_tests PRIVATE SimulationCore Catch2::Catch2WithMain)
target_include_directories(run_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

include(CTest)
include(Catch)
catch_discover_tests(run_tests)

# ==============================================================================
# LINK-TIME OPTIMIZATION (IPO)
# ==============================================================================
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported AND ENABLE_IPO)
    message(STATUS "IPO / LTO enabled: compiler supports it.")
    set_target_properties(SimulationCore PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_target_properties(run_tests PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    if(NOT ENABLE_IPO)
        message(STATUS "IPO / LTO manually disabled via ENABLE_IPO=OFF")
    else()
        message(STATUS "IPO / LTO disabled (unsupported): ${ipo_error}")
    endif()
    
    set_target_properties(SimulationCore PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
    set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
    set_target_properties(run_tests PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)
endif()

# ==============================================================================
# CODE LINTER (CLANG-TIDY)
# ==============================================================================
if(USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
    if(CLANG_TIDY_EXE)
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set_target_properties(SimulationCore PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
        set_target_properties(${PROJECT_NAME} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    else()
        message(WARNING "clang-tidy not found!")
    endif()
endif()
