cmake_minimum_required(VERSION 3.24)

project(ArekoRemesher VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Floating point rules that keep results reproducible: no contraction to FMA,
# no fast-math, no constant folding of libm calls, baseline instruction set
# only (no -march, no /arch:AVX2). Every core target links this one interface target.
add_library(arm_rechenregeln INTERFACE)
if(MSVC)
    # Static runtime: the DLL needs no runtime DLLs besides the system ones.
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    target_compile_options(arm_rechenregeln INTERFACE
        /fp:precise
        /bigobj)
else()
    target_compile_options(arm_rechenregeln INTERFACE
        -ffp-contract=off
        -fno-fast-math
        -fno-builtin
        -Wall -Wextra -Wpedantic)
endif()

add_subdirectory(kern)
