cmake_minimum_required(VERSION 3.24)

project(ArekoRemesher VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Floating point rules that keep results reproducible: no contraction to FMA,
# no fast-math, no constant folding of libm calls, baseline x86-64 only
# (no -march). Every core target links this one interface target.
add_library(arm_rechenregeln INTERFACE)
target_compile_options(arm_rechenregeln INTERFACE
    -ffp-contract=off
    -fno-fast-math
    -fno-builtin
    -Wall -Wextra -Wpedantic)

add_subdirectory(kern)
