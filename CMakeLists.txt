cmake_minimum_required(VERSION 3.13)
project(hard-sphere-simulation CXX)

cmake_policy(SET CMP0076 NEW)

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_library(hardsphere)
target_compile_features(hardsphere
  PUBLIC
    cxx_std_17
  )
target_compile_options(hardsphere
  PUBLIC
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<CXX_COMPILER_ID:Intel>:$<IF:$<PLATFORM_ID:Windows>,/warn:all,-warn all>>
  )

target_link_libraries(hardsphere
  PUBLIC
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
  )
target_include_directories(hardsphere
  PUBLIC
    include
  )

add_subdirectory(include/shirose)
add_subdirectory(src)