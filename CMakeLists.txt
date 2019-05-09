project(lpm)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_C_FLAGS "-Wall -g" )
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")

if (APPLE)
    set(OPENMP_CXX_INCLUDE_DIRS /usr/local/opt/libomp/include)
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    message(STATUS "hey")
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
