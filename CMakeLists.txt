project(lpm)

cmake_minimum_required (VERSION 3.12)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_C_FLAGS "-Wall -g" )

set(CMAKE_LIBRARY_PATH "/home/seongjun/boost_1_67_0/stage/lib/")
set(CMAKE_INCLUDE_PATH "/home/seongjun/boost_1_67_0/boost/")
set(BOOST_ROOT "/home/seongjun/boost_1_67_0/")
set(BOOST_INCLUDEDIR "/home/seongjun/boost_1_67_0/boost")
set(BOOST_LIBRARYDIR "/home/seongjun/boost_1_67_0/stage/lib")
#set(Boost_DEBUG ON)

set(OPENMP_CXX_INCLUDE_DIRS /usr/local/opt/libomp/include)
if (APPLE AND XCODE)
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
endif()

enable_testing()

add_subdirectory(src)
