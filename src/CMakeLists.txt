include_directories(/usr/local/include/ ${OPENMP_CXX_INCLUDE_DIRS})

file( GLOB SRCS *.cpp *.hpp)
add_executable(lpm ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(BOOST 1.66 COMPONENTS regex program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(lpm GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a OpenMP::OpenMP_CXX)
