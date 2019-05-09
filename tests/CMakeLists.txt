include_directories(/usr/local/include/ ${OPENMP_CXX_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/src)

file( GLOB SRCS *.cpp *.hpp)

# add tests
add_executable(tests ${SRCS})
add_test(Tester tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(Boost 1.67 COMPONENTS regex program_options REQUIRED)
find_package(OpenMP REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_libraries(${Boost_LIBRARIES})
endif()
if (OPENMP_FOUND)
    message(STATUS "OpenMP found!")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(tests GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a OpenMP::OpenMP_CXX)

INSTALL(TARGETS tests DESTINATION ${CMAKE_SOURCE_DIR}/bin/)

