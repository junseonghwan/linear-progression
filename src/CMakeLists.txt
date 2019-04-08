include_directories(/usr/local/include/ ${OPENMP_CXX_INCLUDE_DIRS})

file( GLOB SRCS *.cpp *.hpp)

# create executable file
#add_executable(lpm ${SRCS})
# add tests
add_executable(tests ${SRCS})
add_test(Tester tests)

# build as library
add_library(lpm_lib SHARED ${SRCS})

find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
#find_package(Boost 1.67 COMPONENTS regex program_options filesytem REQUIRED)
find_package(Boost 1.67 COMPONENTS regex program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#target_link_libraries(lpm GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a OpenMP::OpenMP_CXX)
target_link_libraries(tests GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a OpenMP::OpenMP_CXX)
target_link_libraries(lpm_lib GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a /usr/local/lib/libboost_program_options.a OpenMP::OpenMP_CXX)

message(${CMAKE_SOURCE_DIR}/bin)
#file(COPY ${CMAKE_SOURCE_DIR}/build/src/liblpm_lib.dylib DESTINATION ${CMAKE_SOURCE_DIR}/bin/)

# somehow WORKING_DIRECTORY command does not work!
#add_test(Tester tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
INSTALL(TARGETS lpm_lib LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/bin/)

