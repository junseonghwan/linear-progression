include_directories(/usr/local/include/ ${OPENMP_CXX_INCLUDE_DIRS} /home/seongjun/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp)

#set (Boost_USE_STATIC_LIBS OFF)
set (Boost_DEBUG ON)
#add_definitions(-DBOOST_ALL_NO_LIB)
#add_definitions(-BOOST_ALL_DYN_LINK)

# create executable file
#add_executable(lpm ${SRCS})
# add tests
add_executable(tests ${SRCS})
add_test(Tester tests)

# build as library
add_library(lpm_lib SHARED ${SRCS})


find_package (PkgConfig REQUIRED)
find_package (GSL REQUIRED)
find_package(Boost 1.67 COMPONENTS regex program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_libraries(${Boost_LIBRARIES})
endif()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message(STATUS "OpenMP found!")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(lpm_lib GSL::gsl GSL::gslcblas /home/seongjun/usr/local/bin/libSPF.a /home/seongjun/boost_1_67_0/stage/lib/libboost_system.a /home/seongjun/boost_1_67_0/stage/lib/libboost_filesystem.a /home/seongjun/boost_1_67_0/stage/lib/libboost_program_options.a OpenMP::OpenMP_CXX)
target_link_libraries(tests GSL::gsl GSL::gslcblas /home/seongjun/usr/local/bin/libSPF.a /home/seongjun/boost_1_67_0/stage/lib/libboost_system.a /home/seongjun/boost_1_67_0/stage/lib/libboost_filesystem.a /home/seongjun/boost_1_67_0/stage/lib/libboost_program_options.a OpenMP::OpenMP_CXX)

message(${CMAKE_SOURCE_DIR}/bin)

# somehow WORKING_DIRECTORY command does not work!
#add_test(Tester tests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
INSTALL(TARGETS lpm_lib LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/bin/)

