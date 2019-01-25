include_directories(/usr/local/include/)

file( GLOB SRCS *.cpp *.hpp)
add_executable(lpm ${SRCS}) # create the executable file

find_package (PkgConfig REQUIRED)
find_package(BOOST 1.66 COMPONENTS regex  REQUIRED)
find_package (GSL REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

target_link_libraries(lpm GSL::gsl GSL::gslcblas /usr/local/bin/libSPF.a /usr/local/lib/libboost_system.a /usr/local/lib/libboost_filesystem.a)
