cmake_minimum_required(VERSION 3.15)
project(pasr)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(HOME /home/ihme/mbonanni)
set(UTILITIES_ROOT ${HOME}/local)

find_package(OpenMP)

set(CANTERA_ROOT ${UTILITIES_ROOT}/cantera)
set(SUNDIALS_ROOT ${UTILITIES_ROOT}/sundials/3.1.0/)
set(EIGEN3_ROOT ${UTILITIES_ROOT}/eigen/3.2.9/eigen)
set(BOOST_ROOT ${HOME}/tools/boost_1_55_0/installationDir)

include_directories(
    ${CANTERA_ROOT}/include
    ${SUNDIALS_ROOT}/include
    ${EIGEN3_ROOT}/include/eigen
    ${BOOST_ROOT}/include)
link_directories(
    ${CANTERA_ROOT}/build/lib
    ${SUNDIALS_ROOT}/lib
    ${})

add_executable(pasr
        src/pasr.cpp
        src/solver.cpp
        src/solver.h)
target_link_libraries(pasr)
if(OpenMP_CXX_FOUND)
    target_link_libraries(pasr PUBLIC cantera OpenMP::OpenMP_CXX)
endif()
