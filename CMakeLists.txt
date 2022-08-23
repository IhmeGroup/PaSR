cmake_minimum_required(VERSION 3.15)
project(pasr)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(HOME /Users/mbonanni)
set(UTILITIES_ROOT ${HOME}/local)

find_package(OpenMP)

set(CANTERA_ROOT ${UTILITIES_ROOT}/cantera)

include_directories(${CANTERA_ROOT}/include)
link_directories(${CANTERA_ROOT}/build/lib ${})

add_executable(pasr
        src/pasr.cpp
        src/solver.cpp
        src/solver.h)
target_link_libraries(pasr)
if(OpenMP_CXX_FOUND)
    target_link_libraries(pasr PUBLIC cantera OpenMP::OpenMP_CXX)
endif()
