cmake_minimum_required(VERSION 3.15)
set(TARGET_NAME pasr)
project(${TARGET_NAME})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(HOME /Users/mbonanni)
set(UTILITIES_ROOT ${HOME}/local)

find_package(OpenMP)

set(CANTERA_ROOT /opt/homebrew/Caskroom/miniforge/base/envs/cantera_dev)

include_directories(
    ${CANTERA_ROOT}/include)

link_directories(
    ${CANTERA_ROOT}/lib
    ${})

add_executable(${TARGET_NAME}
        src/pasr.cpp
        src/common.h
        src/PartiallyStirredReactor.h
        src/PartiallyStirredReactor.cpp
        src/Particle.h
        src/Particle.cpp
        src/Injector.h
        src/Injector.cpp
        src/Histogram.h
        src/Histogram.cpp)

target_link_libraries(${TARGET_NAME})

if(OpenMP_CXX_FOUND)
    target_link_libraries(
        ${TARGET_NAME}
        PUBLIC
        cantera
        OpenMP::OpenMP_CXX)
endif()