cmake_minimum_required(VERSION 3.15)
set(TARGET_NAME pasr)
project(${TARGET_NAME})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

set(HOME /home/ihme/mbonanni)
# set(HOME /Users/mbonanni)
set(UTILITIES_ROOT ${HOME}/local)

find_package(OpenMP)

# set(CANTERA_ROOT /opt/homebrew/Caskroom/miniforge/base/envs/cantera_dev)
set(CANTERA_ROOT ${UTILITIES_ROOT}/cantera)
set(SUNDIALS_ROOT ${UTILITIES_ROOT}/sundials/3.1.0/)
set(EIGEN3_ROOT ${UTILITIES_ROOT}/eigen/3.2.9/eigen)
set(BOOST_ROOT ${HOME}/tools/boost_1_55_0/installationDir)

include_directories(
    ${CANTERA_ROOT}/include
    ${SUNDIALS_ROOT}/include
    ${EIGEN3_ROOT}/include/eigen
    ${BOOST_ROOT}/include)
# include_directories(
#     ${CANTERA_ROOT}/include)
link_directories(
    ${CANTERA_ROOT}/build/lib
    ${SUNDIALS_ROOT}/lib
    ${})
# link_directories(
#     ${CANTERA_ROOT}/lib
#     ${})

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
        sundials_cvodes
        sundials_ida
        sundials_nvecserial
        OpenMP::OpenMP_CXX)
endif()
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(
#         ${TARGET_NAME}
#         PUBLIC
#         cantera
#         OpenMP::OpenMP_CXX)
# endif()