cmake_minimum_required(VERSION 3.15)
set(TARGET_NAME pasr)
project(${TARGET_NAME})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CANTERA_ROOT "" CACHE PATH "Path to Cantera installation directory")
if(NOT CANTERA_ROOT)
    message(FATAL_ERROR "Please provide the path to the Cantera installation directory.")
endif()

find_package(OpenMP)
find_package(BLAS REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(${CANTERA_ROOT}/include)

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

target_include_directories(${TARGET_NAME} PUBLIC ${CANTERA_ROOT}/include)
target_link_directories(${TARGET_NAME} PUBLIC ${CANTERA_ROOT}/build/lib)
target_link_libraries(${TARGET_NAME}
    PUBLIC
    cantera
    ${CMAKE_DL_LIBS}
    ${BLAS_LIBRARIES}
    yaml-cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

install(TARGETS ${TARGET_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})