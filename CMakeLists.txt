cmake_minimum_required(VERSION 2.6)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()
project(StepRecoveryMPC CXX)

find_package(YCM REQUIRED)
include(YCMDefaultDirs)

find_package(iDynTree REQUIRED)
find_package(IPOPT REQUIRED)

find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${IPOPT_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

set(SRMPC_SOURCES src/MPCIpOptSolver.cpp)
set(SRMPC_HEADERS include/MPCIpOptSolver.h)

add_library(${PROJECT_NAME}  SHARED ${SRMPC_HEADERS} ${SRMPC_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} ${iDynTree_LIBRARIES} ${IPOPT_LIBRARIES})
