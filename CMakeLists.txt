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

find_package(iDynTree 0.7.2 REQUIRED)
find_package(IPOPT REQUIRED)
find_package(YARP REQUIRED)

find_package(Eigen3 REQUIRED)

include_directories(SYSTEM ${IPOPT_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
include_directories(${YARP_INCLUDE_DIRS})

set(MPCIPOPT_SOURCES src/MPCIpOptSolver.cpp)
set(MPCIPOPT_HEADERS include/MPCIpOptSolver.h)
add_library(MPCIpOpt  SHARED ${MPCIPOPT_HEADERS} ${MPCIPOPT_SOURCES})
target_include_directories(MPCIpOpt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(MPCIpOpt ${iDynTree_LIBRARIES} ${IPOPT_LIBRARIES})

set(STEPRECOVERY_SOURCES src/StepRecoveryMPC.cpp)
set(STEPRECOVERY_HEADERS include/StepRecoveryMPC.h)
add_library(srMPC  SHARED ${STEPRECOVERY_HEADERS} ${STEPRECOVERY_SOURCES})
target_include_directories(srMPC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(srMPC ${iDynTree_LIBRARIES} ${IPOPT_LIBRARIES} MPCIpOpt)


set(COORDINATOR_SOURCES src/StepRecoveryCoordinator.cpp)
set(COORDINATOR_HEADERS include/StepRecoveryCoordinator.h)
add_executable(${PROJECT_NAME} ${COORDINATOR_SOURCES} ${COORDINATOR_HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} ${iDynTree_LIBRARIES} ${YARP_LIBRARIES} srMPC)
