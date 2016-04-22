cmake_minimum_required(VERSION 2.8)

find_package(Eigen3 REQUIRED)
find_package(SNOPT REQUIRED)
find_package(Matlab REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${MATLAB_INCLUDE_DIR})
include_directories(${SNOPT_DIRS})

include_directories(headers)

if (CMAKE_BUILD_TYPE)
	message("The build type is " ${CMAKE_BUILD_TYPE})
	if (CMAKE_BUILD_TYPE  STREQUAL Debug)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fkeep-inline-functions")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
	endif()
endif()

message(${CMAKE_CXX_FLAGS})

#add_executable(test2 test2.cpp)
#target_link_libraries(test2 ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})

add_executable(quadrotor quadrotor.cpp)
target_link_libraries(quadrotor ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})

#add_executable(car car.cpp)
#target_link_libraries(car ${SNOPT_LIBRARIES})

#add_executable(test test.cpp)
#target_link_libraries(test ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})
