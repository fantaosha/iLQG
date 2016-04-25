cmake_minimum_required(VERSION 2.8)

find_package(Eigen3 REQUIRED)
find_package(SNOPT REQUIRED)
find_package(Matlab REQUIRED)
#find_package(OpenMP REQUIRED)

#if (OPENMP_FOUND)
#	message("OPENMP FOUND")
#	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

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

#add_executable(quadrotor_ddp_est quadrotor_ddp_est.cpp)
#target_link_libraries(quadrotor_ddp_est ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})

add_executable(quadrotor_ddp quadrotor_ddp.cpp)
target_link_libraries(quadrotor_ddp ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})

#add_executable(quadrotor_lqr quadrotor_lqr.cpp)
#target_link_libraries(quadrotor_lqr ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})

#add_executable(car_ddp car_ddp.cpp)
#target_link_libraries(car_ddp)

#add_executable(test test.cpp)
#target_link_libraries(test ${SNOPT_LIBRARIES} ${MATLAB_LIBRARIES})
