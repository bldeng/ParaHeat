cmake_minimum_required(VERSION 3.1)

project(GeodDistSolver)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Surface mesh library	
set(SURFACEMESH_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/surface_mesh")
file(GLOB_RECURSE SURFACEMESH_HEADERS "${SURFACEMESH_SOURCE_DIR}/*.h")
file(GLOB_RECURSE SURFACEMESH_SOURCES "${SURFACEMESH_SOURCE_DIR}/*.cpp")
add_library(SurfaceMesh STATIC ${SURFACEMESH_HEADERS} ${SURFACEMESH_SOURCES})

# Additional compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message("Clang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	target_compile_options(SurfaceMesh PUBLIC -pedantic -Wno-long-long)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
	message("AppleClang compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	target_compile_options(SurfaceMesh PUBLIC -pedantic -Wno-long-long)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	message("GNU compiler found.")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
	target_compile_options(SurfaceMesh PUBLIC -pedantic -Wno-long-long)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message("MSVC compiler found.")
    add_definitions(/DUSE_MSVC)
    add_definitions(/D_USE_MATH_DEFINES)
endif()


# Add the current folder into include path
include_directories(SYSTEM "${CMAKE_CURRENT_SOURCE_DIR}")

# Executable for distance solver
add_executable(GeodDistSolver
	EigenTypes.h
	OMPHelper.h
	Parameters.h
	GeodesicDistanceSolver.h
	DistanceFile.h
	GeodesicDistanceSolver.cpp
	Parameters.cpp
	ComputeDistance.cpp
)

# Executable for distance solver
add_executable(CompareDistance
	EigenTypes.h
	Parameters.h
	DistanceFile.h
	Parameters.cpp
	CompareDistance.cpp
)

# GLFW viewer
set(WITH_VIEWER ON CACHE BOOL "With Viewer")
if(WITH_VIEWER)
	add_executable(ViewScalarField
		DistanceFile.h
		ViewScalarField.cpp
	)
	
	if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/libigl/include")
		message(FATAL_ERROR "Cannot find libigl in the external folder.")
	endif()
	
	set(GLAD_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/libigl/external/glad")
	add_library(glad 
		"${GLAD_PATH}/src/glad.c"
	)
	target_include_directories(glad SYSTEM PUBLIC "${GLAD_PATH}/include")
	if(NOT WIN32)
	  target_link_libraries(glad PUBLIC ${CMAKE_DL_LIBS})
	endif()
	
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
	find_package(GLFW REQUIRED)
	target_include_directories(ViewScalarField SYSTEM PUBLIC "${GLFW_INCLUDE_DIR}")
	target_include_directories(ViewScalarField SYSTEM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/external/libigl/include")
	target_include_directories(ViewScalarField SYSTEM PUBLIC "${GLAD_PATH}/include")
	target_link_libraries(ViewScalarField ${GLFW_LIBRARIES})
	target_link_libraries(ViewScalarField glad)
endif()


# Detect Eigen3 in the external folder
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen/Eigen/Dense)
	message("Found user-provided Eigen.")
	set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen")
	target_include_directories(GeodDistSolver SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(CompareDistance SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	if(WITH_VIEWER)
		target_include_directories(ViewScalarField SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	endif()
else()
	find_package(Eigen3 REQUIRED)
	if(EIGEN3_FOUND)
		message("Found system-installed Eigen")
		target_include_directories(GeodDistSolver SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		target_include_directories(CompareDistance SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		if(WITH_VIEWER)
			target_include_directories(ViewScalarField SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
		endif()
	endif()
endif()


# Linking surface_mesh
target_link_libraries(GeodDistSolver SurfaceMesh)

# Detect OpenMP environment
set(OPENMP ON CACHE BOOL "OpenMP")
if(OPENMP)
  FIND_PACKAGE(OpenMP QUIET)
  if(OPENMP_FOUND)
      message("OpenMP found. OpenMP activated in release.")
      target_compile_options(GeodDistSolver PUBLIC "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
      target_compile_definitions(GeodDistSolver PUBLIC "$<$<CONFIG:RELEASE>:USE_OPENMP>")
      target_link_libraries(GeodDistSolver "$<$<CONFIG:RELEASE>:${OpenMP_CXX_FLAGS}>")
  else()
      message("OpenMP not found.")
  endif()
endif()

