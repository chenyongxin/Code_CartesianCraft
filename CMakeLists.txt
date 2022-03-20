cmake_minimum_required(VERSION 3.12)

SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
# Toggleable build type
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# Find the modules included with CCC by appending the `cmake' folder
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
 
project(Code_CartesianCraft LANGUAGES C CXX)

# CCC requires C++11
SET(CMAKE_CXX_STANDARD 11 CACHE STRING "Default value for CXX_STANDARD property of targets")
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Default C/CXX compilers
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Dependencies
INCLUDE(FindHYPRE)
INCLUDE(FindHDF5)

if(NOT HYPRE_FOUND)
	message(FATAL_ERROR "No HYPRE found!")
endif()

# Compile
set(EXE CCC)
file(GLOB ccfiles "${PROJECT_SOURCE_DIR}/src/*.cc")
add_executable(${EXE} ${ccfiles})
target_include_directories(${EXE} PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(${EXE} PRIVATE ${HYPRE_INCLUDE_DIR})
target_link_libraries(${EXE} PRIVATE ${HYPRE_LIBRARY})
if(HDF5_FOUND)
	target_compile_definitions(${EXE} PRIVATE USE_HDF)
	target_include_directories(${EXE} PRIVATE ${HDF5_INCLUDE_DIR})
	target_link_libraries(${EXE} PRIVATE ${HDF5_LIBRARY})
endif()

# Install
include(GNUInstallDirs)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH "" FORCE)
endif()

install(TARGETS ${EXE}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})