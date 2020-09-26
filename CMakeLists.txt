cmake_minimum_required(VERSION 3.0)
project(rma CXX)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(PEBBL_ROOT ${CMAKE_SOURCE_DIR}/external/pebbl)

#include(ExternalProject)
#set( EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external_install )
#add_subdirectory(external)

set(PEBBL_LIBRARIES ${PEBBL_ROOT}/build/src/pebbl/libpebbl.a)
set(PEBBL_INCLUDES  ${PEBBL_ROOT}/src ${PEBBL_ROOT}/build/src)

#set(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
  # set(CMAKE_BUILD_TYPE "Debug")
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall -Werror=reorder")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

MESSAGE( STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}" )
MESSAGE( STATUS "C++ compiler ID: ${CMAKE_CXX_COMPILER_ID}" )
MESSAGE( STATUS "C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}" )
MESSAGE( STATUS "C++ compiler string: ${CXX_COMPILER_STRING}" )
MESSAGE( STATUS "MPI C++ compiler: ${MPI_CXX_COMPILER}" )
MESSAGE( STATUS "MPI C++ compiler options: ${MPI_CXX_COMPILE_OPTIONS}" )
MESSAGE( STATUS "MPI libraries: ${MPI_LIBRARIES}")
MESSAGE( STATUS "PEBBL include dir: ${PEBBL_INCLUDES}" )
MESSAGE( STATUS "PEBBL libraries: ${PEBBL_LIBRARIES}")

include_directories(
  src
  ${MPI_INCLUDE_PATH}
  ${CMAKE_SOURCE_DIR}
	${PEBBL_INCLUDES}
)

file(GLOB SOURCES "src/*.cpp")
add_executable(rma ${SOURCES})
target_link_libraries(rma ${PEBBL_LIBRARIES})

#file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
add_library(lrma STATIC ${SOURCES})
target_link_libraries(lrma ${PEBBL_LIBRARIES})
