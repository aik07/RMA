cmake_minimum_required(VERSION 3.0)
project(rma)

include_directories(src)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_PREFIX_PATH "../../pebbl/build")

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall" )

set ( PEBBBL_LINK_LIBS libpebbl.a )
link_directories( ../../pebbl/build/src/pebbl  )

FIND_PACKAGE( MPI   REQUIRED )

file(GLOB SOURCES "src/*.cpp")

add_executable(rma ${SOURCES})

add_library( rma STATIC ${SOURCES} )

install( TARGET rma DESTINATION PEBBL_LINK_LIBS )

#target_link_libraries(librma ${PROJECT_LINK_LIBS} )
