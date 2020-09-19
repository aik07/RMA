#########################################################
message( STATUS "Building external PEBBL project.")
#########################################################

set(PEBBL_ROOT ${CMAKE_SOURCE_DIR}/external/pebbl)

ExternalProject_Add(pebbl_external

                    # PREFIX ${PEBBL_ROOT}

                    GIT_REPOSITORY "https://github.com/PEBBL/pebbl.git"
                    GIT_TAG  "master"

                    UPDATE_COMMAND ""
                    PATCH_COMMAND ""

                    # BINARY_DIR ${PEBBL_ROOT}
                    SOURCE_DIR ${PEBBL_ROOT}
                    # INSTALL_DIR ${PEBBL_ROOT}

                    # CMAKE_COMMAND cmake .

                    CMAKE_ARGS -Denable_mpi=ON -Denable_examples=OFF
                               -DCMAKE_INSTALL_PREFIX=${PEBBL_ROOT}
                               PEBBL_ROOT=${PEBBL_ROOT}

                    # BUILD_COMMAND make
                    # BUILD_BYPRODUCTS ${PEBBL_LIBRARIES}

                    TEST_COMMAND ""
                    )


ExternalProject_Add_Step(pebbl_external
                         bootstrap
                         COMMAND make
                         DEPENDEES download
                         DEPENDERS configure
                         WORKING_DIRECTORY ${PEBBL_ROOT})


# ExternalProject_Get_Property(pebbl_external source_dir)

message(${CMAKE_BINARY_DIR})

set(PEBBL_INCLUDES ${PEBBL_ROOT}/src) # ${PEBBL_ROOT}/build/src)
set(PEBBL_LIBRARIES ${PEBBL_ROOT}/src/pebbl/libpebbl.a)

add_library(pebbl STATIC IMPORTED)
set_target_properties(pebbl PROPERTIES IMPORTED_LOCATION ${PEBBL_ROOT}/build/src/pebbl/libpebbl.a)
add_dependencies(pebbl pebbl_external)
