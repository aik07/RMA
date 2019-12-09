cmake_minimum_required(VERSION 3.5)
project(RMA CXX)
include(CMakePackageConfigHelpers)

set(include_install_dir "include" CACHE PATH "Relative location where to install headers.")
set(rma_header_install_dir ${include_install_dir}/pebbl)
set(lib_install_dir "lib" CACHE PATH "Relative location where to install libraries.")
set(package_install_dir ${lib_install_dir}/rma)

message("MPI enabled.")
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})


add_subdirectory(src)

configure_package_config_file(rma-config.cmake.in rma-config.cmake
                              INSTALL_DESTINATION ${package_install_dir}
                              PATH_VARS include_install_dir)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/rma-config.cmake DESTINATION ${package_install_dir})
