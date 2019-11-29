cmake_minimum_required(VERSION 2.8.12)

find_package(pybind11 REQUIRED)

set(Mesh_DIR "${PROJECT_BINARY_DIR}/cpp/mesh")
find_package(Mesh REQUIRED)

add_library(pymesh MODULE pymesh.cc)
target_compile_features(pymesh PRIVATE cxx_std_17)

target_link_libraries(pymesh PRIVATE pybind11::module)
target_link_libraries(pymesh PRIVATE Mesh::mesh)
set_target_properties(pymesh PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
#        SUFFIX "${PYTHON_MODULE_EXTENSION}")