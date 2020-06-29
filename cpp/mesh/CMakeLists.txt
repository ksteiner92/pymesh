add_library(Mesh SHARED mesh.cc segment.cc elements.cc system.cc meshing.cc)
add_library(Mesh::Mesh ALIAS Mesh)

target_compile_features(Mesh PRIVATE cxx_std_17)

find_package(Eigen3 REQUIRED)

add_subdirectory(tetgen)
add_subdirectory(triangle)


target_include_directories(Mesh
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        ${EIGEN3_INCLUDE_DIR}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
        $<INSTALL_INTERFACE:../include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/cpp/include
        ${CMAKE_SOURCE_DIR}/cpp/io/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../../lib
)

# If we have compiler requirements for this library, list them
# here
#target_compile_features(lib
#        PUBLIC cxx_auto_type
#        PRIVATE cxx_variadic_templates)

target_link_libraries(mesh PRIVATE tetgen triangle)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS Mesh
        EXPORT MeshConfig
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(EXPORT MeshConfig NAMESPACE Mesh:: DESTINATION share/cmake/mesh)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

# This makes the project importable from the build directory
export(TARGETS mesh FILE MeshConfig.cmake NAMESPACE Mesh::)
