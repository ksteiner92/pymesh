add_library(Math INTERFACE logger.cc)
add_library(Mesh::Math ALIAS Math)

target_compile_features(Math PRIVATE cxx_std_17)
target_compile_definitions(Math INTERFACE Math=1)
target_include_directories(Math
        PUBLIC
        $<INSTALL_INTERFACE:include/expression.hh>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/expression.hh>
        )

set(Math_Version 1.0)

include(GNUInstallDirs)
install(TARGETS Math
        EXPORT MathTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT MathTargets
        FILE MathTargets.cmake
        NAMESPACE Mesh::
        DESTINATION lib/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("MathConfigVersion.cmake"
        VERSION ${Math_Version}
        COMPATIBILITY SameMajorVersion)
install(FILES "cmake/MathConfig.cmake" "MathConfigVersion.cmake"
        DESTINATION lib/cmake/mesh)

export(TARGETS Math FILE cmake/MathConfig.cmake NAMESPACE Mesh::)