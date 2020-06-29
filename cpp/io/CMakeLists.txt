add_library(IO SHARED logger.cc)
add_library(Mesh::IO ALIAS IO)

#target_compile_definitions(IO INTERFACE IO=1)
set(IO_Version 1.0)

include(GNUInstallDirs)
install(TARGETS IO
        EXPORT IOTargets
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT IOTargets
        FILE IOTargets.cmake
        NAMESPACE Mesh::
        DESTINATION lib/cmake/Mesh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("IOConfigVersion.cmake"
        VERSION ${IO_Version}
        COMPATIBILITY SameMajorVersion)
install(FILES "cmake/IOConfig.cmake" "IOConfigVersion.cmake"
        DESTINATION lib/cmake/mesh)

export(TARGETS IO FILE cmake/IOConfig.cmake NAMESPACE Mesh::)