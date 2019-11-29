cmake_minimum_required(VERSION 3.13)
project(pyulb)

set(CMAKE_CXX_STANDARD 14)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

option(USE_OMP "use OpenMP" ON)
option(USE_MPI "use MPI" ON)
option(BUILD_TEST "build tests" OFF)

#find_package(MKL)
#if(MKL_FOUND)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DMKL_LP64")
#    target_include_directories(pyulb PRIVATE ${MKL_INCLUDE_DIRS})
#    target_link_libraries(pyulb ${MKL_LIBRARIES} pthread dl)
#else()
#    find_package(BLAS REQUIRED)
#    if (BLAS_FOUND)
#        message("MKL is not found using blas instead")
#        target_link_libraries(pyulb blas)
#    else()
#        message(STATUS "No blas library could be found.")
#    endif()
#endif()
#
#if(USE_MPI)
#    find_package(MPI)
#    if(MPI_CXX_FOUND)
#        target_include_directories(pyulb PRIVATE ${MPI_CXX_INCLUDE_PATH})
#        set_target_properties(pyulb PROPERTIES COMPILE_FLAGS "-DMPI_ENABLED")
#        if(MPI_CXX_COMPILE_FLAGS)
#            set_target_properties(pyulb PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
#        endif()
#        if (MPI_CXX_LINK_FLAGS)
#            set_target_properties(pyulb PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
#        endif()
#        if (MPI_CXX_COMPILER)
#            set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#        endif()
#        message(STATUS "Configured MPI compilation.")
#    else()
#        message(STATUS "MPI could not be found.")
#    endif()
#    target_include_directories(pyulb PRIVATE ${MPI_CXX_INCLUDE_PATH})
#    target_link_libraries(pyulb ${MPI_CXX_LIBRARIES})
#endif()

if (USE_OMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(cpp/mesh)
add_subdirectory(py/mesh)

#find_package(PythonInterp REQUIRED)
#add_subdirectory(swig)



if(BUILD_TEST)
    add_subdirectory(test)
endif()


#include_directories(lib)

#install(TARGETS pyulb LIBRARY DESTINATION pyulb)
