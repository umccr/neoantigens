cmake_minimum_required(VERSION 2.8)

project(readsSelectForHLAminer_cpp)

if(NOT CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr")
endif()
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/vendor;${CMAKE_SOURCE_DIR}/build-common/cmake")

set(EXE_VERSION_SUFFIX ${FULL_VERSION})

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)


#if (CMAKE_COMPILER_IS_GNUCC)
#    set(CMAKE_C_FLAGS "-std=c++0x -Wall -ansi -pedantic")
#endif ()

#if (CMAKE_COMPILER_IS_GNUCXX)
#    set(CMAKE_CXX_FLAGS "-std=gnu++0x -Wall -ansi -pedantic -Wno-long-long -Wno-variadic-macros")
#endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR CMAKE_COMPILER_IS_GNUCXX)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(--std=c++11 SUPPORTS_STD_CXX11)
    check_cxx_compiler_flag(--std=c++0x SUPPORTS_STD_CXX01)
    if(SUPPORTS_STD_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -Wall -pthread")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --std=c++11 -Wall -pthread")
    elseif(SUPPORTS_STD_CXX01)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++0x -Wall -pthread")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --std=c++0x -Wall -pthread")
    else()
        message(ERROR "Compiler does not support --std=c++11 or --std=c++0x.")
    endif()
endif()


# make sure to pick up headers from library dirs
include_directories("src")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
enable_testing(true)


add_subdirectory(src)

