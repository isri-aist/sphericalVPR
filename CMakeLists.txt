# ############################################################################
#
# This file is part of the Spherical VPR software.
# Copyright (C) 2024 by CNRS-AIST JRL. All rights reserved.
#
# This software was developed at:
# CNRS - AIST JRL (Joint Robotics Laboratory)
# 1-1-1 Umezono, Tsukuba, Ibaraki
# Japan
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# Short script to compile the main.cpp file with the OpenSeqSLAM class and the
# OpenCV library. The script also uses the Boost library to handle the file
# system and the Eigen library to handle the matrix operations. The script
# also uses the libPeR library to handle the data acquisition from the camera and
# project the images on the sphere (with various representations available).
#
# Authors:
# Antoine ANDRÉ, Guillaume CARON
#
# ############################################################################

project(main)

cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(PER REQUIRED per_core per_io per_acquisition_representation)

if(PER_FOUND)
  include(${PER_USE_FILE})
endif(PER_FOUND)

find_package(OpenCV REQUIRED)

find_package(Eigen3 REQUIRED)

# Boost
FIND_PACKAGE(Boost REQUIRED)

# find_package(Fitsio REQUIRED)
set(header_files "./header/OpenSeqSLAM.h" "./header/ssmlSeqSLAM.hpp")
set(src_files "./src/OpenSeqSLAM.cpp" "./main.cpp")
add_executable(main ${header_files} ${my_cpps})
include_directories("./header")

include_directories(${Boost_INCLUDE_DIRS}) # /opt/local/include/ might be needed as well under MacOS
link_directories(${Boost_LIBRARY_DIRS}) # /opt/local/lib/ might be needed as well under MacOS # similar /Users/guillaume/Development/libraries/visp-3.0.1/build/lib/Release/ might be needed as well under MacOS

foreach(cpp ${main_cpp})
  per_add_target(${cpp})

  if(COMMAND pr_add_dependency)
    pr_add_dependency(${cpp} "main")
  endif()
endforeach()

target_link_libraries(main libboost_system.so libboost_filesystem.so libboost_regex.so)
target_link_libraries(main ${OpenCV_LIBS})
