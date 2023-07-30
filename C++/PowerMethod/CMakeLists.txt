# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(ZLIB REQUIRED)
find_package(CURL REQUIRED)

add_executable(POSSUM_power)
target_sources(POSSUM_power PRIVATE GWevIntra_new.cpp hgFlipSign.cpp thdMul.cpp theBestEigen.cpp "${straw_SOURCE_DIR}/C++/straw.cpp")
target_include_directories(POSSUM_power PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} "${straw_SOURCE_DIR}/C++/")
target_link_libraries(POSSUM_power PUBLIC CURL::libcurl ZLIB::ZLIB)

# Workaround missing #include <cstdint> in straw
if(MSVC)
  add_definitions(/FI"cstdint")
else()
  add_definitions(-include cstdint)
endif()

