# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


add_executable(POSSUM_power)
target_sources(POSSUM_power PRIVATE GWevIntra_new.cpp hgFlipSign.cpp thdMul.cpp theBestEigen.cpp)
target_link_libraries(POSSUM_power PUBLIC hictk::file)
target_include_directories(POSSUM_power PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(POSSUM_power INTERFACE "cxx_std_${CMAKE_CXX_STANDARD}")
