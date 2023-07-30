# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.27)

project(POSSUM LANGUAGES C CXX)


include(FetchContent)

# add_subdirectory(C++/Lanczos)
add_subdirectory(C++/PowerMethod)
