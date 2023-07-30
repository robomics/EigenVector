# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.27)

project(POSSUM LANGUAGES C CXX)


include(FetchContent)
FetchContent_Declare(
    straw
    GIT_REPOSITORY https://github.com/aidenlab/straw.git
    GIT_TAG        2525edc29bbb48463799cad94cbd6e5e810210a0
)

FetchContent_MakeAvailable(straw)

# add_subdirectory(C++/Lanczos)
add_subdirectory(C++/PowerMethod)
