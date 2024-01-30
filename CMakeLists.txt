# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.25)
cmake_policy(VERSION 3.25...3.27)

project(POSSUM LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

include(FetchContent)

FetchContent_Declare(
        hictk
        GIT_REPOSITORY https://github.com/paulsengroup/hictk.git
        GIT_TAG v0.0.8
        EXCLUDE_FROM_ALL SYSTEM)

set(HICTK_ENABLE_GIT_VERSION_TRACKING OFF)
set(HICTK_BUILD_TOOLS OFF)
set(HICTK_ENABLE_TESTING OFF)
FetchContent_MakeAvailable(hictk)

# add_subdirectory(C++/Lanczos)
add_subdirectory(C++/PowerMethod)
