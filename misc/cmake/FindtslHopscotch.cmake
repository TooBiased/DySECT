################################################################################
# misc/cmake/FindXXHash.cmake
#
# Finds xxHash directory
# Looks in PATH and XXHASH_ROOT environment variables
#
# Part of Project growt - https://github.com/TooBiased/growt.git
#
# Copyright (C) 2015-2016 Tobias Maier <t.maier@kit.edu>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

message(STATUS "Looking for xxHash Implementation")

find_path(TSLHOPSCOTCH_PATH hopscotch-map/src/hopscotch_map.h
  PATHS ENV PATH ENV TSLHOPSCOTCH_ROOT)

if(TSLHOPSCOTCH_PATH)
  message(STATUS "Looking for tslHopscotch - found path: ${TSLHOPSCOTCH_PATH}/hopscotch-map")
  set(TSLHOPSCOTCH_FOUND ON)

  set(TSLHOPSCOTCH_INCLUDE_DIRS "${TSLHOPSCOTCH_PATH}/hopscotch-map/src/")
else()
  message(STATUS "failed finding tslHopscotch - try setting the TSL_HOPSCOTCH_ROOT environment variable or append it to PATH")
  if (TSLHOPSCOTCH_FIND_REQUIRED)
    message(FATAL_ERROR "Required package tslHopscotchg missing!")
  endif()
endif()
