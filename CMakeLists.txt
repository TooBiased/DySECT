cmake_minimum_required (VERSION 2.8)

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

project (space-grow-table)

set (CMAKE_CXX_FLAGS "-std=c++14 -msse4.2 -Wall -Wextra -O3")


set(GROWT_HASHFCT XXHASH CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE GROWT_HASHFCT PROPERTY STRINGS XXHASH MURMUR2 MURMUR3 CRC)


include_directories(.)


if (GROWT_HASHFCT STREQUAL XXHASH)
  find_package(XXHash)
  if (XXHASH_FOUND)
    include_directories(${XXHASH_INCLUDE_DIRS})
  else()
    message("Cannot find xxHash. Use MURMUR2 instead.")
    set(GROWT_HASHFCT MURMUR2)
  endif()
endif()

if (GROWT_HASHFCT STREQUAL MURMUR3)
  find_package(smhasher)
  if (SMHASHER_FOUND)
    include_directories(${SMHASHER_INCLUDE_DIRS})
  else()
    message("Cannot find SMHasher which implements MURMUR3. Use MURMUR2 instead.")
  endif()
endif()



add_executable(in source/in_test.cpp)
