cmake_minimum_required (VERSION 2.8)

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

project (space-grow-table)

set (CMAKE_CXX_FLAGS "-std=c++14 -msse4.2 -Wall -Wextra -O3 -g")


set(GROWS_HASHFCT XXHASH CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE GROWS_HASHFCT PROPERTY STRINGS XXHASH MURMUR2 MURMUR3 CRC)


include_directories(.)


if (GROWS_HASHFCT STREQUAL XXHASH)
  find_package(XXHash)
  if (XXHASH_FOUND)
    include_directories(${XXHASH_INCLUDE_DIRS})
  else()
    message("Cannot find xxHash. Use MURMUR2 instead.")
    set(GROWS_HASHFCT MURMUR2)
  endif()
endif()

if (GROWS_HASHFCT STREQUAL MURMUR3)
  find_package(smhasher)
  if (SMHASHER_FOUND)
    include_directories(${SMHASHER_INCLUDE_DIRS})
  else()
    message("Cannot find SMHasher which implements MURMUR3. Use MURMUR2 instead.")
    set(GROWS_HASHFCT MURMUR2)
  endif()
endif()



add_executable(in_grows   source/in_test.cpp)
target_compile_definitions(in_grows PRIVATE    -D GROWS -D ${GROWS_HASHFCT})

add_executable(step_grows source/step_test.cpp)
target_compile_definitions(step_grows PRIVATE  -D GROWS -D ${GROWS_HASHFCT})

add_executable(in_sim     source/in_test.cpp)
target_compile_definitions(in_sim PRIVATE   -D CUCKOO -D ${GROWS_HASHFCT})

add_executable(step_sim   source/step_test.cpp)
target_compile_definitions(step_sim PRIVATE -D CUCKOO -D ${GROWS_HASHFCT})

add_executable(in_hom     source/in_test.cpp)
target_compile_definitions(in_hom PRIVATE   -D HOM2LVL -D ${GROWS_HASHFCT})

add_executable(step_hom   source/step_test.cpp)
target_compile_definitions(step_hom PRIVATE -D HOM2LVL -D ${GROWS_HASHFCT})
