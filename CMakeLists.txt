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


foreach(t in step time)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h grows cuckoo hom2lvl triv_grow linprob)
    string(TOUPPER ${h} h_uc)
    add_executable(${t}_${h} source/${t}_test.cpp)
    target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${GROWS_HASHFCT})
  endforeach()
endforeach()
