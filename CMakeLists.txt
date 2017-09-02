cmake_minimum_required (VERSION 2.8)

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

project (space-grow-table)

#### USER DEFINED ##############################################################

set(DYSECT_HASHFCT XXHASH CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE DYSECT_HASHFCT PROPERTY STRINGS XXHASH MURMUR2 MURMUR3 CRC)

set(DYSECT_MALLOC_COUNT OFF CACHE BOOL
  "Display the amount of allocated memory! Needs the malloc_count submodule.")

#### BASIC SETTINGS ############################################################

include_directories(.)

set (CMAKE_CXX_FLAGS "-std=c++14 -msse4.2 -Wall -Wextra -O3 -g")

#### HASH FUNCTIONS ############################################################

if (DYSECT_HASHFCT STREQUAL XXHASH)
  message(STATUS "Looking for xxhash")
  find_path(XXHASH_DIR xxhash/xxhash.h
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/module)
  if (XXHASH_DIR)
    include_directories(${XXHASH_DIR}/xxhash)
    message(STATUS "Looking for xxhash -- found")
  else()
    message(STATUS "Looking for xxhash -- not found")
    message("Cannot find xxHash. Use MURMUR2 instead.")
    set(GROWT_HASHFCT MURMUR2)
  endif()
endif()

if (DYSECT_HASHFCT STREQUAL MURMUR3)
  find_package(smhasher)
  if (SMHASHER_FOUND)
    include_directories(${SMHASHER_INCLUDE_DIRS})
  else()
    message("Cannot find SMHasher which implements MURMUR3. Use MURMUR2 instead.")
    set(DYSECT_HASHFCT MURMUR2)
  endif()
endif()

#### HASH TABLES ###############################################################

set(HASH_TABLES_LIST "cuckoo_standard;cuckoo_standard_inplace;cuckoo_deamortized;dysect;dysect_inplace;cuckoo_independent_2lvl;cuckoo_overlap;cuckoo_overlap_inplace;hopscotch;hopscotch_inplace;robin;robin_inplace;linear_doubling;linear;linear_inplace")

#### THE ONE STEP GROW EXECUTABLE ##############################################

if (DYSECT_MALLOC_COUNT)
  message(STATUS "Looking for malloc_count")
  find_path(MALLOC_COUNT_DIR malloc_count/malloc_count.c
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/module)
  if (MALLOC_COUNT_DIR)
    include_directories(${MALLOC_COUNT_DIR}/malloc_count/)
    message(STATUS "Looking for malloc_count -- found")
  else()
    message(STATUS "Looking for malloc_count -- not found")
    message("Cannot find malloc_count. Turned it off instead.")
    set(DYSECT_MALLOC_COUNT FALSE)
  endif()
endif()

#### CONSTRUCT EXECUTABLE ######################################################
#add_library(mallocc /home/maier/RANDOM/malloc_count/malloc_count.c)

foreach(t time del eps mix crawl mixd)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h ${HASH_TABLES_LIST})
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()

foreach(t sing)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h dysect dysect_inplace)
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()

foreach(t mxls)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h cuckoo_standard cuckoo_standard_inplace)
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()
