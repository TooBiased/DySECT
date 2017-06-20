cmake_minimum_required (VERSION 2.8)

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

project (space-grow-table)

#### USER DEFINED ##############################################################

set(GROWS_HASHFCT XXHASH CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE GROWS_HASHFCT PROPERTY STRINGS XXHASH MURMUR2 MURMUR3 CRC)

#### BASIC SETTINGS ############################################################

include_directories(.)

set (CMAKE_CXX_FLAGS "-std=c++14 -msse4.2 -Wall -Wextra -O0 -g")

#### HASH FUNCTIONS ############################################################

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

#### HASH TABLES ###############################################################

set(DYSECT_LIST "ceg2l;ceipg2l")
set(HASH_TABLES_LIST "ceg2l;csimple;cover;coverip;chomogeneous2l;cindependent2l;cinplace;ceipg2l;strctcuckoo")
set(NON_CUCKOO_LIST  "linprob;spaceprob;robinprob;ripprob;hopprob;hipprob;std_unordered;inplaceprob;testtab;strctprob;strctrobin")

find_package(tslHopscotch)
if (TSLHOPSCOTCH_FOUND)
  include_directories(${TSLHOPSCOTCH_INCLUDE_DIRS})
  set(NON_CUCKOO_LIST "${NON_CUCKOO_LIST};hopscotch;spacehopscotch")
endif()

#### THE ONE STEP GROW EXECUTABLE ##############################################
include_directories("/home/maier/RANDOM/malloc_count/")

#Currently not up to date due to multichange + ... time

#add_executable(step_all source/step_test.cpp)
#target_compile_definitions(step_all PRIVATE -D GROWS -D ${GROWS_HASHFCT})



#### CONSTRUCT EXECUTABLE ######################################################
#add_library(mallocc /home/maier/RANDOM/malloc_count/malloc_count.c)

foreach(t time del mix reg eps crawl mixd)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h ${HASH_TABLES_LIST})  #grows cuckoo hom2lvl triv_grow linprob)
    string(TOUPPER ${h} h_uc)
    add_executable(${t}_${h} source/${t}_test.cpp /home/maier/RANDOM/malloc_count/malloc_count.c)
    #add_dependencies(${t}_${h} mallocc)
    target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${GROWS_HASHFCT} -D MALLOC_COUNT)
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()

foreach(t time del mix reg eps crawl mixd)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h ${NON_CUCKOO_LIST})
    string(TOUPPER ${h} h_uc)
    add_executable(${t}_${h} source/${t}_test.cpp /home/maier/RANDOM/malloc_count/malloc_count.c)
    target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${GROWS_HASHFCT} -D MALLOC_COUNT)
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()

foreach(t sing)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h ${DYSECT_LIST})
    string(TOUPPER ${h} h_uc)
    add_executable(${t}_${h} source/${t}_test.cpp /home/maier/RANDOM/malloc_count/malloc_count.c)
    #add_dependencies(${t}_${h} mallocc)
    target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${GROWS_HASHFCT} -D MALLOC_COUNT)
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()

foreach(t mxls)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h csimple cinplace)
    string(TOUPPER ${h} h_uc)
    add_executable(${t}_${h} source/${t}_test.cpp /home/maier/RANDOM/malloc_count/malloc_count.c)
    #add_dependencies(${t}_${h} mallocc)
    target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${GROWS_HASHFCT} -D MALLOC_COUNT)
    target_link_libraries(${t}_${h} dl)
  endforeach()
endforeach()
