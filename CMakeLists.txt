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

set(DYSECT_RSS_COUNT OFF CACHE BOOL
  "Display the amount of physically used memory pages (resident set size must be linux + no malloc_count)!")

#### BASIC SETTINGS ############################################################

include_directories(.)

set (CMAKE_CXX_FLAGS "-std=c++17 -msse4.2 -Wall -Wextra -O3 -g")

#### FIND PTHREAD (less variance after pin to core) ############################

find_package(Threads REQUIRED)
set(TEST_DEP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${TEST_DEP_LIBRARIES})

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

set(HASH_TABLES_LIST "dysect;dysect_inplace;cuckoo_standard;cuckoo_standard_inplace;cuckoo_deamortized;cuckoo_independent_2lvl;cuckoo_overlap;cuckoo_overlap_inplace;hopscotch;hopscotch_inplace;robin;robin_inplace;multitable_robin;linear_doubling;linear;linear_inplace;multitable_linear;quadratic")

#### LOOKS FOR THE MALLOC COUNTING LIB #########################################

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

#### BUILD THE EXAMPLE #########################################################
add_executable(example example/example.cpp)


#### CONSTRUCT EXECUTABLE ######################################################
#add_library(mallocc /home/maier/RANDOM/malloc_count/malloc_count.c)

foreach(t time del eps mix crawl mixd displ)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h ${HASH_TABLES_LIST})
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    elseif (DYSECT_RSS_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D RSS_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
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
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
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
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
  endforeach()
endforeach()

add_custom_target(all_time)
add_dependencies(all_time
  time_dysect
  time_dysect_inplace
  time_cuckoo_standard
  time_cuckoo_standard_inplace
  time_cuckoo_deamortized
  time_cuckoo_independent_2lvl
  time_cuckoo_overlap
  time_cuckoo_overlap_inplace
  time_hopscotch
  time_hopscotch_inplace
  time_robin
  time_robin_inplace
  time_multitable_robin
  time_linear_doubling
  time_linear
  time_linear_inplace
  time_multitable_linear
  time_quadratic)

add_custom_target(all_displ)
add_dependencies(all_displ
  displ_dysect
  displ_dysect_inplace
  displ_cuckoo_standard
  displ_cuckoo_standard_inplace
  displ_cuckoo_deamortized
  displ_cuckoo_independent_2lvl
  displ_cuckoo_overlap
  displ_cuckoo_overlap_inplace
  displ_hopscotch
  displ_hopscotch_inplace
  displ_robin
  displ_robin_inplace
  displ_multitable_robin
  displ_linear_doubling
  displ_linear
  displ_linear_inplace
  displ_multitable_linear
  displ_quadratic)


add_custom_target(all_del)
add_dependencies(all_del
  del_dysect
  del_dysect_inplace
  del_cuckoo_standard
  del_cuckoo_standard_inplace
  del_cuckoo_deamortized
  del_cuckoo_independent_2lvl
  del_cuckoo_overlap
  del_cuckoo_overlap_inplace
  del_hopscotch
  del_hopscotch_inplace
  del_robin
  del_robin_inplace
  del_multitable_robin
  del_linear_doubling
  del_linear
  del_linear_inplace
  del_multitable_linear
  del_quadratic)

add_custom_target(all_eps)
add_dependencies(all_eps
  eps_dysect
  eps_dysect_inplace
  eps_cuckoo_standard
  eps_cuckoo_standard_inplace
  eps_cuckoo_deamortized
  eps_cuckoo_independent_2lvl
  eps_cuckoo_overlap
  eps_cuckoo_overlap_inplace
  eps_hopscotch
  eps_hopscotch_inplace
  eps_robin
  eps_robin_inplace
  eps_multitable_robin
  eps_linear_doubling
  eps_linear
  eps_linear_inplace
  eps_multitable_linear
  eps_quadratic)

add_custom_target(all_mix)
add_dependencies(all_mix
  mix_dysect
  mix_dysect_inplace
  mix_cuckoo_standard
  mix_cuckoo_standard_inplace
  mix_cuckoo_deamortized
  mix_cuckoo_independent_2lvl
  mix_cuckoo_overlap
  mix_cuckoo_overlap_inplace
  mix_hopscotch
  mix_hopscotch_inplace
  mix_robin
  mix_robin_inplace
  mix_multitable_robin
  mix_linear_doubling
  mix_linear
  mix_linear_inplace
  mix_multitable_linear
  mix_quadratic)

add_custom_target(all_crawl)
add_dependencies(all_crawl
  crawl_dysect
  crawl_dysect_inplace
  crawl_cuckoo_standard
  crawl_cuckoo_standard_inplace
  crawl_cuckoo_deamortized
  crawl_cuckoo_independent_2lvl
  crawl_cuckoo_overlap
  crawl_cuckoo_overlap_inplace
  crawl_hopscotch
  crawl_hopscotch_inplace
  crawl_robin
  crawl_robin_inplace
  crawl_multitable_robin
  crawl_linear_doubling
  crawl_linear
  crawl_linear_inplace
  crawl_multitable_linear
  crawl_quadratic)

add_custom_target(all_mixd)
add_dependencies(all_mixd
  mixd_dysect
  mixd_dysect_inplace
  mixd_cuckoo_standard
  mixd_cuckoo_standard_inplace
  mixd_cuckoo_deamortized
  mixd_cuckoo_independent_2lvl
  mixd_cuckoo_overlap
  mixd_cuckoo_overlap_inplace
  mixd_hopscotch
  mixd_hopscotch_inplace
  mixd_robin
  mixd_robin_inplace
  mixd_multitable_robin
  mixd_linear_doubling
  mixd_linear
  mixd_linear_inplace
  mixd_multitable_linear
  mixd_quadratic)

add_custom_target(all_linear)
add_dependencies(all_linear
  time_linear
  del_linear
  eps_linear
  mix_linear
  mixd_linear
  crawl_linear
  displ_linear)

add_custom_target(all_linear_inplace)
add_dependencies(all_linear_inplace
  time_linear_inplace
  del_linear_inplace
  eps_linear_inplace
  mix_linear_inplace
  mixd_linear_inplace
  crawl_linear_inplace
  displ_linear_inplace)

add_custom_target(all_cuckoo_standard)
add_dependencies(all_cuckoo_standard
  time_cuckoo_standard
  del_cuckoo_standard
  eps_cuckoo_standard
  mix_cuckoo_standard
  mixd_cuckoo_standard
  crawl_cuckoo_standard
  displ_cuckoo_standard)

add_custom_target(all_cuckoo_standard_inplace)
add_dependencies(all_cuckoo_standard_inplace
  time_cuckoo_standard_inplace
  del_cuckoo_standard_inplace
  eps_cuckoo_standard_inplace
  mix_cuckoo_standard_inplace
  mixd_cuckoo_standard_inplace
  crawl_cuckoo_standard_inplace
  displ_cuckoo_standard_inplace)

add_custom_target(all_robin)
add_dependencies(all_robin
  time_robin
  del_robin
  eps_robin
  mix_robin
  mixd_robin
  crawl_robin
  displ_robin)

add_custom_target(all_robin_inplace)
add_dependencies(all_robin_inplace
  time_robin_inplace
  del_robin_inplace
  eps_robin_inplace
  mix_robin_inplace
  mixd_robin_inplace
  crawl_robin_inplace
  displ_robin_inplace)

add_custom_target(all_dysect)
add_dependencies(all_dysect
  time_dysect
  del_dysect
  eps_dysect
  mix_dysect
  mixd_dysect
  crawl_dysect
  displ_dysect)

add_custom_target(all_dysect_inplace)
add_dependencies(all_dysect_inplace
  time_dysect_inplace
  del_dysect_inplace
  eps_dysect_inplace
  mix_dysect_inplace
  mixd_dysect_inplace
  crawl_dysect_inplace
  displ_dysect_inplace)

add_custom_target(all_quadratic)
add_dependencies(all_quadratic
  time_quadratic
  del_quadratic
  eps_quadratic
  mix_quadratic
  mixd_quadratic
  crawl_quadratic
  displ_quadratic)
