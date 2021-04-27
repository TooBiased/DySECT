cmake_minimum_required (VERSION 2.8)

# custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

project (space-grow-table)

#### USER DEFINED ##############################################################

set(DYSECT_BUILD_MODE RELEASE CACHE STRING
  "the build mode of the current tests (controls optimizations and symbols)!")
set_property(CACHE DYSECT_BUILD_MODE PROPERTY STRINGS RELEASE DEBUG REL_W_SYMBOLS )

set(DYSECT_HASHFCT XXH3 CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE DYSECT_HASHFCT PROPERTY STRINGS XXH3 XXHASH MURMUR2 MURMUR3 CRC)

set(DYSECT_MALLOC_COUNT OFF CACHE BOOL
  "Display the amount of allocated memory! Needs the malloc_count submodule.")

set(DYSECT_RSS_COUNT OFF CACHE BOOL
  "Display the amount of physically used memory pages (resident set size must be linux + no malloc_count)!")

set(DYSECT_CUCKOO_PREFETCH ON CACHE BOOL
  "Use prefetching in conjunction with accessing cuckoo buckets")

#### BASIC SETTINGS ############################################################

include_directories(.)

#set (CMAKE_CXX_FLAGS "-std=c++17 -msse4.2 -Wall -Wextra -O3 -g -march=native")
set (FLAGS "-std=c++17 -msse4.2 -Wall -Wextra")

if (DYSECT_BUILD_MODE STREQUAL DEBUG)
set (FLAGS "${FLAGS} -g3 -ggdb -O0")
endif()
if (DYSECT_BUILD_MODE STREQUAL RELEASE)
set (FLAGS "${FLAGS} -march=native -O3")
endif()
if (DYSECT_BUILD_MODE STREQUAL REL_W_SYMBOLS)
set (FLAGS "${FLAGS} -g3 -ggdb -march=native -O3")
endif()


#### FIND PTHREAD (less variance after pin to core) ############################

find_package(Threads REQUIRED)
set(TEST_DEP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${TEST_DEP_LIBRARIES})

#### HASH FUNCTIONS ############################################################

if ((DYSECT_HASHFCT STREQUAL XXHASH) OR
    (DYSECT_HASHFCT STREQUAL XXH3))
  message(STATUS "Looking for xxhash")
  find_path(XXHASH_DIR xxhash/xxhash.h
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/module)
  if (XXHASH_DIR)
    include_directories(${XXHASH_DIR}/xxhash)
    message(STATUS "Looking for xxhash -- found")
  else()
    message(STATUS "Looking for xxhash -- not found")
    message("Cannot find xxHash. Use MURMUR2 instead.")
    set(DYSECT_HASHFCT MURMUR2)
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

set(HASH_TABLES_LIST "multi_dysect;multi_dysect_inplace;multi_cuckoo_standard;multi_cuckoo_standard_inplace;multi_cuckoo_deamortized;multi_cuckoo_independent_2lvl;multi_cuckoo_overlap;multi_cuckoo_overlap_inplace;hop_hopscotch;hop_hopscotch_inplace;triv_robin;triv_robin_inplace;triv_multitable_robin;triv_linear;triv_linear_inplace;triv_multitable_linear;triv_quadratic;triv_multitable_quadratic;triv_quadratic_inplace;triv_chaining")

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
set_target_properties(example PROPERTIES COMPILE_FLAGS "${FLAGS}")



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
    set_target_properties(${t}_${h} PROPERTIES COMPILE_FLAGS "${FLAGS}")
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
    if (DYSECT_CUCKOO_PREFETCH)
      target_compile_definitions(${t}_${h} PRIVATE -D PREFETCH)
    endif()
  endforeach()
endforeach()

foreach(t sing)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h multi_dysect multi_dysect_inplace)
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    set_target_properties(${t}_${h} PROPERTIES COMPILE_FLAGS "${FLAGS}")
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
  endforeach()
endforeach()

foreach(t mxls)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${t})
  foreach(h multi_cuckoo_standard multi_cuckoo_standard_inplace)
    string(TOUPPER ${h} h_uc)
    if (DYSECT_MALLOC_COUNT)
      add_executable(${t}_${h} source/${t}_test.cpp ${MALLOC_COUNT_DIR}/malloc_count/malloc_count.c)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT} -D MALLOC_COUNT)
    else()
      add_executable(${t}_${h} source/${t}_test.cpp)
      target_compile_definitions(${t}_${h} PRIVATE -D ${h_uc} -D ${DYSECT_HASHFCT})
    endif()
    set_target_properties(${t}_${h} PROPERTIES COMPILE_FLAGS "${FLAGS}")
    target_link_libraries(${t}_${h} ${TEST_DEP_LIBRARIES} dl)
  endforeach()
endforeach()

add_custom_target(all_time)
add_dependencies(all_time
  time_multi_dysect
  time_multi_dysect_inplace
  time_multi_cuckoo_standard
  time_multi_cuckoo_standard_inplace
  time_multi_cuckoo_deamortized
  time_multi_cuckoo_independent_2lvl
  time_multi_cuckoo_overlap
  time_multi_cuckoo_overlap_inplace
  time_hop_hopscotch
  time_hop_hopscotch_inplace
  time_triv_robin
  time_triv_robin_inplace
  time_triv_multitable_robin
  time_triv_linear
  time_triv_linear_inplace
  time_triv_multitable_linear
  time_triv_quadratic
  time_triv_chaining)

add_custom_target(all_displ)
add_dependencies(all_displ
  displ_multi_dysect
  displ_multi_dysect_inplace
  displ_multi_cuckoo_standard
  displ_multi_cuckoo_standard_inplace
  displ_multi_cuckoo_deamortized
  displ_multi_cuckoo_independent_2lvl
  displ_multi_cuckoo_overlap
  displ_multi_cuckoo_overlap_inplace
  displ_hop_hopscotch
  displ_hop_hopscotch_inplace
  displ_triv_robin
  displ_triv_robin_inplace
  displ_triv_multitable_robin
  displ_triv_linear
  displ_triv_linear_inplace
  displ_triv_multitable_linear
  displ_triv_quadratic
  displ_triv_chaining)

add_custom_target(all_del)
add_dependencies(all_del
  del_multi_dysect
  del_multi_dysect_inplace
  del_multi_cuckoo_standard
  del_multi_cuckoo_standard_inplace
  del_multi_cuckoo_deamortized
  del_multi_cuckoo_independent_2lvl
  del_multi_cuckoo_overlap
  del_multi_cuckoo_overlap_inplace
  del_hop_hopscotch
  del_hop_hopscotch_inplace
  del_triv_robin
  del_triv_robin_inplace
  del_triv_multitable_robin
  del_triv_linear
  del_triv_linear_inplace
  del_triv_multitable_linear
  del_triv_quadratic
  del_triv_chaining)

add_custom_target(all_eps)
add_dependencies(all_eps
  eps_multi_dysect
  eps_multi_dysect_inplace
  eps_multi_cuckoo_standard
  eps_multi_cuckoo_standard_inplace
  eps_multi_cuckoo_deamortized
  eps_multi_cuckoo_independent_2lvl
  eps_multi_cuckoo_overlap
  eps_multi_cuckoo_overlap_inplace
  eps_hop_hopscotch
  eps_hop_hopscotch_inplace
  eps_triv_robin
  eps_triv_robin_inplace
  eps_triv_multitable_robin
  eps_triv_linear
  eps_triv_linear_inplace
  eps_triv_multitable_linear
  eps_triv_quadratic
  eps_triv_chaining)

add_custom_target(all_mix)
add_dependencies(all_mix
  mix_multi_dysect
  mix_multi_dysect_inplace
  mix_multi_cuckoo_standard
  mix_multi_cuckoo_standard_inplace
  mix_multi_cuckoo_deamortized
  mix_multi_cuckoo_independent_2lvl
  mix_multi_cuckoo_overlap
  mix_multi_cuckoo_overlap_inplace
  mix_hop_hopscotch
  mix_hop_hopscotch_inplace
  mix_triv_robin
  mix_triv_robin_inplace
  mix_triv_multitable_robin
  mix_triv_linear
  mix_triv_linear_inplace
  mix_triv_multitable_linear
  mix_triv_quadratic
  mix_triv_chaining)

add_custom_target(all_crawl)
add_dependencies(all_crawl
  crawl_multi_dysect
  crawl_multi_dysect_inplace
  crawl_multi_cuckoo_standard
  crawl_multi_cuckoo_standard_inplace
  crawl_multi_cuckoo_deamortized
  crawl_multi_cuckoo_independent_2lvl
  crawl_multi_cuckoo_overlap
  crawl_multi_cuckoo_overlap_inplace
  crawl_hop_hopscotch
  crawl_hop_hopscotch_inplace
  crawl_triv_robin
  crawl_triv_robin_inplace
  crawl_triv_multitable_robin
  crawl_triv_linear
  crawl_triv_linear_inplace
  crawl_triv_multitable_linear
  crawl_triv_quadratic
  crawl_triv_chaining)

add_custom_target(all_mixd)
add_dependencies(all_mixd
  mixd_multi_dysect
  mixd_multi_dysect_inplace
  mixd_multi_cuckoo_standard
  mixd_multi_cuckoo_standard_inplace
  mixd_multi_cuckoo_deamortized
  mixd_multi_cuckoo_independent_2lvl
  mixd_multi_cuckoo_overlap
  mixd_multi_cuckoo_overlap_inplace
  mixd_hop_hopscotch
  mixd_hop_hopscotch_inplace
  mixd_triv_robin
  mixd_triv_robin_inplace
  mixd_triv_multitable_robin
  mixd_triv_linear
  mixd_triv_linear_inplace
  mixd_triv_multitable_linear
  mixd_triv_quadratic
  mixd_triv_chaining)

add_custom_target(all_linear)
add_dependencies(all_linear
  time_triv_linear
  del_triv_linear
  eps_triv_linear
  mix_triv_linear
  mixd_triv_linear
  crawl_triv_linear
  displ_triv_linear)

add_custom_target(all_linear_inplace)
add_dependencies(all_linear_inplace
  time_triv_linear_inplace
  del_triv_linear_inplace
  eps_triv_linear_inplace
  mix_triv_linear_inplace
  mixd_triv_linear_inplace
  crawl_triv_linear_inplace
  displ_triv_linear_inplace)

add_custom_target(all_cuckoo_standard)
add_dependencies(all_cuckoo_standard
  time_multi_cuckoo_standard
  del_multi_cuckoo_standard
  eps_multi_cuckoo_standard
  mix_multi_cuckoo_standard
  mixd_multi_cuckoo_standard
  crawl_multi_cuckoo_standard
  displ_multi_cuckoo_standard)

add_custom_target(all_cuckoo_standard_inplace)
add_dependencies(all_cuckoo_standard_inplace
  time_multi_cuckoo_standard_inplace
  del_multi_cuckoo_standard_inplace
  eps_multi_cuckoo_standard_inplace
  mix_multi_cuckoo_standard_inplace
  mixd_multi_cuckoo_standard_inplace
  crawl_multi_cuckoo_standard_inplace
  displ_multi_cuckoo_standard_inplace)

add_custom_target(all_robin)
add_dependencies(all_robin
  time_triv_robin
  del_triv_robin
  eps_triv_robin
  mix_triv_robin
  mixd_triv_robin
  crawl_triv_robin
  displ_triv_robin)

add_custom_target(all_robin_inplace)
add_dependencies(all_robin_inplace
  time_triv_robin_inplace
  del_triv_robin_inplace
  eps_triv_robin_inplace
  mix_triv_robin_inplace
  mixd_triv_robin_inplace
  crawl_triv_robin_inplace
  displ_triv_robin_inplace)

add_custom_target(all_dysect)
add_dependencies(all_dysect
  time_multi_dysect
  del_multi_dysect
  eps_multi_dysect
  mix_multi_dysect
  mixd_multi_dysect
  crawl_multi_dysect
  displ_multi_dysect)

add_custom_target(all_dysect_inplace)
add_dependencies(all_dysect_inplace
  time_multi_dysect_inplace
  del_multi_dysect_inplace
  eps_multi_dysect_inplace
  mix_multi_dysect_inplace
  mixd_multi_dysect_inplace
  crawl_multi_dysect_inplace
  displ_multi_dysect_inplace)

add_custom_target(all_quadratic)
add_dependencies(all_quadratic
  time_triv_quadratic
  del_triv_quadratic
  eps_triv_quadratic
  mix_triv_quadratic
  mixd_triv_quadratic
  crawl_triv_quadratic
  displ_triv_quadratic)

add_custom_target(diss_targets)
add_dependencies(diss_targets
  time_multi_dysect
  time_multi_dysect_inplace
  time_multi_cuckoo_standard_inplace
  time_multi_cuckoo_independent_2lvl
  time_triv_robin_inplace
  time_triv_multitable_robin
  time_triv_linear_inplace
  time_triv_multitable_linear
  time_triv_quadratic_inplace
  time_triv_multitable_quadratic

  displ_multi_cuckoo_standard
  displ_hop_hopscotch
  displ_triv_robin
  displ_triv_linear
  displ_triv_quadratic
  displ_triv_chaining

  del_multi_dysect
  del_multi_dysect_inplace
  del_multi_cuckoo_standard_inplace
  del_multi_cuckoo_independent_2lvl
  del_triv_robin_inplace
  del_triv_multitable_robin
  del_triv_linear_inplace
  del_triv_multitable_linear

  eps_multi_dysect
  eps_multi_dysect_inplace
  eps_multi_cuckoo_standard_inplace
  eps_triv_robin_inplace
  eps_triv_linear_inplace
  eps_triv_quadratic_inplace

  crawl_multi_dysect
  crawl_multi_dysect_inplace
  crawl_multi_cuckoo_standard_inplace
  crawl_multi_cuckoo_independent_2lvl
  crawl_triv_robin_inplace
  crawl_triv_multitable_robin
  crawl_triv_linear_inplace
  crawl_triv_multitable_linear
  crawl_triv_quadratic_inplace
  crawl_triv_multitable_quadratic

  sing_multi_dysect)
