cmake_minimum_required(VERSION 2.6)
include(FindPkgConfig)
project(montecarlo2)

set(Boost_ADDITIONAL_VERSIONS "1.41.0" "1.34.1" "1.35.0" "1.35.1" "1.36.0" "1.38.0" "1.39.0")
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.38.0 COMPONENTS thread filesystem program_options system date_time iostreams REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
link_libraries(${Boost_LIBRARIES})

#pkg_check_modules(SQLITE sqlite3)
find_package(ZLIB REQUIRED)

link_libraries(${ZLIB_LIBRARIES})
link_libraries(dl)
link_libraries(${CMAKE_THREAD_LIBS_INIT})
#link_libraries(${SQLITE_LIBRARIES})
#link_libraries("/home/trojanowski/usr/lib/libsqlite3.a")
#link_libraries("/home/trojanowski/usr/lib/libsqlite3.a")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSQLITE_TEMP_STORE=3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSQLITE_TEMP_STORE=3")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARIY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
link_directories(${PROJECT_SOURCE_DIR}/bin)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -msse4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -msse4.1")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -msse4.1 -std=c99")

add_subdirectory(boostbase)
add_subdirectory(mc2)
add_subdirectory(tools)
