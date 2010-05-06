include_directories(${CMAKE_SOURCE_DIR}/boostbase/include)
include_directories(${CMAKE_SOURCE_DIR}/boostbase/include/hashlib)
include_directories(${CMAKE_SOURCE_DIR}/mc2/include)
link_directories(${CMAKE_SOURCE_DIR}/lib)

link_libraries(boostbase mc2)
add_executable(simulation simulation.cpp)
add_executable(dbmine mc2mine.cpp)
add_executable(recalculate recalculate.cpp)
