include_directories(../include)
include_directories(${CMAKE_SOURCE_DIR}/boostbase/include)
include_directories(${CMAKE_SOURCE_DIR}/boostbase/include/hashlib)

link_libraries(boostbase)
add_library(mc2 SHARED 4DSphereRW.cpp  Contractions.cpp  Lattice.cpp  Particle.cpp  Random.cpp  Statistical.cpp  valarray_external.cpp)
