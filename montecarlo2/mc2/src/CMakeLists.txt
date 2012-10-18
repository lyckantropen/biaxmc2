include_directories(../include)
include_directories(../include/3x3-C)
include_directories(${CMAKE_SOURCE_DIR}/boostbase/include)

include_directories(${CMAKE_SOURCE_DIR}/boostbase/include/hashlib)

link_libraries(boostbase)
link_libraries(sqlite3)
add_library(mc2 STATIC mother.cpp MarsagliaRNG.cpp RuntimePropertiesServer.cpp FIFOInterface.cpp 4DSphereRW.cpp  Contractions.cpp  Lattice.cpp  Particle.cpp  Random.cpp  Statistical.cpp  valarray_external.cpp 3x3-C/dsyev2.cpp   3x3-C/dsyevd3.cpp  3x3-C/dsyevj3.cpp  3x3-C/dsyevv3.cpp  3x3-C/slvsec3.cpp 
	3x3-C/dsyevc3.cpp  3x3-C/dsyevh3.cpp  3x3-C/dsyevq3.cpp  3x3-C/dsytrd3.cpp  
	)
