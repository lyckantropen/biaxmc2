include_directories(../include)
include_directories(../include/3x3-C)
include_directories(${CMAKE_SOURCE_DIR}/boostbase/include)

include_directories(${CMAKE_SOURCE_DIR}/boostbase/include/hashlib)

if(CMAKE_USE_PTHREADS_INIT)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
endif(CMAKE_USE_PTHREADS_INIT)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")

add_library(mc2 SHARED omp_debug.cpp serialhash.cpp hl_md5.cpp hl_md5wrapper.cpp SimulationDBFind.cpp Metropolis.cpp AutoCorrelationTimeCalculator.cpp ILoggable.cpp mother.cpp MarsagliaRNG.cpp RuntimePropertiesServer.cpp FIFOInterface.cpp 4DSphereRW.cpp  Contractions.cpp  Lattice.cpp  Particle.cpp  Random.cpp  Statistical.cpp  Maths.cpp 3x3-C/dsyev2.cpp   3x3-C/dsyevd3.cpp  3x3-C/dsyevj3.cpp  3x3-C/dsyevv3.cpp  3x3-C/slvsec3.cpp 
	3x3-C/dsyevc3.cpp  3x3-C/dsyevh3.cpp  3x3-C/dsyevq3.cpp  3x3-C/dsytrd3.cpp  
	)
