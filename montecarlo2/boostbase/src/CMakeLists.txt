include_directories(../include)
include_directories(../include/hashlib)

link_libraries(${SQLITE_LIBRARIES})
add_library(boostbase STATIC serialhash.cpp hl_md5.cpp hl_md5wrapper.cpp)
