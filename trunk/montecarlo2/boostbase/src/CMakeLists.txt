include_directories(../include)
include_directories(../include/hashlib)

add_library(boostbase SHARED serialhash.cpp hl_md5.cpp hl_md5wrapper.cpp)
