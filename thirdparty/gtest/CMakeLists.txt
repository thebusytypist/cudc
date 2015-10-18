add_library(gtest STATIC gtest-all.cc)

target_include_directories(gtest PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../)
