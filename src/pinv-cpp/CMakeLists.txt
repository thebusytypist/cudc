add_library(pinv-cpp STATIC pinv.cpp)

target_include_directories(pinv-cpp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(tests)
