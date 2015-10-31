add_library(pinv STATIC pinv.cpp)

target_include_directories(pinv PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(tests)
