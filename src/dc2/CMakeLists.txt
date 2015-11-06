add_library(dc2 STATIC dc2.cpp)

target_include_directories(dc2 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(dc2
    pinv)

add_subdirectory(demos)
add_subdirectory(tests)
