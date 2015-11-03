add_executable(dc2-tests
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/gtest/gtest_main.cc
    simple.cpp)

if(UNIX)
    set(PTHREAD "pthread")
else()
    set(PTHREAD "")
endif()

target_link_libraries(dc2-tests
    dc2
    gtest
    ${PTHREAD})
