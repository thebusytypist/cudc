add_executable(pinv-tests
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/gtest/gtest_main.cc
    simple.cpp)

if(UNIX)
    set(PTHREAD "pthread")
else()
    set(PTHREAD "")
endif()

target_link_libraries(pinv-tests
    pinv
    gtest
    ${PTHREAD})
