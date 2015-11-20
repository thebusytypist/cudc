add_executable(cudc2-tests
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/gtest/gtest_main.cc
    simple.cpp)

if(UNIX)
    set(PTHREAD "pthread")
else()
    set(PTHREAD "")
endif()

target_link_libraries(cudc2-tests
    cudc2
    gtest
    ${PTHREAD})
