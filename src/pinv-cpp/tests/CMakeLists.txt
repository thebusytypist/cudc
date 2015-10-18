add_executable(pinv-cpp-tests
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../thirdparty/gtest/gtest_main.cc
    simple.cpp)

target_link_libraries(pinv-cpp-tests
    pinv-cpp
    gtest)
