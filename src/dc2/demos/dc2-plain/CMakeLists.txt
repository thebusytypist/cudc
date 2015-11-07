add_executable(dc2-plain main.cpp)

if(WIN32)
    target_compile_definitions(dc2-plain
        PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

target_link_libraries(dc2-plain
    dc2)
