add_executable(cudc2-plain main.cpp)

if(WIN32)
    target_compile_definitions(cudc2-plain
        PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

target_link_libraries(cudc2-plain
    cudc2)
