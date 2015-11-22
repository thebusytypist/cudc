find_package(CUDA QUIET REQUIRED)

if(WIN32)
    set(DW "-Xcompiler /wd4819")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}${DW})
endif()

cuda_add_library(cudc2 STATIC cudc2.cpp kcudc2.cu)

target_include_directories(cudc2 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(tests)
