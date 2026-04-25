cmake_minimum_required(VERSION 3.10)
project(libctkem VERSION 2.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -O3 -fstack-protector-strong)
endif()

add_library(ctkem STATIC 
    src/rlwe_kem.c
)

target_include_directories(ctkem PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

add_executable(ctkem_test src/rlwe_kem.c)
target_compile_definitions(ctkem_test PRIVATE RLWE_TESTING) 

install(TARGETS ctkem
    EXPORT ctkemTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES src/rlwe_kem.h DESTINATION include)

install(EXPORT ctkemTargets
    FILE ctkemTargets.cmake
    NAMESPACE ctkem::
    DESTINATION lib/cmake/ctkem
)
