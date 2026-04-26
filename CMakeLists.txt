cmake_minimum_required(VERSION 3.10)
project(libctkem VERSION 2.0.0 LANGUAGES C ASM)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(include src)

if(NOT DEFINED ARCH)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(CTKEM_SOURCES 
    src/rlwe_kem.c 
    src/randombytes.c
)

if(ARCH MATCHES "aarch64")
    list(APPEND CTKEM_SOURCES asm/arm64/ntt_arm64.S)
else()
    list(APPEND CTKEM_SOURCES ref/ntt_ref.c)
    add_compile_definitions(RLWE_NO_ASM)
endif()

if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -O3 -fstack-protector-strong)
    if(ARCH MATCHES "aarch64")
        add_compile_options(-march=armv8-a+simd)
    endif()
endif()

add_library(ctkem STATIC ${CTKEM_SOURCES})

install(TARGETS ctkem
    EXPORT ctkemTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include)

install(EXPORT ctkemTargets
    FILE ctkemTargets.cmake
    NAMESPACE ctkem::
    DESTINATION lib/cmake/ctkem
)
