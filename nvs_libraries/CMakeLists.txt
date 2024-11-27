cmake_minimum_required (VERSION 3.15)

project(nvs_libraries)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(nvs_filters
	src/nvs_filters.cpp
)

add_library(nvs_LFO
	src/nvs_LFO.cpp
)

target_include_directories(nvs_filters PUBLIC include)
target_include_directories(nvs_LFO PUBLIC include)

add_library(nvs_gen INTERFACE)
add_library(nvs_logpow INTERFACE)
add_library(nvs_matrix2x2 INTERFACE)
add_library(nvs_memoryless INTERFACE)

target_include_directories(nvs_gen INTERFACE include)
target_include_directories(nvs_logpow INTERFACE include)
target_include_directories(nvs_matrix2x2 INTERFACE include)
target_include_directories(nvs_memoryless INTERFACE include)

enable_testing()

add_executable(nvs_gen_test
		tests/nvs_gen_test.cpp
)

target_link_libraries(nvs_gen_test PRIVATE nvs_gen)

add_test(NAME nvs_gen_test COMMAND nvs_gen_test)