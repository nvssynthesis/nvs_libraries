cmake_minimum_required (VERSION 3.15)

project(nvs_libraries LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_definitions(
		_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
)

include_directories(
		include
		external/sprout
)

add_library(nvs_libraries INTERFACE)
target_include_directories(nvs_libraries INTERFACE
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${CMAKE_CURRENT_SOURCE_DIR}/external/sprout
)

enable_testing()
add_executable(nvs_gen_test tests/nvs_gen_test.cpp)
target_link_libraries(nvs_gen_test PRIVATE nvs_libraries)
