cmake_minimum_required (VERSION 3.15)

project(nvs_libraries LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_definitions(
		_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
)

include(FetchContent)

FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG        v3.8.1
)

FetchContent_MakeAvailable(Catch2)

include_directories(
		include
		external/sprout
)

add_library(nvs_libraries INTERFACE)
target_include_directories(nvs_libraries INTERFACE
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${CMAKE_CURRENT_SOURCE_DIR}/external/sprout
)

# Check compiler and add constexpr limits
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	target_compile_options(nvs_libraries INTERFACE -fconstexpr-ops-limit=10000000)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	target_compile_options(nvs_libraries INTERFACE -fconstexpr-steps=10000000)
endif()

enable_testing()

add_executable(nvs_gen_test tests/nvs_gen_test.cpp)
add_executable(nvs_LUT_test tests/nvs_LUT_test.cpp)
target_link_libraries(nvs_gen_test PRIVATE nvs_libraries Catch2::Catch2WithMain)
target_link_libraries(nvs_LUT_test PRIVATE nvs_libraries Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(nvs_gen_test)