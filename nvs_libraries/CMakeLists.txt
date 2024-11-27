cmake_minimum_required (VERSION 3.15)

project(nvs_libraries)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_definitions(
		_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
		#H_ONLY
)

include_directories(
		include
		external/sprout
)

set(SOURCES
		src/nvs_filters.cpp
		src/nvs_LFO.cpp
)

set(HEADERS
		include/nvs_delayFilters.h
		include/nvs_denormal.h
		include/nvs_filters.h
		include/nvs_gen.h
		include/nvs_LFO.h
		include/nvs_logpow.h
		include/nvs_matrix2x2.h
		include/nvs_memoryless.h
)

add_library(nvs_libraries STATIC ${SOURCES} ${HEADERS})

target_include_directories(nvs_libraries PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${CMAKE_CURRENT_SOURCE_DIR}/external/sprout
)

enable_testing()

add_executable(nvs_gen_test tests/nvs_gen_test.cpp)

target_link_libraries(nvs_gen_test PRIVATE nvs_libraries)