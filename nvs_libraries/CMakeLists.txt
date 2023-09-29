cmake_minimum_required (VERSION 3.22)

add_library(nvs_filters INTERFACE)

target_sources(nvs_filters INTERFACE include/nvs_filters.h)	# bad: using include here

target_include_directories(nvs_filters 
						INTERFACE 
								include
)