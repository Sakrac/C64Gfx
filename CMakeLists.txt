cmake_minimum_required(VERSION 3.16)
project(C64Gfx LANGUAGES C)

add_executable(C64Gfx src/C64Gfx.c src/gpx_load.c src/gpx_load.h)

target_compile_features(C64Gfx PRIVATE c_std_11)

if(MSVC)
	target_compile_options(C64Gfx PRIVATE /W4)
else()
	target_compile_options(C64Gfx PRIVATE -Wall -Wextra)
endif()

if(UNIX AND NOT APPLE)
	target_link_libraries(C64Gfx PRIVATE m)
endif()