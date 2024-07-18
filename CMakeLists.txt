cmake_minimum_required(VERSION 3.10)

# Project name
project(RepairTools)

# Set the C standard and compilation flags
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -g -m64 -O2")

# Optionally add malloc_count support
option(USE_MALLOC_COUNT "Use malloc_count for memory tracking" OFF)
if(USE_MALLOC_COUNT)
    add_definitions(-DMALLOC_COUNT)
    link_libraries(dl)
    include_directories(tools)
    set(MALLOC_FILES tools/malloc_count.c tools/malloc_count.h)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR})

# Source files
set(SOURCES
    array.c
    hash.c
    heap.c
    records.c
    basics.c
)

# Header files for dependency tracking
set(HEADERS
    hash.h
    heap.h
    records.h
    basics.h
)

# Executable xrepair.x
add_executable(xrepair.x xrepair.c ${SOURCES} ${MALLOC_FILES})

# Executable despair.x
add_executable(despair.x despair.c basics.c ${MALLOC_FILES})

# Test target
add_custom_target(test
    COMMAND ./xrepair.x README
    COMMAND ./despair.x README
    COMMAND cmp README README.out
    DEPENDS xrepair.x despair.x
)

# Clean target not needed, "make clean" provided by CMake