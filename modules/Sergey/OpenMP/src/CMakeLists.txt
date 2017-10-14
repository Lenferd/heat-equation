set(target ${LIBRARY})

file(GLOB hrcs "../include/*.h")
file(GLOB srcs "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
add_definitions(-mavx -o2)
add_library(${target} STATIC ${hrcs} ${srcs})
