# ЯВНЫЙ МЕТОД
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})

set_target_properties(${target} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Sergey_Debug
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Sergey_Release
)

# НЕЯВНЫЙ МЕТОД
set(target "${APP}_Implicit")

file(GLOB src_implicit Implicit.c)

add_executable(${target} ${src_implicit})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})

set_target_properties(${target} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Sergey_Debug
                      RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Sergey_Release
)
