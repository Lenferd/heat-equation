# ЯВНЫЙ МЕТОД
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})

# Рунхе Хутта
set(target "${APP}_Runge")

file(GLOB src_runge Runge.cpp)

add_executable(${target} ${src_runge})
target_link_libraries(${target} ${LIBRARY})


# Неяvная сhема
set(target "${APP}_Implicit")

file(GLOB src_implicit Implicit.cpp)

add_executable(${target} ${src_implicit})
target_link_libraries(${target} ${LIBRARY})

set_target_properties(
    "${APP}_Euler"
    "${APP}_Runge"
    "${APP}_Implicit"
        PROPERTIES
     RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Sergey_Debug
     RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Sergey_Release
)
