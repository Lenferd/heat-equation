# Явный метод Эйлера
set(target_euler_name "${APP}_Euler")
set(target "${APP}_Euler")

file(GLOB src_euler Euler.cpp)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})

## Euler AVX_SHUFFLE
set(target "${target_euler_name}_AVX")
add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_compile_definitions(${target} PRIVATE AVX_SHUFFLE)

## Euler NAIVE
set(target "${target_euler_name}_NAIVE")
add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_compile_definitions(${target} PRIVATE NAIVE)

## Euler AVX_NAIVE
set(target "${target_euler_name}_AVX_NAIVE")
add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_compile_definitions(${target} PRIVATE AVX_NAIVE)

# Метод Рунге-Кутты 4 порядка
set(target "${APP}_Runge")

file(GLOB src_runge Runge.cpp)

add_executable(${target} ${src_runge})
target_link_libraries(${target} ${LIBRARY})


# Неявная схема, метод Эйлера
set(target "${APP}_Implicit")

file(GLOB src_implicit Implicit.cpp)

add_executable(${target} ${src_implicit})
target_link_libraries(${target} ${LIBRARY})

set_target_properties(
    "${APP}_Euler"
    "${APP}_Runge"
    "${APP}_Implicit"
    "${target_euler_name}_AVX"
    "${target_euler_name}_NAIVE"
    "${target_euler_name}_AVX_NAIVE"
        PROPERTIES
     RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Sergey_Debug
     RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Sergey_Release
)
