# ЯВНЫЙ МЕТОД
set(target "${APP}_Euler")

file(GLOB src_euler Euler.c)

add_executable(${target} ${src_euler})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})


# НЕЯВНЫЙ МЕТОД
set(target "${APP}_Implicit")

file(GLOB src_implicit Implicit.c)

add_executable(${target} ${src_implicit})
target_link_libraries(${target} ${LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES})
