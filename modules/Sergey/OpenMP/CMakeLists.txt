include_directories(include)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)

find_package(AVX REQUIRED)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native ")

set(DIR_NAME ${DIR_AUTHOR}_${DIR_NAME})
set(MODULE      "${DIR_NAME}")
set(APP         "app_${MODULE}")
set(LIBRARY     "lib_${MODULE}")

add_subdirectory(src)
add_subdirectory(app)
