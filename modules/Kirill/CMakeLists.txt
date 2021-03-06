include_directories(include)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_LIST_DIR} NAME)
set(MODULE "${DIR_NAME}")

set(MODULE      "${DIR_NAME}")
set(APP         "app_${MODULE}")
set(LIBRARY     "lib_${MODULE}")

set(THREAD_RUN  0)
set(MKL_RUN     0)
set(AVX2_RUN    0)
set(FPGA_RUN 	0)
set(CPU_CL_RUN  0)
set(GPU_CL_RUN  0)

set(TYPE        DOUBLE)

message(STATUS "------LIB------")
message(STATUS "${LIB}")
message(STATUS "------ARCH-----")
message(STATUS "${ARCH}")

if (LIB)
	if (${LIB} MATCHES "mkl")
		set(MKL_RUN      1)
	elseif(${LIB} MATCHES "avx")
		set(AVX2_RUN     1)
	elseif(${LIB} MATCHES "omp")
		set(THREAD_RUN   1)
	elseif(ARCH)
		if(${LIB} MATCHES "opencl" AND ${ARCH} MATCHES "gpu")
			set(GPU_CL_RUN   1)
		elseif(${LIB} MATCHES "opencl" AND ${ARCH} MATCHES "fpga")
			set(FPGA_RUN     1)
		elseif(${LIB} MATCHES "opencl" AND ${ARCH} MATCHES "cpu")
			set(CPU_CL_RUN   1)
		endif()
	endif()
endif()

	#   Setting paths
set(INPUT_EULER_SETTING_PATH    "${CMAKE_SOURCE_DIR}/initial/setting3.ini")
set(INPUT_EULER_FUNCTION_PATH   "${CMAKE_SOURCE_DIR}/initial/function3.txt")
set(RESULT_EULER_PATH           "${CMAKE_SOURCE_DIR}/result/Kirill/${OUT_FILE}")

if(${FPGA_RUN})
	set(KERNEL_CL_PATH    		"${CMAKE_BINARY_DIR}/modules/Kirill/src/kernel.aocx")
elseif(${CPU_CL_RUN} OR ${GPU_CL_RUN})
	set(KERNEL_CL_PATH    		"${CMAKE_BINARY_DIR}/modules/Kirill/src/kernel.cl")
endif()
#

find_package(AVX)
find_package(MKL)
if(${GPU_CL_RUN} OR ${CPU_CL_RUN})
	find_package(OpenCL)
elseif(${FPGA_RUN})
	find_package(AlteraCL)
endif()

if(${AVX2_RUN} AND AVX_FOUND)
    set(AVX2_RUN, 1)
    message(STATUS "AVX2 RUN")
elseif(${MKL_RUN} AND MKL_FOUND)
    set(AVX2_RUN, 0)
    set(MKL_RUN, 1)
    message(STATUS "MKL RUN")
elseif(${GPU_CL_RUN} AND OPENCL_FOUND)
    set(MKL_RUN, 0)
	set(GPU_CL_RUN, 1)
	message(STATUS "GPU_CL_RUN")
elseif(${CPU_CL_RUN} AND OPENCL_FOUND)
	set(MKL_RUN, 0)
	set(CPU_CL_RUN, 1)
	message(STATUS "CPU_CL_RUN")
else()
	set(CPU_CL_RUN, 0)
	set(GPU_CL_RUN, 0)
endif()

configure_file("include/utils/ts.h.in" "${CMAKE_CURRENT_LIST_DIR}/include/utils/ts.h")

add_subdirectory(src)
add_subdirectory(app)
