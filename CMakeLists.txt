 #
 # breeze2d - a model to experiment with aerosol emission impact on dynamics.
 #
 # Copyright (C) 2012 Dmitry Mikushin
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #

cmake_minimum_required (VERSION 2.6)
project (breeze2d)

set (breeze2d_VERSION_MAJOR 1)
set (breeze2d_VERSION_MINOR 0)

INCLUDE(CheckIncludeFile)
INCLUDE(CheckLibraryExists)

find_package(LAPACK REQUIRED)
 
option(HAVE_SINGLE "Use single precision floating-point computations" ON)
option(HAVE_DOUBLE "Use double precision floating-point computations" OFF)

if (HAVE_SINGLE AND HAVE_DOUBLE)
	message(FATAL_ERROR "The single or double floating-precision must be set")
endif (HAVE_SINGLE AND HAVE_DOUBLE)

if ((NOT HAVE_SINGLE) AND (NOT HAVE_DOUBLE))
	message(FATAL_ERROR "The single or double floating-precision must be set")
endif ((NOT HAVE_SINGLE) AND (NOT HAVE_DOUBLE))

if (HAVE_SINGLE)
	add_definitions(-DHAVE_SINGLE)
endif (HAVE_SINGLE)

if (HAVE_DOUBLE)
	add_definitions(-DHAVE_DOUBLE)
endif (HAVE_DOUBLE)

option(HAVE_FFTW "Use FFTW Fast Fourier Transform library" OFF)
option(HAVE_FFTW_MKL "Use MKL FFTW-compatible Fast Fourier Transform library" OFF)

if ((NOT HAVE_FFTW) AND (NOT HAVE_FFTW_MKL))
	message(FATAL_ERROR "The FFTW or MKL library must be turned on as FFT backend")
endif ((NOT HAVE_FFTW) AND (NOT HAVE_FFTW_MKL))

if ((HAVE_FFTW AND HAVE_FFTW_MKL))
	message(FATAL_ERROR "The FFTW or MKL library must be turned on as FFT backend")
endif ((HAVE_FFTW AND HAVE_FFTW_MKL))

if (HAVE_FFTW)
	CHECK_INCLUDE_FILE(fftw3.h HAVE_FFTW_HEADER)
	if (NOT HAVE_FFTW_HEADER)
   		message(FATAL_ERROR "Cannot find fftw3.h")
	endif (NOT HAVE_FFTW_HEADER)
	include_directories(HAVE_FFTW_HEADER)
	if (HAVE_SINGLE)
		find_library(HAVE_FFTWF_THREADS_LIBRARY NAMES fftw3f_threads HINTS "${HAVE_FFTWF_THREADS_LIBRARY}")
		find_library(HAVE_FFTWF_OMP_LIBRARY NAMES fftw3f_omp HINTS "${HAVE_FFTWF_OMP_LIBRARY}" )
		if (HAVE_FFTWF_THREADS_LIBRARY OR HAVE_FFTWF_OMP_LIBRARY)
			if (HAVE_FFTWF_THREADS_LIBRARY)
				set(FFT_LIBRARY fftw3f_threads fftw3f)
				link_directories(${HAVE_FFTWF_THREADS_LIBRARY})
			else (HAVE_FFTWF_THREADS_LIBRARY)
				set(FFT_LIBRARY fftw3f_omp fftw3f)
				link_directories(${HAVE_FFTWF_OMP_LIBRARY})
			endif (HAVE_FFTWF_THREADS_LIBRARY)
			add_definitions(-DHAVE_FFTW_THREADS)
		else (HAVE_FFTWF_THREADS_LIBRARY OR HAVE_FFTWF_OMP_LIBRARY)
			message(WARNING "Cannot find fftw3_threads or fftw3_omp, will use single-threaded version")
			CHECK_LIBRARY_EXISTS(fftw3f fftwf_plan_many_r2r "" HAVE_FFTWF_LIBRARY)
			if (NOT HAVE_FFTWF_LIBRARY)
				message(FATAL_ERROR "Cannot find fftw3f library")
			endif (NOT HAVE_FFTWF_LIBRARY)
			set(FFT_LIBRARY fftw3f)
			link_directories(${HAVE_FFTWF_LIBRARY})
		endif (HAVE_FFTWF_THREADS_LIBRARY OR HAVE_FFTWF_OMP_LIBRARY)
	endif (HAVE_SINGLE)
	if (HAVE_DOUBLE)
		find_library(HAVE_FFTW_THREADS_LIBRARY fftw3_threads)
		find_library(HAVE_FFTW_OMP_LIBRARY fftw3_omp)
		if (HAVE_FFTW_THREADS_LIBRARY OR HAVE_FFTW_OMP_LIBRARY)
			if (HAVE_FFTW_THREADS_LIBRARY)
				set(FFT_LIBRARY fftw3_threads fftw3)
				link_directories(${HAVE_FFTW_THREADS_LIBRARY})
			else (HAVE_FFTW_THREADS_LIBRARY)
				set(FFT_LIBRARY fftw3_omp fftw3)
				link_directories(${HAVE_FFTW_OMP_LIBRARY})
			endif (HAVE_FFTW_THREADS_LIBRARY)
			add_definitions(-DHAVE_FFTW_THREADS)
		else (HAVE_FFTW_THREADS_LIBRARY OR HAVE_FFTW_OMP_LIBRARY)
			message(WARNING "Cannot find fftw3_threads or fftw3_omp, will use single-threaded version")
			CHECK_LIBRARY_EXISTS(fftw3 fftw_plan_many_r2r "" HAVE_FFTW_LIBRARY)
			if (NOT HAVE_FFTW_LIBRARY)
				message(FATAL_ERROR "Cannot find fftw3 library")
			endif (NOT HAVE_FFTW_LIBRARY)
			set(FFT_LIBRARY fftw3)
			link_directories(${HAVE_FFTWF_LIBRARY})
		endif (HAVE_FFTW_THREADS_LIBRARY OR HAVE_FFTW_OMP_LIBRARY)
	endif (HAVE_DOUBLE)
	add_definitions(-DHAVE_FFTW)
endif (HAVE_FFTW)

if (HAVE_FFTW_MKL)
	CHECK_INCLUDE_FILE(fftw3_mkl.h HAVE_FFTW_MKL_HEADER)
	if (NOT HAVE_FFTW_MKL_HEADER)
		message(FATAL_ERROR "Cannot find fftw3_mkl.h")
	endif (NOT HAVE_FFTW_MKL_HEADER)
	include_directories(${HAVE_FFTW_MKL_HEADER})
	find_library(HAVE_FFTW_MKL_THREADS_LIBRARY NAMES mkl_intel_thread HINTS "${HAVE_FFTW_MKL_THREADS_LIBRARY}")
	if (HAVE_FFTW_MKL_THREADS_LIBRARY)
		set(FFT_LIBRARY mkl_intel_thread mkl_rt)
		link_directories(${HAVE_FFTW_MKL_THREADS_LIBRARY})
		add_definitions(-DHAVE_FFTW_THREADS)
	else (HAVE_FFTW_MKL_THREADS_LIBRARY)
		message(WARNING "Cannot find mkl_intel_thread, will use single-threaded version")
		find_library(HAVE_FFTW_MKL_LIBRARY NAMES mkl_sequential HINTS "${HAVE_FFTW_MKL_LIBRARY}")
		if (NOT HAVE_FFTW_MKL_LIBRARY)
			message(FATAL_ERROR "Cannot find mkl library")
		endif (NOT HAVE_FFTW_MKL_LIBRARY)
		set(FFT_LIBRARY mkl_sequential mkl_rt)
		link_directories(${HAVE_FFTW_MKL_LIBRARY})
	endif (HAVE_FFTW_MKL_THREADS_LIBRARY)
	add_definitions(-DHAVE_FFTW_MKL)
endif (HAVE_FFTW_MKL)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/tests/poisson2d_fft")
include_directories("${PROJECT_SOURCE_DIR}/tests/poisson2d_fft")

install(FILES breeze2d.h DESTINATION include)
install(FILES breeze2d_advection.h DESTINATION include)
install(FILES breeze2d_interop.h DESTINATION include)
install(FILES breeze2d_poisson.h DESTINATION include)
install(FILES breeze2d_status.h DESTINATION include)
install(FILES breeze2d_timing.h DESTINATION include)

add_library(poisson2d
	poisson2d/poisson2d.c
	poisson2d/fft/fft.c poisson2d/fft/fft.h
	poisson2d/fft/shutter.h poisson2d/fft/shutter_c.cpp poisson2d/fft/shutter_r.c
	poisson2d/fft/wrapper.c poisson2d/fft/wrapper.h)
add_subdirectory(poisson2d)

add_library(interop
	interop/dump2db.c interop/grads.c)
add_subdirectory(interop)

add_library(timing
	timing/timing.c)

add_executable(poisson2d_fft tests/poisson2d_fft/poisson2d_fft.c)
target_link_libraries(poisson2d_fft
	poisson2d interop timing lapack ${FFT_LIBRARY})
install(TARGETS poisson2d_fft DESTINATION bin)

file(COPY ${PROJECT_SOURCE_DIR}/cbarm.gs DESTINATION ${PROJECT_BINARY_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -O0")

