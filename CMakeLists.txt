cmake_minimum_required(VERSION 3.18)
project(fesom2_port C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Phase 4: MPI is a hard dependency.
find_package(MPI REQUIRED COMPONENTS C)

# Detect netCDF via nc-config. On Levante: `module load netcdf-c/4.8.1-gcc-11.2.0`
# puts the right nc-config first in PATH.
find_program(NC_CONFIG nc-config)
if(NOT NC_CONFIG)
  message(FATAL_ERROR "nc-config not found. On Levante: module load netcdf-c/4.8.1-gcc-11.2.0")
endif()
execute_process(COMMAND ${NC_CONFIG} --includedir
                OUTPUT_VARIABLE NC_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${NC_CONFIG} --libs
                OUTPUT_VARIABLE NC_LIBS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(NC_LIBS UNIX_COMMAND "${NC_LIBS_RAW}")
execute_process(COMMAND ${NC_CONFIG} --libdir
                OUTPUT_VARIABLE NC_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "netCDF include: ${NC_INCDIR}")
message(STATUS "netCDF libdir:  ${NC_LIBDIR}")
message(STATUS "netCDF link:    ${NC_LIBS}")
# Embed the netCDF lib path into the binary so it runs without `module load`.
set(CMAKE_INSTALL_RPATH "${NC_LIBDIR}")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Warnings tight, frame pointer kept for profiling.
add_compile_options(
  -Wall -Wextra -Wshadow -Wstrict-prototypes -Wmissing-prototypes
  -Wno-unused-parameter
  -fno-omit-frame-pointer
)
add_compile_options(
  "$<$<CONFIG:Debug>:-O0>"  "$<$<CONFIG:Debug>:-g3>"
  "$<$<CONFIG:Release>:-O3>" "$<$<CONFIG:Release>:-g>"
)

add_executable(fesom_port
  src/fesom_main.c
  src/fesom_mesh.c
  src/fesom_mpi.c
  src/fesom_dyn.c
  src/fesom_tracers.c
  src/fesom_aux.c
  src/fesom_forcing.c
  src/fesom_ic.c
  src/fesom_eos.c
  src/fesom_ale.c
  src/fesom_ssh.c
  src/fesom_momentum.c
  src/fesom_tracer_adv.c
  src/fesom_pp.c
  src/fesom_step.c
  src/fesom_forcing_analytical.c
  src/fesom_io.c
  src/fesom_constants.c
  src/fesom_phc.c
  src/fesom_jra55.c
  src/fesom_bulk.c
  src/fesom_tracer_diff.c
  src/fesom_sss_runoff.c
  src/fesom_partit.c
  src/fesom_halo.c
  src/fesom_ice.c
  src/fesom_ice_thermo.c
  src/fesom_ice_coupling.c
  src/fesom_ice_evp.c
  src/fesom_ice_fct.c
  src/fesom_gm.c
  src/fesom_kpp.c
  src/fesom_calendar.c
  src/fesom_io_stream.c
  src/fesom_io_stream_dispatch.c
  src/fesom_io_config.c
)
target_include_directories(fesom_port PRIVATE src ${NC_INCDIR})
target_link_libraries(fesom_port PRIVATE m ${NC_LIBS} MPI::MPI_C)

# --- Tests --------------------------------------------------------------
# Calendar is pure logic — no MPI, no netCDF — so test_calendar links only
# fesom_calendar.c and the math library. CTest invocation: `ctest -V`.
enable_testing()
add_executable(test_calendar
  tests/test_calendar.c
  src/fesom_calendar.c
)
target_include_directories(test_calendar PRIVATE src)
target_link_libraries(test_calendar PRIVATE m)
add_test(NAME calendar COMMAND test_calendar)

# Stream dispatcher unit test: no netCDF, no fesom_io.c, just MPI for
# the FESOM_DIE → MPI_Abort symbol on the never-taken error path.
add_executable(test_io_stream_unit
  tests/test_io_stream_unit.c
  src/fesom_io_stream_dispatch.c
)
target_include_directories(test_io_stream_unit PRIVATE src)
target_link_libraries(test_io_stream_unit PRIVATE m MPI::MPI_C)
add_test(NAME io_stream_unit COMMAND test_io_stream_unit)

# Config-file parser test: pure C, no MPI/netCDF.
add_executable(test_io_config
  tests/test_io_config.c
  src/fesom_io_config.c
)
target_include_directories(test_io_config PRIVATE src)
target_compile_definitions(test_io_config PRIVATE _POSIX_C_SOURCE=200809L)
target_link_libraries(test_io_config PRIVATE m)
add_test(NAME io_config COMMAND test_io_config)
