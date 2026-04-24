cmake_minimum_required(VERSION 3.18)
project(fesom2_port C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
)
target_include_directories(fesom_port PRIVATE src ${NC_INCDIR})
target_link_libraries(fesom_port PRIVATE m ${NC_LIBS})
