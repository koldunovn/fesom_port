cmake_minimum_required(VERSION 3.18)
project(fesom2_port C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
)
target_include_directories(fesom_port PRIVATE src)
target_link_libraries(fesom_port PRIVATE m)
