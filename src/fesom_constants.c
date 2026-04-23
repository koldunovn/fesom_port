#include "fesom_constants.h"

/* Definition (sole) of the runtime-mutable timestep. Default = pi step_per_day=36
   → dt = 86400/36 = 2400 s. Override from CLI in fesom_main.c. */
real_t fesom_phase1_dt = 2400.0;
