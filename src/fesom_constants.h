#ifndef FESOM_CONSTANTS_H
#define FESOM_CONSTANTS_H

/* Physical constants — FRESH_START.md §17. */
#define FESOM_PI         3.14159265358979
#define FESOM_RAD        (FESOM_PI / 180.0)
#define FESOM_DENSITY_0  1030.0                    /* reference density, kg/m³ */
#define FESOM_G          9.81                      /* m/s² */
#define FESOM_R_EARTH    6367500.0                 /* m */
#define FESOM_OMEGA      (2.0 * FESOM_PI / 86400.0) /* earth rotation, rad/s */
#define FESOM_VCPW       4.2e6                     /* J/(m³·K) */

/*
 * Mesh rotation. FESOM2 default: pole at (lon=50°, lat=15°), rotation
 * gamma=-90°. work_pi/namelist.config and work_core/namelist.config
 * both set these defaults. Almost all computations are performed in
 * rotated coordinates — see FRESH_START.md §2.
 */
#define FESOM_ALPHA_EULER_DEG   50.0
#define FESOM_BETA_EULER_DEG    15.0
#define FESOM_GAMMA_EULER_DEG  (-90.0)
#define FESOM_FORCE_ROTATION    1   /* 1 = apply g2r at mesh load */

/*
 * Cyclic domain length in rotated radians. namelist.config sets
 * cyclic_length=360 [degrees], then gen_model_setup.F90:97 multiplies
 * by `rad` so internal value is 2π. Mirror that here.
 */
#define FESOM_CYCLIC_LENGTH_RAD  (2.0 * FESOM_PI)

/*
 * Hardcoded Phase 1 namelist defaults, derived from FRESH_START.md §14.7
 * with the Phase 1 simplifications:
 *   - which_ALE        = 'linfs' (simplest free surface)
 *   - mix_scheme       = 'PP'    (simpler than KPP)
 *   - GM/Redi          = .false. (off until Phase 4)
 *
 * No namelist parser yet — these are intentionally compile-time
 * constants until we have a working Phase 1 run.
 */
#define FESOM_PHASE1_K_VER     1.0e-5   /* background tracer vertical diff [m²/s]   */
#define FESOM_PHASE1_A_VER     1.0e-4   /* background momentum vertical visc [m²/s] */
#define FESOM_PHASE1_C_D       0.0025   /* bottom-drag coefficient                  */

/* PP mixing — defaults from oce_modules.F90:25-78 */
#define FESOM_PHASE1_MIX_COEFF_PP   0.01    /* PP scaling coefficient                  */
#define FESOM_PHASE1_INSTABMIX_KV   0.1     /* convective-adjustment Kv [m²/s]         */
/* Phase 1 disables Monin-Obukhov mixing (Fortran default is .true. but it
   needs sea-ice + forcing arrays we don't carry yet; FRESH_START.md §10
   describes only PP + convective, not momix). */
#define FESOM_PHASE1_USE_MOMIX      0
#define FESOM_PHASE1_USE_INSTABMIX  1
#define FESOM_PHASE1_USE_WINDMIX    0       /* matches Fortran default */
#define FESOM_PHASE1_AB_ORDER  2        /* Adams-Bashforth order for Coriolis       */

/* SSH solver implicitness — Fortran defaults from oce_modules.F90:90 are
   alpha=1.0 and theta=1.0. Mirror those literally. */
#define FESOM_PHASE1_ALPHA     1.0
#define FESOM_PHASE1_THETA     1.0

/* Timestep — pi run uses step_per_day=36 → dt = 86400/36 = 2400 s (per
   work_pi/namelist.config). Phase 1 keeps it as a compile-time constant. */
#define FESOM_PHASE1_DT        2400.0

/* CG SSH solver defaults — MOD_DYN.F90:13-25. */
#define FESOM_PHASE1_SOLTOL    1.0e-5
#define FESOM_PHASE1_MAXITER   2000

#endif /* FESOM_CONSTANTS_H */
