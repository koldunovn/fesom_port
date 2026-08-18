/*
 * fesom_phc.c — PHC3.0 (and similar) 3D initial-condition reader.
 * Literal port of gen_ic3d.F90 + extrap_nod3D + insitu2pot. See header for
 * source-line references.
 */
#include "fesom_phc.h"
#include "fesom_constants.h"
#include "fesom_halo.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"
#include "fesom_tracers.h"

#include <math.h>
#include <mpi.h>
#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NC_CHECK(call) do {                                             \
    int _ec = (call);                                                   \
    if (_ec != NC_NOERR) {                                              \
        FESOM_DIE("NetCDF error: %s", nc_strerror(_ec));                \
    }                                                                   \
} while (0)

/* g_config dummy value (gen_modules_config.F90:174). Land marker. */
#define PHC_DUMMY  1.0e10

/*--- binarysearch (gen_ic3d.F90:666-703) ---------------------------------
 * Returns 0-based index of element <= value. Returns -1 if value < arr[0]
 * (Fortran returns 0 — we use -1 in C for "not found"). Matches Fortran's
 * `right` final value when value is in-range.
 *
 * Note: Fortran's binarysearch does an exact-match short-circuit with a
 * 1e-9 tolerance, then returns `right` after the loop. Mirror it.
 */
static int binarysearch_d(int length, const double *arr, double value)
{
    const double d = 1e-9;
    int left = 0;
    int right = length - 1;
    while (left <= right) {
        int middle = (left + right + 1) / 2;
        if (fabs(arr[middle] - value) <= d) return middle;
        if (arr[middle] > value) right = middle - 1;
        else                     left  = middle + 1;
    }
    return right;
}

/*--- atg + ptheta (Bryden 1973) — oce_ale_pressure_bv.F90:2674-2746 -----
 * adiabatic temperature gradient and potential-T integration via 4th-order
 * Runge-Kutta. Identical numerical recipe; verified against Fortran.
 */
static double atg(double s, double t, double p)
{
    double ds = s - 35.0;
    return (((-2.1687e-16*t + 1.8676e-14)*t - 4.6206e-13)*p
           +((2.7759e-12*t - 1.1351e-10)*ds + ((-5.4481e-14*t
           + 8.733e-12)*t - 6.7795e-10)*t + 1.8741e-8))*p
           +(-4.2393e-8*t + 1.8932e-6)*ds
           +((6.6228e-10*t - 6.836e-8)*t + 8.5258e-6)*t + 3.5803e-5;
}

static double ptheta(double s, double t, double p, double pr)
{
    double h, xk, q;
    h = pr - p;
    xk = h * atg(s, t, p);
    t  = t + 0.5 * xk;
    q  = xk;
    p  = p + 0.5 * h;
    xk = h * atg(s, t, p);
    t  = t + 0.29289322 * (xk - q);
    q  = 0.58578644 * xk + 0.121320344 * q;
    xk = h * atg(s, t, p);
    t  = t + 1.707106781 * (xk - q);
    q  = 3.414213562 * xk - 4.121320344 * q;
    p  = p + 0.5 * h;
    xk = h * atg(s, t, p);
    return t + (xk - 2.0 * q) / 6.0;
}

/*--- One-variable load + bilinear-vertical interpolation -----------------
 * Mirror of gen_ic3d.F90:303-491 (getcoeffld) restricted to non-cavity.
 *
 *   Args:
 *     ncid:        already-open NetCDF file
 *     varname:     "temp" or "salt"
 *     mesh:        target FESOM mesh
 *     nc_lon, nc_lat, nc_depth, Nlon, Nlat, Ndepth:  PHC grid (lon already cyclic-wrapped)
 *     bilin_i, bilin_j: per-node bracket indices from nc_ic3d_ini
 *     fill_value:  optional from file's _FillValue (else NaN test)
 *     out:         tracers->data[*].values  [N * nl] — written
 */
static void load_one_variable(int ncid, const char *varname,
                              const struct fesom_mesh *mesh,
                              const double *nc_lon, const double *nc_lat, const double *nc_depth,
                              int Nlon, int Nlat, int Ndepth,
                              const int *bilin_i, const int *bilin_j,
                              double fill_value, int has_fill,
                              real_t *out)
{
    /* Allocate full ncdata[Nlon][Nlat][Ndepth] (C row-major; first index Nlon). */
    size_t ntot = (size_t)Nlon * (size_t)Nlat * (size_t)Ndepth;
    double *ncdata = malloc(ntot * sizeof(double));
    FESOM_CHECK(ncdata, "phc: oom (ncdata for %s)", varname);

    /* PHC stores temp/salt as (depth, lat, lon). nf_get_vara reads in dimension
       order. We want ncdata[ilon, ilat, idepth] for cache locality on the lon
       axis (since Fortran reads lon as fastest-varying). To match Fortran
       semantics on the wrap (interior = indices [1..Nlon-2], wrap from 1
       and Nlon-2), we read into a temporary (depth, lat, Nlon-2) buffer then
       transpose + wrap. */
    size_t nint = (size_t)(Nlon - 2);
    double *tmp = malloc((size_t)Ndepth * (size_t)Nlat * nint * sizeof(double));
    FESOM_CHECK(tmp, "phc: oom (tmp for %s)", varname);

    int varid;
    NC_CHECK(nc_inq_varid(ncid, varname, &varid));

    /* Nlon is the PADDED count (real+2); read all `nint = Nlon-2` real columns
       from file lon index 0 into the padded interior ncdata[1..Nlon-2]. The two
       halo columns [0] and [Nlon-1] are filled by the wrap copies below. */
    size_t start[3] = { 0, 0, 0 };          /* depth, lat, lon */
    size_t count[3] = { (size_t)Ndepth, (size_t)Nlat, nint };
    NC_CHECK(nc_get_vara_double(ncid, varid, start, count, tmp));

    /* Transpose to ncdata[ilon, ilat, idepth] AND apply wrap.
       After transpose:
         ncdata at lon-index 1..Nlon-2 = file data at lon-index 1..Nlon-2
         ncdata at lon-index 0          = ncdata at lon-index Nlon-2  (Fortran: ncdata(1) = ncdata(nc_Nlon-1))
         ncdata at lon-index Nlon-1     = ncdata at lon-index 1       (Fortran: ncdata(nc_Nlon) = ncdata(2))
       Replace NaN/FillValue with PHC_DUMMY in the same pass. */
    for (int ilon_int = 0; ilon_int < (int)nint; ++ilon_int) {
        int ilon = ilon_int + 1;            /* 0-based real-data index */
        for (int ilat = 0; ilat < Nlat; ++ilat) {
            for (int idep = 0; idep < Ndepth; ++idep) {
                /* tmp index: (depth, lat, lon) row-major in C with stride
                   (Nlat * nint, nint, 1). */
                size_t tidx = ((size_t)idep * (size_t)Nlat + (size_t)ilat) * nint + (size_t)ilon_int;
                double v = tmp[tidx];
                if (isnan(v)) v = PHC_DUMMY;
                if (has_fill && v == fill_value) v = PHC_DUMMY;
                if (v < -0.99 * PHC_DUMMY || v > PHC_DUMMY) v = PHC_DUMMY;
                size_t didx = ((size_t)ilon * (size_t)Nlat + (size_t)ilat) * (size_t)Ndepth + (size_t)idep;
                ncdata[didx] = v;
            }
        }
    }
    /* Wrap copies (Fortran lines 370-371). */
    for (int ilat = 0; ilat < Nlat; ++ilat) {
        for (int idep = 0; idep < Ndepth; ++idep) {
            size_t src_w = ((size_t)(Nlon - 2) * (size_t)Nlat + (size_t)ilat) * (size_t)Ndepth + (size_t)idep;
            size_t dst_w = ((size_t)0           * (size_t)Nlat + (size_t)ilat) * (size_t)Ndepth + (size_t)idep;
            ncdata[dst_w] = ncdata[src_w];
            size_t src_e = ((size_t)1          * (size_t)Nlat + (size_t)ilat) * (size_t)Ndepth + (size_t)idep;
            size_t dst_e = ((size_t)(Nlon - 1) * (size_t)Nlat + (size_t)ilat) * (size_t)Ndepth + (size_t)idep;
            ncdata[dst_e] = ncdata[src_e];
        }
    }
    free(tmp);

    const int N  = mesh->myDim_nod2D;
    const int nl = mesh->nl;

    /* Initialise out to dummy (Fortran line 342). */
    for (int n = 0; n < N; ++n) {
        for (int nz = 0; nz < nl; ++nz) {
            out[FESOM_NODE3D(n, nz, nl)] = PHC_DUMMY;
        }
    }

    /* Per-node bilinear horizontal + linear vertical interpolation
       (Fortran lines 392-485). */
    double *data1d = malloc((size_t)Ndepth * sizeof(double));
    FESOM_CHECK(data1d, "phc: oom (data1d)");

    for (int ii = 0; ii < N; ++ii) {
        int nl1 = mesh->nlevels_nod2D[ii] - 1;   /* number of valid layers   */
        int ul1 = mesh->ulevels_nod2D[ii] - 1;   /* topmost layer (0 typically) */
        int i = bilin_i[ii];
        int j = bilin_j[ii];
        if (i < 0 || j < 0) continue;             /* out of forcing domain    */
        int ip1 = i + 1, jp1 = j + 1;

        /* Skip if the surface (depth index 0) of any of the 4 corners is
           dummy — Fortran line 415 (any(ncdata(i:ip1, j:jp1, 1) > 0.99*dummy)). */
        double s00 = ncdata[((size_t)i  *(size_t)Nlat + (size_t)j  )*(size_t)Ndepth];
        double s10 = ncdata[((size_t)ip1*(size_t)Nlat + (size_t)j  )*(size_t)Ndepth];
        double s01 = ncdata[((size_t)i  *(size_t)Nlat + (size_t)jp1)*(size_t)Ndepth];
        double s11 = ncdata[((size_t)ip1*(size_t)Nlat + (size_t)jp1)*(size_t)Ndepth];
        if (s00 > 0.99*PHC_DUMMY || s10 > 0.99*PHC_DUMMY ||
            s01 > 0.99*PHC_DUMMY || s11 > 0.99*PHC_DUMMY) continue;

        /* Geographic node coords for interp weights. */
        double x = (double)mesh->geo_coord_nod2D[2*ii + 0] / FESOM_RAD;
        double y = (double)mesh->geo_coord_nod2D[2*ii + 1] / FESOM_RAD;
        if (x < 0.0)   x += 360.0;
        if (x > 360.0) x -= 360.0;

        double x1 = nc_lon[i];
        double x2 = nc_lon[ip1];
        double y1 = nc_lat[j];
        double y2 = nc_lat[jp1];
        double denom = (x2 - x1) * (y2 - y1);

        /* Bilinear interp per depth — entire column (Fortran lines 423-428). */
        for (int k = 0; k < Ndepth; ++k) {
            double v00 = ncdata[((size_t)i  *(size_t)Nlat + (size_t)j  )*(size_t)Ndepth + k];
            double v10 = ncdata[((size_t)ip1*(size_t)Nlat + (size_t)j  )*(size_t)Ndepth + k];
            double v01 = ncdata[((size_t)i  *(size_t)Nlat + (size_t)jp1)*(size_t)Ndepth + k];
            double v11 = ncdata[((size_t)ip1*(size_t)Nlat + (size_t)jp1)*(size_t)Ndepth + k];
            data1d[k] = (v00 * (x2 - x) * (y2 - y)
                       + v10 * (x  - x1) * (y2 - y)
                       + v01 * (x2 - x) * (y  - y1)
                       + v11 * (x  - x1) * (y  - y1)) / denom;
            if (v00 > 0.99*PHC_DUMMY || v10 > 0.99*PHC_DUMMY ||
                v01 > 0.99*PHC_DUMMY || v11 > 0.99*PHC_DUMMY) {
                data1d[k] = PHC_DUMMY;
            }
        }

        /* Vertical interpolation onto FESOM Z[k] (lines 463-482, non-cavity). */
        for (int k = ul1; k < nl1; ++k) {
            /* Fortran: -Z_3d_n[k, ii] (positive depth in metres). */
            double depth_pos = -(double)mesh->Z[k];
            int d_indx = binarysearch_d(Ndepth, nc_depth, depth_pos);
            if (d_indx >= 0 && d_indx < Ndepth - 1) {
                int d_indx_p1 = d_indx + 1;
                double delta_d = nc_depth[d_indx + 1] - nc_depth[d_indx];
                double d1 = data1d[d_indx];
                double d2 = data1d[d_indx_p1];
                if (d1 < 0.99*PHC_DUMMY && d2 < 0.99*PHC_DUMMY) {
                    double cf_a = (d2 - d1) / delta_d;
                    double cf_b = d1 - cf_a * nc_depth[d_indx];
                    /* Fortran: -cf_a * Z_3d_n[k, ii] + cf_b */
                    out[FESOM_NODE3D(ii, k, nl)] =
                        (real_t)(-cf_a * (double)mesh->Z[k] + cf_b);
                }
            } else if (d_indx == -1) {
                /* Above the topmost PHC depth — use surface value. */
                if (data1d[0] < 0.99*PHC_DUMMY) {
                    out[FESOM_NODE3D(ii, k, nl)] = (real_t)data1d[0];
                }
            }
            /* d_indx == Ndepth-1 (below deepest PHC) → stays dummy; extrap_nod
               will fill it from above. */
        }
    }

    free(data1d);
    free(ncdata);
}

/*--- extrap_nod3D (gen_support.F90:400-507) ---------------
 * Literal port. MPI variant: each rank fills its own myDim nodes, then
 * exchange_nod is called between sweeps so other ranks' progress propagates
 * across partition boundaries. The outer-loop continuation test is a
 * GLOBAL max via MPI_Allreduce — without this, ranks stop independently
 * and a coastal node whose only valid PHC source lives on another rank
 * keeps its dummy → cleanup converts to 0 → wrong density → eta blow-up. */
static void extrap_nod3D(const struct fesom_mesh *mesh,
                         struct fesom_partit     *partit,
                         real_t                  *arr)
{
    const int N       = mesh->myDim_nod2D;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl      = mesh->nl;

    /* work_array sized for full local extent — Fortran line 418. Reads at
     * `work[enodes(j)]` may index halo nodes (eDim) which need valid values. */
    real_t *work = malloc((size_t)N_alloc * sizeof(real_t));
    FESOM_CHECK(work, "extrap_nod3D: oom");

    /* Initial halo exchange so halo entries of arr reflect owner state.
     * Fortran line 419: `call exchange_nod(arr, partit)` before the outer loop. */
    if (partit && partit->npes > 1) {
        for (int nz = 0; nz < nl; ++nz) {
            /* Pack a layer slice, exchange, unpack. arr is laid out as
             * [node][nz] in row-major C; we'd prefer a strided exchange
             * but our halo module operates on contiguous per-entity buffers.
             * Here we exchange the whole 3D array layer-by-layer using a
             * temporary scalar field, which is cheap (PHC IC is one-shot). */
            real_t *layer = malloc((size_t)N_alloc * sizeof(real_t));
            for (int n = 0; n < N_alloc; ++n)
                layer[n] = arr[FESOM_NODE3D(n, nz, nl)];
            fesom_halo_exchange(layer, FESOM_HALO_NOD2D, 1, 1, partit);
            for (int n = 0; n < N_alloc; ++n)
                arr[FESOM_NODE3D(n, nz, nl)] = layer[n];
            free(layer);
        }
    }

    int iter_outer = 0;
    while (iter_outer < 200) {
        ++iter_outer;
        real_t loc_max = arr[FESOM_NODE3D(0, 0, nl)];
        for (int n = 1; n < N; ++n) {
            real_t v = arr[FESOM_NODE3D(n, 0, nl)];
            if (v > loc_max) loc_max = v;
        }
        real_t glob_max = loc_max;
        if (partit && partit->npes > 1) {
            MPI_Allreduce(&loc_max, &glob_max, 1, MPI_DOUBLE, MPI_MAX,
                          partit->MPI_COMM_FESOM);
        }
        if (glob_max <= 0.99 * PHC_DUMMY) break;

        for (int nz = 0; nz < nl - 1; ++nz) {
            for (int n = 0; n < N_alloc; ++n) work[n] = arr[FESOM_NODE3D(n, nz, nl)];
            int success = 1;
            int sweep = 0;
            while (success && sweep < 200) {
                success = 0;
                ++sweep;
                for (int n = 0; n < N; ++n) {
                    if (!(work[n] > 0.99 * PHC_DUMMY) ||
                        !(mesh->nlevels_nod2D[n] - 1 > nz)) continue;
                    real_t val = 0.0;
                    int cnt = 0;
                    int o0 = mesh->nod_in_elem2D_offsets[n];
                    int o1 = mesh->nod_in_elem2D_offsets[n + 1];
                    for (int kk = o0; kk < o1; ++kk) {
                        int el = mesh->nod_in_elem2D[kk];
                        if (nz > mesh->nlevels[el] - 1) continue;
                        for (int j = 0; j < 3; ++j) {
                            int v = mesh->elem_nodes[3*el + j];
                            if (v < 0 || v >= N_alloc) continue;
                            if (work[v] < 0.99 * PHC_DUMMY &&
                                mesh->nlevels_nod2D[v] - 1 > nz) {
                                val += work[v];
                                ++cnt;
                            }
                        }
                    }
                    if (cnt > 0) {
                        work[n] = val / (real_t)cnt;
                        success = 1;
                    }
                }
            }
            for (int n = 0; n < N; ++n) arr[FESOM_NODE3D(n, nz, nl)] = work[n];
        }

        /* Cross-rank propagation between outer iterations — Fortran line 484. */
        if (partit && partit->npes > 1) {
            for (int nz = 0; nz < nl; ++nz) {
                real_t *layer = malloc((size_t)N_alloc * sizeof(real_t));
                for (int n = 0; n < N_alloc; ++n)
                    layer[n] = arr[FESOM_NODE3D(n, nz, nl)];
                fesom_halo_exchange(layer, FESOM_HALO_NOD2D, 1, 1, partit);
                for (int n = 0; n < N_alloc; ++n)
                    arr[FESOM_NODE3D(n, nz, nl)] = layer[n];
                free(layer);
            }
        }
    }

    /* Vertical fill — Fortran lines 491-498. */
    for (int n = 0; n < N; ++n) {
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = 1; nz < nl1; ++nz) {
            if (arr[FESOM_NODE3D(n, nz, nl)] > 0.99 * PHC_DUMMY) {
                arr[FESOM_NODE3D(n, nz, nl)] = arr[FESOM_NODE3D(n, nz - 1, nl)];
            }
        }
    }
    /* Final exchange — Fortran line 502. */
    if (partit && partit->npes > 1) {
        for (int nz = 0; nz < nl; ++nz) {
            real_t *layer = malloc((size_t)N_alloc * sizeof(real_t));
            for (int n = 0; n < N_alloc; ++n)
                layer[n] = arr[FESOM_NODE3D(n, nz, nl)];
            fesom_halo_exchange(layer, FESOM_HALO_NOD2D, 1, 1, partit);
            for (int n = 0; n < N_alloc; ++n)
                arr[FESOM_NODE3D(n, nz, nl)] = layer[n];
            free(layer);
        }
    }

    free(work);
}

/*--------------------------------------------------------------------------
 * BACK-PORT, and against the grain of the project's direction.
 *
 * Everything else in this tree is a literal Fortran -> C port and the Kokkos
 * tree is a port OF THIS ONE. This block runs the other way: extrap_nod3D_det
 * is an M13 invention with no Fortran original, written in the Kokkos tree to
 * fix a defect that is present upstream in gen_support.F90 and was inherited
 * here. It is copied back so the C reference can produce the same initial
 * condition, and so the C-vs-Kokkos twin can be gated byte for byte.
 *
 * Kept verbatim from port_kokkos_int/src/fesom_phc.cpp on purpose: byte
 * identity of the twin depends on the tolerance, the iteration cap AND the
 * summation order, so the two copies must not drift. There is nothing
 * Kokkos-specific in it — it is host code over raw pointers — so "verbatim" is
 * achievable rather than aspirational. tests/ic_extrap_twin/check_det_drift.py
 * fails if the two ever diverge.
 *
 * The knob defaults to legacy, so every byte-certified baseline in this tree is
 * untouched.
 *------------------------------------------------------------------------*/
/*--- extrap_nod3D_det — M13 CG-blowups fix (FESOM_IC_EXTRAP=det) ----------
 * The legacy extrap_nod3D is ORDER- AND PARTITION-DEPENDENT: it fills each
 * dummy node from already-valid neighbours by in-place Gauss-Seidel, runs each
 * rank's fill to local exhaustion between halo exchanges, and never revisits a
 * filled node (first-fill-wins). Marginal seas (e.g. the whole Sea of Marmara:
 * any land corner of the 1° PHC cell marks a node dummy) therefore get filled
 * with DIFFERENT WATER MASSES on different partitions — measured 2026-08-14 on
 * NG5 elem 8241930 nz=36: S = 18.611/18.611/18.611 on dist_4096 vs
 * 21.9/25.5/27.9 on dist_8192 → Δrho ≈ 10 kg/m³ across one 2-km element →
 * pgf 5.5e-3 m/s² → linear |uv| ramp 0.5 m/s per step from step 1 → CFL blow-up
 * (the "CG blow-ups on specific partitions" class; same defect upstream in
 * gen_support.F90:400-507).
 *
 * det variant: ring-Jacobi. Per sweep, validity is FROZEN to the previous
 * sweep's state; every still-dummy node with ≥1 valid neighbour is filled with
 * the multiplicity-weighted mean of those neighbours (same stencil and
 * eligibility tests as legacy), accumulated in ASCENDING GLOBAL-ID order so
 * the FP sum is bitwise partition-independent; halo exchange after every
 * sweep; terminate when a global Allreduce reports no progress on any layer.
 * The fill thus depends only on the global mesh + ring distance, not on the
 * decomposition or the local numbering. NOTE on np=1 this differs from legacy
 * (Jacobi ≠ Gauss-Seidel), so the knob defaults OFF and byte-certified
 * baselines are untouched.
 */
static void extrap_det_exchange_all(const struct fesom_mesh *mesh,
                                    struct fesom_partit *partit, real_t *arr)
{
    if (!partit || partit->npes <= 1) return;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl = mesh->nl;
    real_t *layer = (real_t *)malloc((size_t)N_alloc * sizeof(real_t));
    FESOM_CHECK(layer, "extrap_det: oom (layer)");
    for (int nz = 0; nz < nl; ++nz) {
        for (int n = 0; n < N_alloc; ++n)
            layer[n] = arr[FESOM_NODE3D(n, nz, nl)];
        fesom_halo_exchange(layer, FESOM_HALO_NOD2D, 1, 1, partit);
        for (int n = 0; n < N_alloc; ++n)
            arr[FESOM_NODE3D(n, nz, nl)] = layer[n];
    }
    free(layer);
}

/* Layer-local mean of a node's eligible non-dummy neighbours, accumulated in
 * ascending-gid order (partition-independent FP sum). Returns count. */
static int extrap_det_ring_mean(const struct fesom_mesh *mesh,
                                struct fesom_partit *partit,
                                const real_t *layer, int n, int nz,
                                int N_alloc, real_t *out)
{
    int    ngid[64];
    real_t nval[64];
    int    cnt = 0;
    int o0 = mesh->nod_in_elem2D_offsets[n];
    int o1 = mesh->nod_in_elem2D_offsets[n + 1];
    for (int kk = o0; kk < o1; ++kk) {
        int el = mesh->nod_in_elem2D[kk];
        if (nz > mesh->nlevels[el] - 1) continue;
        for (int j = 0; j < 3; ++j) {
            int v = mesh->elem_nodes[3*el + j];
            if (v < 0 || v >= N_alloc || v == n) continue;
            if (layer[v] < 0.99 * PHC_DUMMY &&
                mesh->nlevels_nod2D[v] - 1 > nz) {
                FESOM_CHECK(cnt < 64, "extrap_det: ring overflow");
                ngid[cnt] = partit->myList_nod2D ? partit->myList_nod2D[v] : v + 1;
                nval[cnt] = layer[v];
                ++cnt;
            }
        }
    }
    if (cnt == 0) return 0;
    for (int a = 1; a < cnt; ++a) {           /* insertion sort by gid */
        int    g = ngid[a]; real_t w = nval[a]; int b = a - 1;
        while (b >= 0 && ngid[b] > g) {
            ngid[b+1] = ngid[b]; nval[b+1] = nval[b]; --b;
        }
        ngid[b+1] = g; nval[b+1] = w;
    }
    real_t val = 0.0;
    for (int a = 0; a < cnt; ++a) val += nval[a];
    *out = val / (real_t)cnt;
    return cnt;
}

static void extrap_det_exchange_layer(const struct fesom_mesh *mesh,
                                      struct fesom_partit *partit, real_t *layer)
{
    if (!partit || partit->npes <= 1) return;
    fesom_halo_exchange(layer, FESOM_HALO_NOD2D, 1, 1, partit);
}

static void extrap_nod3D_det(const struct fesom_mesh *mesh,
                             struct fesom_partit     *partit,
                             real_t                  *arr)
{
    const int N       = mesh->myDim_nod2D;
    const int N_alloc = mesh->myDim_nod2D + mesh->eDim_nod2D;
    const int nl      = mesh->nl;

    /* Relaxation tolerance: after the ring fill, filled nodes are relaxed
     * (Jacobi mean over eligible neighbours; ORIGINAL valid nodes stay fixed
     * as boundary data) until the global max per-sweep change < tol — this
     * takes the fill to the smooth quasi-harmonic interior interpolation.
     * Without it, ring filling from two water masses meets in an equidistant
     * front of maximal contrast across ~one ring (measured on NG5: step-1
     * |uv| 3.06 m/s and blow-up at step 9, partition-INdependent). */
    static real_t s_tol = -1.0;
    if (s_tol < 0.0) {
        const char *e = getenv("FESOM_IC_EXTRAP_TOL");
        s_tol = (e && e[0]) ? (real_t)atof(e) : (real_t)1e-3;
        if (s_tol <= 0.0) s_tol = (real_t)1e-3;
    }

    real_t *cur  = (real_t *)malloc((size_t)N_alloc * sizeof(real_t));
    real_t *nxt  = (real_t *)malloc((size_t)N_alloc * sizeof(real_t));
    char   *fill = (char *)malloc((size_t)N_alloc);
    FESOM_CHECK(cur && nxt && fill, "extrap_det: oom (layer buffers)");

    long tot_fill_sweeps = 0, tot_relax_sweeps = 0;
    for (int nz = 0; nz < nl - 1; ++nz) {
        for (int n = 0; n < N_alloc; ++n) cur[n] = arr[FESOM_NODE3D(n, nz, nl)];
        extrap_det_exchange_layer(mesh, partit, cur);
        for (int n = 0; n < N_alloc; ++n)
            fill[n] = (cur[n] > 0.99 * PHC_DUMMY) ? 1 : 0;   /* dummy at entry */

        /* Phase 1 — ring fill (Jacobi, validity frozen per sweep). */
        for (;;) {
            ++tot_fill_sweeps;
            FESOM_CHECK(tot_fill_sweeps < 1000000, "extrap_det: fill sweep cap");
            memcpy(nxt, cur, (size_t)N_alloc * sizeof(real_t));
            int progress = 0;
            for (int n = 0; n < N; ++n) {
                if (!(cur[n] > 0.99 * PHC_DUMMY) ||
                    !(mesh->nlevels_nod2D[n] - 1 > nz)) continue;
                real_t v;
                if (extrap_det_ring_mean(mesh, partit, cur, n, nz, N_alloc, &v)) {
                    nxt[n] = v; progress = 1;
                }
            }
            if (partit && partit->npes > 1)
                MPI_Allreduce(MPI_IN_PLACE, &progress, 1, MPI_INT, MPI_MAX,
                              partit->MPI_COMM_FESOM);
            if (!progress) break;
            memcpy(cur, nxt, (size_t)N_alloc * sizeof(real_t));
            extrap_det_exchange_layer(mesh, partit, cur);
        }

        /* Phase 2 — relaxation of the FILLED nodes toward the harmonic
         * interpolation of the fixed original data. Still-dummy nodes
         * (no source at this layer) are left for the vertical fill. */
        for (int n = 0; n < N; ++n)
            if (cur[n] > 0.99 * PHC_DUMMY) fill[n] = 0;      /* unfillable */
        long relax_cap = tot_relax_sweeps + 20000;
        for (;;) {
            real_t dmax = 0.0;
            memcpy(nxt, cur, (size_t)N_alloc * sizeof(real_t));
            for (int n = 0; n < N; ++n) {
                if (!fill[n] || !(mesh->nlevels_nod2D[n] - 1 > nz)) continue;
                real_t v;
                if (extrap_det_ring_mean(mesh, partit, cur, n, nz, N_alloc, &v)) {
                    nxt[n] = v;
                    real_t d = fabs(v - cur[n]);
                    if (d > dmax) dmax = d;
                }
            }
            if (partit && partit->npes > 1)
                MPI_Allreduce(MPI_IN_PLACE, &dmax, 1, MPI_DOUBLE, MPI_MAX,
                              partit->MPI_COMM_FESOM);
            if (dmax <= s_tol) break;
            ++tot_relax_sweeps;
            if (tot_relax_sweeps >= relax_cap) {
                if (!partit || partit->mype == 0)
                    printf("[fesom_phc] det extrap: relax cap at nz=%d (dmax=%.3e) — proceeding\n",
                           nz, (double)dmax);
                break;
            }
            memcpy(cur, nxt, (size_t)N_alloc * sizeof(real_t));
            extrap_det_exchange_layer(mesh, partit, cur);
        }

        for (int n = 0; n < N; ++n) arr[FESOM_NODE3D(n, nz, nl)] = cur[n];
    }
    free(cur); free(nxt); free(fill);
    if (!partit || partit->mype == 0) {
        printf("[fesom_phc] det extrap: %ld fill sweeps + %ld relax sweeps (tol %.1e)\n",
               tot_fill_sweeps, tot_relax_sweeps, (double)s_tol);
    }

    /* Vertical fill — identical to legacy (order-independent). */
    for (int n = 0; n < N; ++n) {
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = 1; nz < nl1; ++nz) {
            if (arr[FESOM_NODE3D(n, nz, nl)] > 0.99 * PHC_DUMMY) {
                arr[FESOM_NODE3D(n, nz, nl)] = arr[FESOM_NODE3D(n, nz - 1, nl)];
            }
        }
    }
    extrap_det_exchange_all(mesh, partit, arr);   /* final exchange */
}

/* FESOM_IC_EXTRAP knob: "det" = deterministic fill; unset/"legacy" = the
 * literal-port Gauss-Seidel. Any other value refuses loudly (house rule:
 * a knob announces itself or aborts — never dies silent). */
static int ic_extrap_det_on(struct fesom_partit *partit)
{
    static int s_det = -1;
    if (s_det < 0) {
        const char *e = getenv("FESOM_IC_EXTRAP");
        if (!e || !e[0] || strcmp(e, "legacy") == 0) {
            s_det = 0;
        } else if (strcmp(e, "det") == 0) {
            s_det = 1;
            if (!partit || partit->mype == 0) {
                printf("[fesom_phc] FESOM_IC_EXTRAP=det — deterministic "
                       "partition-independent hole fill (ring-Jacobi, gid-ordered sums)\n");
            }
        } else {
            FESOM_DIE("FESOM_IC_EXTRAP=%s not recognised (want: det | legacy)", e);
        }
    }
    return s_det;
}

/*--- Top-level entry -----------------------------------------------------*/
void fesom_phc_load_ic(const char                  *path,
                       const struct fesom_mesh     *mesh,
                       struct fesom_tracers        *tracers,
                       int                          t_insitu,
                       struct fesom_partit         *partit)
{
    int ncid;
    NC_CHECK(nc_open(path, NC_NOWRITE, &ncid));

    /* nc_readGrid (gen_ic3d.F90:68-217). */
    int Nlon = 0, Nlat = 0, Ndepth = 0;
    {
        int dimid;
        size_t len;
        /* Try LAT, lat, latitude (Fortran 95-103). */
        if (nc_inq_dimid(ncid, "LAT", &dimid) != NC_NOERR
         && nc_inq_dimid(ncid, "lat", &dimid) != NC_NOERR
         && nc_inq_dimid(ncid, "latitude", &dimid) != NC_NOERR)
            FESOM_DIE("phc: no lat dimension in %s", path);
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len)); Nlat = (int)len;

        if (nc_inq_dimid(ncid, "LON", &dimid) != NC_NOERR
         && nc_inq_dimid(ncid, "lon", &dimid) != NC_NOERR
         && nc_inq_dimid(ncid, "longitude", &dimid) != NC_NOERR)
            FESOM_DIE("phc: no lon dimension in %s", path);
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len)); Nlon = (int)len;

        if (nc_inq_dimid(ncid, "DEPTH", &dimid) != NC_NOERR
         && nc_inq_dimid(ncid, "depth", &dimid) != NC_NOERR)
            FESOM_DIE("phc: no depth dimension in %s", path);
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &len)); Ndepth = (int)len;
    }
    printf("[fesom_phc] %s : Nlon=%d  Nlat=%d  Ndepth=%d\n", path, Nlon, Nlat, Ndepth);

    /* Cyclic halo padding (Fortran gen_ic3d.F90:167,186-189,213-215): the lon
       axis is periodic, so allocate Nlon+2 and keep ALL `Nlon` real columns in
       nc_lon[1..Nlon], with wrap-around halo cells at [0] (= last real − 360)
       and [Nlon+1] (= first real + 360). ncdata is padded identically in
       load_one_variable. (The previous code dropped the first/last real columns
       and shifted them by ±360 in place — breaking the interpolation within ~2°
       of lon 0/360, e.g. the western Med.) */
    double *nc_lon   = malloc((size_t)(Nlon + 2) * sizeof(double));
    double *nc_lat   = malloc((size_t)Nlat   * sizeof(double));
    double *nc_depth = malloc((size_t)Ndepth * sizeof(double));
    FESOM_CHECK(nc_lon && nc_lat && nc_depth, "phc: oom (coord arrays)");
    {
        int varid;
        if (nc_inq_varid(ncid, "lon", &varid) != NC_NOERR) NC_CHECK(nc_inq_varid(ncid, "LON", &varid));
        NC_CHECK(nc_get_var_double(ncid, varid, nc_lon + 1));   /* real lons -> [1..Nlon] */
        if (nc_inq_varid(ncid, "lat", &varid) != NC_NOERR) NC_CHECK(nc_inq_varid(ncid, "LAT", &varid));
        NC_CHECK(nc_get_var_double(ncid, varid, nc_lat));
        if (nc_inq_varid(ncid, "depth", &varid) != NC_NOERR) NC_CHECK(nc_inq_varid(ncid, "DEPTH", &varid));
        NC_CHECK(nc_get_var_double(ncid, varid, nc_depth));
    }
    nc_lon[0]        = nc_lon[Nlon] - 360.0;   /* west wrap halo */
    nc_lon[Nlon + 1] = nc_lon[1]    + 360.0;   /* east wrap halo */
    Nlon += 2;   /* Nlon is now the PADDED count; real cells at [1..Nlon-2] */

    /* nc_ic3d_ini: per-node bilin_indx (Fortran 281-299, non-cavity branch). */
    int *bilin_i = malloc((size_t)mesh->myDim_nod2D * sizeof(int));
    int *bilin_j = malloc((size_t)mesh->myDim_nod2D * sizeof(int));
    FESOM_CHECK(bilin_i && bilin_j, "phc: oom (bilin)");

    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        double x = (double)mesh->geo_coord_nod2D[2*n + 0] / FESOM_RAD;
        double y = (double)mesh->geo_coord_nod2D[2*n + 1] / FESOM_RAD;
        if (x < 0.0)   x += 360.0;
        if (x > 360.0) x -= 360.0;
        bilin_i[n] = (x <= nc_lon[Nlon - 1] && x >= nc_lon[0])
                       ? binarysearch_d(Nlon, nc_lon, x) : -1;
        bilin_j[n] = (y <= nc_lat[Nlat - 1] && y >= nc_lat[0])
                       ? binarysearch_d(Nlat, nc_lat, y) : -1;
        /* Reject the rightmost / topmost bracket (would index past array). */
        if (bilin_i[n] >= Nlon - 1) bilin_i[n] = -1;
        if (bilin_j[n] >= Nlat - 1) bilin_j[n] = -1;
    }

    /* getcoeffld for T and S. PHC variable names: "temp" and "salt". */
    real_t *T = tracers->data[FESOM_TRACER_T].values;
    real_t *S = tracers->data[FESOM_TRACER_S].values;

    /* Look up _FillValue if present. */
    int has_fill_T = 0, has_fill_S = 0;
    double fill_T = 0.0, fill_S = 0.0;
    {
        int vid;
        if (nc_inq_varid(ncid, "temp", &vid) == NC_NOERR) {
            if (nc_get_att_double(ncid, vid, "_FillValue", &fill_T) == NC_NOERR) has_fill_T = 1;
        }
        if (nc_inq_varid(ncid, "salt", &vid) == NC_NOERR) {
            if (nc_get_att_double(ncid, vid, "_FillValue", &fill_S) == NC_NOERR) has_fill_S = 1;
        }
    }

    load_one_variable(ncid, "temp", mesh, nc_lon, nc_lat, nc_depth,
                      Nlon, Nlat, Ndepth, bilin_i, bilin_j,
                      fill_T, has_fill_T, T);
    load_one_variable(ncid, "salt", mesh, nc_lon, nc_lat, nc_depth,
                      Nlon, Nlat, Ndepth, bilin_i, bilin_j,
                      fill_S, has_fill_S, S);

    NC_CHECK(nc_close(ncid));

    /* Debug dump: post-bilin pre-extrap surface T,S. */
    {
        const char *dd = getenv("FESOM_EVP_DUMP_DIR");
        if (dd && partit) {
            char dpath[1024];
            snprintf(dpath, sizeof(dpath),
                     "%s/phc_dump_preextrap_rank%d.txt", dd, partit->mype);
            FILE *df = fopen(dpath, "w");
            if (df) {
                fprintf(df, "# var=T_surface,S_surface,bilin_i,bilin_j,lon,lat rank=%d N=%d\n",
                        partit->mype, mesh->myDim_nod2D);
                for (int n = 0; n < mesh->myDim_nod2D; ++n) {
                    int gid = partit->myList_nod2D[n];
                    size_t k = FESOM_NODE3D(n, 0, mesh->nl);
                    double xlon = (double)mesh->geo_coord_nod2D[2*n + 0]
                                / FESOM_RAD;
                    double xlat = (double)mesh->geo_coord_nod2D[2*n + 1]
                                / FESOM_RAD;
                    fprintf(df, "%d %.17g %.17g %d %d %.17g %.17g\n", gid,
                            (double)T[k], (double)S[k],
                            bilin_i[n], bilin_j[n], xlon, xlat);
                }
                fclose(df);
            }
        }
    }

    /* Extrapolate to fill any remaining dummy-tagged ocean nodes. */
    if (ic_extrap_det_on(partit)) {
        extrap_nod3D_det(mesh, partit, T);
        extrap_nod3D_det(mesh, partit, S);
    } else {
        extrap_nod3D(mesh, partit, T);
        extrap_nod3D(mesh, partit, S);
    }

    /* Final cleanup (Fortran do_ic3d lines 536-557): replace remaining dummy
     * with 0, zero-out below nlevels_nod2D, K → C if T > 100.
     * Fortran iterates myDim+eDim on each sweep (the `where (...)` applies
     * to the full array). If we only clean myDim, halo entries keep raw PHC
     * values (in-situ, possibly Kelvin, possibly dummy) and downstream code
     * that reads T/S at halo positions sees garbage — partition-dependent. */
    const int N_full_ic = mesh->myDim_nod2D + mesh->eDim_nod2D;
    for (int n = 0; n < N_full_ic; ++n) {
        int nl1 = mesh->nlevels_nod2D[n] - 1;
        for (int nz = 0; nz < mesh->nl; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, mesh->nl);
            if (T[k] > 0.9 * PHC_DUMMY) T[k] = 0.0;
            if (S[k] > 0.9 * PHC_DUMMY) S[k] = 0.0;
            if (nz >= nl1) {
                T[k] = 0.0;
                S[k] = 0.0;
            }
            if (T[k] > 100.0) T[k] -= 273.15;
        }
    }

    /* insitu→pot via ptheta (Fortran insitu2pot, lines 3074-3118).
     * Fortran insitu2pot loops myDim+eDim (oce_ale_pressure_bv.F90:3102). */
    if (t_insitu) {
        printf("[fesom_phc] converting in-situ → potential temperature\n");
        for (int n = 0; n < N_full_ic; ++n) {
            int nzmin = mesh->ulevels_nod2D[n] - 1;
            int nzmax = mesh->nlevels_nod2D[n] - 1;
            for (int nz = nzmin; nz < nzmax; ++nz) {
                size_t k = FESOM_NODE3D(n, nz, mesh->nl);
                double tt = (double)T[k];
                double ss = (double)S[k];
                double pp = fabs((double)mesh->Z[nz]);
                T[k] = (real_t)ptheta(ss, tt, pp, 0.0);
            }
        }
    }

    free(bilin_i);
    free(bilin_j);
    free(nc_lon);
    free(nc_lat);
    free(nc_depth);

    /* Stats summary (Fortran lines 586-590). */
    real_t Tmin = T[0], Tmax = T[0], Smin = S[0], Smax = S[0];
    for (int n = 0; n < mesh->myDim_nod2D; ++n) {
        int nzmax = mesh->nlevels_nod2D[n] - 1;
        for (int nz = 0; nz < nzmax; ++nz) {
            size_t k = FESOM_NODE3D(n, nz, mesh->nl);
            if (T[k] < Tmin) Tmin = T[k];
            if (T[k] > Tmax) Tmax = T[k];
            if (S[k] < Smin) Smin = S[k];
            if (S[k] > Smax) Smax = S[k];
        }
    }
    printf("[fesom_phc] T ∈ [%.4f, %.4f] °C   S ∈ [%.4f, %.4f] PSU\n",
           (double)Tmin, (double)Tmax, (double)Smin, (double)Smax);

    /* Debug dump: post-PHC-load surface T per gid (gated on FESOM_EVP_DUMP_DIR). */
    {
        const char *dd = getenv("FESOM_EVP_DUMP_DIR");
        if (dd && partit) {
            char dpath[1024];
            snprintf(dpath, sizeof(dpath),
                     "%s/phc_dump_postload_rank%d.txt", dd, partit->mype);
            FILE *df = fopen(dpath, "w");
            if (df) {
                fprintf(df, "# var=T_surface,S_surface rank=%d N=%d\n",
                        partit->mype, mesh->myDim_nod2D);
                for (int n = 0; n < mesh->myDim_nod2D; ++n) {
                    int gid = partit->myList_nod2D[n];
                    size_t k = FESOM_NODE3D(n, 0, mesh->nl);
                    fprintf(df, "%d %.17g %.17g\n", gid,
                            (double)T[k], (double)S[k]);
                }
                fclose(df);
            }
        }
    }
}
