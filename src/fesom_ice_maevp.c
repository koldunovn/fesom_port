#include "fesom_ice_maevp.h"
#include "fesom_mesh.h"
#include "fesom_partit.h"

#include <mpi.h>
#include <stdio.h>

/*
 * M1 skeleton: dispatch target for whichEVP == 1. The literal port of
 * EVPdynamics_m (ice_maEVP.F90:429-882) lands in Phase M2; until then any
 * run that reaches this point aborts loudly rather than silently running
 * the wrong rheology.
 */
void fesom_ice_evp_dynamics_m(fesom_ice            *ice,
                              struct fesom_partit  *partit,
                              struct fesom_mesh    *mesh)
{
    (void)ice; (void)mesh;
    if (partit->mype == 0)
        fprintf(stderr, "fesom_ice_maevp: mEVP port pending (M2) — "
                        "whichEVP=1 selected but EVPdynamics_m is not ported yet\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
}
