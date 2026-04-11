#include <string.h>

#include "prj.h"

int main(void)
{
    prj_sim sim;

    memset(&sim, 0, sizeof(sim));
    sim.coord.x1min = -1.0;
    sim.coord.x1max = 1.0;
    sim.coord.x2min = -1.0;
    sim.coord.x2max = 1.0;
    sim.coord.x3min = -1.0;
    sim.coord.x3max = 1.0;
    if (prj_mesh_init(&sim.mesh, 1, 1, 1, 0, &sim.coord) != 0) {
        return 1;
    }
    prj_gravity_init(&sim);
    if (sim.monopole_grav.nbins <= 0 || sim.monopole_grav.rf == 0 || sim.monopole_grav.accel == 0) {
        prj_mesh_destroy(&sim.mesh);
        return 2;
    }
    prj_gravity_free(&sim.monopole_grav);
    prj_mesh_destroy(&sim.mesh);
    return prj_gravity_apply();
}
