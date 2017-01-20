#include <test_tdpic.h>

using std::cout;
using std::endl;

namespace TEST_TDPIC {
    TDPICTest1::TDPICTest1() {
        env = new Environment;

        // test input
        env->jobtype = "new";
        env->nx = 16;
        env->ny = 16;
        env->nz = 16;
        env->proc_x = 1;
        env->proc_y = 1;
        env->proc_z = 1;
        env->solver_type = "EM";
        env->boundary = "DDDDDD";
        env->dimension = "3D";
        env->dx = 0.1;
        env->dt = 1e-8;
        env->cell_x = env->nx/env->proc_x;
        env->cell_y = env->ny/env->proc_y;
        env->cell_z = env->nz/env->proc_z;

        env->num_of_particle_types = 2;
        ptype = new ParticleType[env->num_of_particle_types];

        ptype[0].setId(0);
        ptype[0].setName("Electron");
        ptype[0].setType("ambient");
        ptype[0].setMass(1.0);
        ptype[0].setCharge(-1.0);
        ptype[0].setTemperature(1.0);
        ptype[0].setDensity(1.0e6);
        ptype[0].setPcell(20);
        ptype[0].calcTotalNumber(env);
        ptype[0].calcSize(env);

        ptype[1].setId(1);
        ptype[1].setName("Proton");
        ptype[1].setType("ambient");
        ptype[1].setMass(1836.0 * 1.0);
        ptype[1].setCharge(1.0);
        ptype[1].setTemperature(1.0);
        ptype[1].setDensity(1.0e6);
        ptype[1].setPcell(20);
        ptype[1].calcTotalNumber(env);
        ptype[1].calcSize(env);

        // initialize normalizer
        Utils::Normalizer::x_unit = env->dx;
        Utils::Normalizer::t_unit = env->dt;
        Utils::Normalizer::e_unit = e;

        env->ptype = ptype;
        root_grid = Initializer::initializeGrid(env);
    }

    // テスト毎に実行される
    TDPICTest1::~TDPICTest1() {
        Grid::resetNextID();
    }
}
