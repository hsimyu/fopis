#include <tdpic.h>
#include <gtest/gtest.h>

using std::cout;
using std::endl;

namespace {

    class FieldTest : public ::testing::Test {
        protected:
            FieldTest() {
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

                env->ptype = ptype;

                root_grid = Initializer::initializeGrid(env);
                Initializer::initializeRootField(env, root_grid);

                cout << env << endl;
            }

            // members
            Environment* env;
            ParticleType* ptype;
            Grid* root_grid;
    };

    TEST_F(FieldTest, FDTD)
    {
        cout << env << endl;

        cout << "--  Begin A Loop  --" << endl;
        //
        // // particle -> space charge
        // root_grid->updateRho(env);
        //
        // // space charge -> potential
        // root_grid->getField()->solvePoisson(env);
        //
        // // potential -> efield
        // root_grid->getField()->updateEfield(env);
        //
        // cout << "--  End A Loop  --" << endl;

        ASSERT_EQ(1, 1);
        IO::outputParticlePositions( env, root_grid->particles );

#ifdef DEBUG
        // cout << "-- SPACE CHARGE --" << endl;
        // IO::print3DArray( root_grid->getField()->getRho(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
        // cout << "-- POTENTIAL --" << endl;
        // IO::print3DArray( root_grid->getField()->getPhi(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
        cout << "-- Ex --" << endl;
        IO::print3DArray( root_grid->getField()->getEx(), env->cell_x + 1, env->cell_y + 2, env->cell_z + 2);
        cout << "-- Ey --" << endl;
        IO::print3DArray( root_grid->getField()->getEy(), env->cell_x + 2, env->cell_y + 1, env->cell_z + 2);
        cout << "-- Ez --" << endl;
        IO::print3DArray( root_grid->getField()->getEz(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 1);
#endif
    }
}
