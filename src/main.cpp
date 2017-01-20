/**
 * @mainpage
 * TDPIC: 三次元Electromagnetic Full Particle-in-cell code
 * 帯電解析用 & マルチグリッド解析用のPICコード
 **/

#include <tdpic.h>

#ifndef BUILD_TEST
int main(int argc, char* argv[]){
    //! MPI Environmentを初期化
    //! コンストラクタでMPI_Initを呼んで
    //! mainが終わったらdestructされる

    MPI::Environment mpiEnv(argc, argv);
    MPI::Communicator world; // MPI_COMM_WORLD

    Environment* env;
    ParticleType* ptype;
    Grid* root_grid;

    Initializer::initTDPIC(env, ptype, root_grid);

    if( env->isRootNode ) {
        cout << "--  Begin A Loop  --" << endl;
    }

    // particle -> space charge
    root_grid->updateRho(env);

    // space charge -> potential
    root_grid->getField()->solvePoisson(env);

    // potential -> efield
    root_grid->getField()->updateEfield(env);

    if( env->isRootNode ) {
        cout << "--  End A Loop  --" << endl;

#ifdef DEBUG
        // cout << "-- SPACE CHARGE --" << endl;
        // IO::print3DArray( root_grid->getField()->getRho(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
        // cout << "-- POTENTIAL --" << endl;
        // IO::print3DArray( root_grid->getField()->getPhi(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
        // cout << "-- Ex --" << endl;
        // IO::print3DArray( root_grid->getField()->getEx(), env->cell_x + 1, env->cell_y + 2, env->cell_z + 2);
        // cout << "-- Ey --" << endl;
        // IO::print3DArray( root_grid->getField()->getEy(), env->cell_x + 2, env->cell_y + 1, env->cell_z + 2);
        // cout << "-- Ez --" << endl;
        // IO::print3DArray( root_grid->getField()->getEz(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 1);
        // cout << "-- Bx --" << endl;
        // IO::print3DArray( root_grid->getField()->getBx(), env->cell_x + 2, env->cell_y + 1, env->cell_z + 1);
        // cout << "-- By --" << endl;
        // IO::print3DArray( root_grid->getField()->getBy(), env->cell_x + 1, env->cell_y + 2, env->cell_z + 1);
        // cout << "-- Bz --" << endl;
        // IO::print3DArray( root_grid->getField()->getBz(), env->cell_x + 1, env->cell_y + 1, env->cell_z + 2);
        IO::outputParticlePositions( env, root_grid->particles );
#endif

        // Level 2まで
        root_grid->makeChild(2, 2, 2, 5, 5, 5);
        root_grid->makeChild(4, 8, 8, 9, 9, 9);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);
        cout << root_grid << endl;
    }

    // IO::writeDataInParallel(root_grid, 0, "potential");
    return 0;
}
#endif
