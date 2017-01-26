/**
 * @mainpage
 * TDPIC: 三次元Electromagnetic Full Particle-in-cell code
 * 帯電解析用 & マルチグリッド解析用のPICコード
 **/

#include <tdpic.h>

#ifndef BUILD_TEST
int main(int argc, char* argv[]){
    //! MPI Environmentを初期化
    MPIw::Environment mpiEnv(argc, argv);
    MPIw::Communicator world;

    Grid* root_grid;
    Initializer::initTDPIC(root_grid);

    if( Environment::isRootNode ) {
        cout << "--  Begin A Loop  --" << endl;
    }

    // particle -> space charge
    root_grid->updateRho();

    // space charge -> potential
    root_grid->solvePoisson();

    // potential -> efield
    root_grid->updateEfield();

    if( Environment::isRootNode ) {
        cout << "--  End A Loop  --" << endl;

#ifdef DEBUG
        IO::outputParticlePositions(root_grid->particles);
#endif

        // Level 2まで
        root_grid->makeChild(2, 2, 2, 8, 8, 8);
        // root_grid->makeChild(4, 8, 8, 9, 9, 9);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 10, 10, 10);
        cout << root_grid << endl;

        root_grid->copyScalarToChildren("potential");
        root_grid->getChildren()[0]->copyScalarToChildren("potential");

        // root_grid->getChildren()[0]->makeChild(8, 8, 8, 7, 7, 7);
    }

    IO::writeDataInParallel(root_grid, 0, "potential");
    IO::writeDataInParallel(root_grid, 0, "rho");
    IO::writeDataInParallel(root_grid, 0, "efield");
    IO::writeDataInParallel(root_grid, 0, "bfield");
#ifdef DEBUG
    world.barrier();
#endif
    return 0;
}
#endif
