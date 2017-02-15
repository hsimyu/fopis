/**
 * @mainpage
 * TDPIC: 三次元Electromagnetic Full Particle-in-cell code
 * 帯電解析用 & マルチグリッド解析用のPICコード
 **/

#include "global.hpp"
#include "initialize.hpp"
#include "mpiw.hpp"
#include "grid.hpp"
#include "environment.hpp"
#include "dataio.hpp"

#ifndef BUILD_TEST
int main(int argc, char* argv[]){
    //! MPI Environmentを初期化
    MPIw::Environment mpiEnv(argc, argv);

    Grid* root_grid;
    Initializer::initTDPIC(root_grid);

    if( Environment::isRootNode ) {
        cout << "--  Begin A Loop  --" << endl;
    }

    // first update
    root_grid->updateRho();
    root_grid->solvePoisson();
    root_grid->updateEfield();

    for(; Environment::timestep < Environment::max_iteration; ++Environment::timestep) {
        if( Environment::isRootNode ) {
            cout << "--  Iteration " << Environment::timestep << "  --" << endl;
        }
        // new particle position
        root_grid->updateParticleVelocity();
        root_grid->updateParticlePosition();
        root_grid->updateRho();
        root_grid->solvePoisson();
        root_grid->updateEfield();
        root_grid->updateBfield();

        if(Environment::plotPotential())    IO::writeDataInParallel(root_grid, Environment::timestep, "potential");
        if(Environment::plotRho())          IO::writeDataInParallel(root_grid, Environment::timestep, "rho");
        if(Environment::plotEfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "efield");
        if(Environment::plotBfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "bfield");
        if(Environment::plotEnergy())       IO::plotEnergy(*root_grid, Environment::timestep);
        if(Environment::plotEnergyDist())   IO::plotParticleEnergyDistribution(root_grid->particles);
        if(Environment::plotVelocityDist()) IO::plotParticleVelocityDistribution(root_grid->particles);
    }

    /*
    if( !Environment::isRootNode ) {
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
    }*/
    return 0;
}
#endif
