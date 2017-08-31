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

    std::shared_ptr<RootGrid> root_grid = Initializer::initTDPIC();

    if( Environment::isRootNode ) {
        cout << "--  Begin Main Loop  --" << endl;
    }

    // initialized rho, phi, efield
    root_grid->updateRho();
    root_grid->solvePoisson();
    root_grid->updateEfield();
    root_grid->updateBfield();

    for(; Environment::timestep <= Environment::max_timestep; ++Environment::timestep) {
        if( Environment::isRootNode ) {
            cout << "--  Iteration " << Environment::timestep << "  --" << endl;
        }
        root_grid->resetObjects();

        // timing: t + 0.5 dt
        root_grid->updateParticleVelocity(); // 速度更新
        root_grid->updateParticlePosition(); // jx, jy, jz もここで update される
        root_grid->emitParticlesFromObjects();

        if ( Environment::solver_type == "EM" ) {
            // 電磁計算の場合
            // timing: t + dt

            if ( Environment::timestep % 1 == 0 ) {
                root_grid->updateRho(); // 新しい位置に対応する電荷密度算出
                root_grid->solvePoisson(); // Poisson を解く (FDTDの場合はたまにでいい?)
                root_grid->updateEfield(); // 電場更新
            } else {
                root_grid->updateEfieldFDTD(); // FDTDで電場更新
            }
            root_grid->updateBfield(); // 磁場更新
            if(Environment::plotBfield()) IO::writeDataInParallel(root_grid, Environment::timestep, "bfield");
        } else {
            // 静電計算の場合
            // timing: t + dt
            root_grid->updateRho(); // 新しい位置に対応する電荷密度算出
            root_grid->solvePoisson(); // Poisson を解く
            root_grid->updateEfield(); // 電場更新
        }

        if(Environment::plotPotential())    IO::writeDataInParallel(root_grid, Environment::timestep, "potential");
        if(Environment::plotRho())          IO::writeDataInParallel(root_grid, Environment::timestep, "rho");
        if(Environment::plotEfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "efield");
        if(Environment::plotDensity())      IO::writeDataInParallel(root_grid, Environment::timestep, "density");
        if(Environment::plotEnergy())       IO::plotEnergy(root_grid, Environment::timestep);
        if(Environment::plotEnergyDist())   IO::plotParticleEnergyDistribution(root_grid->getParticles());
        if(Environment::plotVelocityDist()) IO::plotParticleVelocityDistribution(root_grid->getParticles());

        IO::plotObjectsData(root_grid);
        IO::plotValidParticleNumber(root_grid);
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
