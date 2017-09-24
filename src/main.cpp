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
#include "utils.hpp"

#ifndef BUILD_TEST
int main(int argc, char* argv[]){
    //! MPI Environmentを初期化
    MPIw::Environment mpiEnv(argc, argv);

    std::shared_ptr<RootGrid> root_grid = Initializer::initTDPIC();
    Utils::TimeCounter time_counter;

    if( Environment::isRootNode ) {
        cout << "===== Begin Main Loop =====" << endl;
        time_counter.enableReport();
    }

    // root_grid->makeChild(2, 2, 2, 8, 8, 15);
    // root_grid->makeChild(15, 15, 8, 18, 20, 25);
    // root_grid->printInfo();

    // initialized rho, phi, efield
    time_counter.begin("updateRho");
        root_grid->updateRho();
    time_counter.switchTo("solvePoisson");
        root_grid->solvePoisson();
    time_counter.switchTo("updateEfield");
        root_grid->updateEfield();
    time_counter.switchTo("updateBfield");
        root_grid->updateBfield();
    time_counter.end();

    for(; Environment::timestep <= Environment::max_timestep; ++Environment::timestep) {
        if( Environment::isRootNode ) {
            cout << "====== Iteration " << Environment::timestep << " =====" << endl;
        }
        time_counter.begin("resetObjects");
        root_grid->resetObjects();

        // -- timing: t + 0.5 dt --
        // 速度更新
        time_counter.switchTo("updateParticleVelocity");
        root_grid->updateParticleVelocity();

        // 位置更新
        time_counter.switchTo("updateParticlePosition");
        root_grid->updateParticlePosition(); // jx, jy, jz もここで update される

        // 粒子注入
        time_counter.switchTo("injectParticleFromBoundary");
        root_grid->injectParticlesFromBoundary();

        // 粒子放出
        time_counter.switchTo("emitParticlesFromObjects");
        root_grid->emitParticlesFromObjects();

        // -- timing: t + dt --
        if ( Environment::solver_type == "EM" ) {
            // 電磁計算の場合はFDTDとPoisson解くのを分ける
            if ( Environment::timestep % 1 == 0 ) {
                // 新しい位置に対応する電荷密度算出
                time_counter.switchTo("updateRho");
                root_grid->updateRho();

                // Poisson を解く (FDTDの場合はたまにでいい?)
                time_counter.switchTo("solvePoisson");
                root_grid->solvePoisson();

                // 電場更新
                time_counter.switchTo("updateEfield");
                root_grid->updateEfield();
            } else {
                // FDTDで電場更新
                time_counter.switchTo("updateEfieldFDTD");
                root_grid->updateEfieldFDTD();
            }

            // 磁場更新
            time_counter.switchTo("updateBfield");
            root_grid->updateBfield();

            // 磁場プロット
            time_counter.switchTo("plotData");
            if(Environment::plotBfield()) {
                IO::writeDataInParallel(root_grid, Environment::timestep, "bfield");
            }
        } else {
            // 静電計算の場合
            // 新しい位置に対応する電荷密度算出
            time_counter.switchTo("updateRho");
            root_grid->updateRho();

            // Poisson を解く
            time_counter.switchTo("solvePoisson");
            root_grid->solvePoisson();

            // 電場更新
            time_counter.switchTo("updateEfield");
            root_grid->updateEfield();
        }

        time_counter.switchTo("plotData");
        if(Environment::plotPotential())    IO::writeDataInParallel(root_grid, Environment::timestep, "potential");
        if(Environment::plotRho())          IO::writeDataInParallel(root_grid, Environment::timestep, "rho");
        if(Environment::plotEfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "efield");
        if(Environment::plotDensity())      IO::writeDataInParallel(root_grid, Environment::timestep, "density");
        if(Environment::plotEnergy())       IO::plotEnergy(root_grid, Environment::timestep);
        if(Environment::plotEnergyDist())   IO::plotParticleEnergyDistribution(root_grid->getParticles());
        if(Environment::plotVelocityDist()) IO::plotParticleVelocityDistribution(root_grid->getParticles());

        IO::plotObjectsData(root_grid);
        IO::plotValidParticleNumber(root_grid);
        time_counter.end();
    }

    if( Environment::isRootNode ) {
        cout << "===== End Main Loop =====" << endl;
    }

    return 0;
}
#endif
