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
    auto time_counter = Utils::TimeCounter::getInstance();

    if( Environment::isRootNode ) {
        cout << "===== Begin Main Loop =====" << endl;
        time_counter->enableReport();
    }

    root_grid->makeChild(2, 2, 2, 8, 8, 15);
    // root_grid->printInfo();

    // initialized rho, phi, efield
    time_counter->begin("updateRho");
        root_grid->updateRho();
    time_counter->switchTo("solvePoisson");
        root_grid->solvePoisson();
    time_counter->switchTo("updateEfield");
        root_grid->updateEfield();
    time_counter->switchTo("updateBfield");
        root_grid->updateBfield();
    time_counter->end();

    for(; Environment::timestep <= Environment::max_timestep; ++Environment::timestep) {
        if( Environment::isRootNode ) {
            cout << "====== Iteration " << Environment::timestep << " =====" << endl;
        }
        root_grid->mainLoop();

        time_counter->begin("plotData");
        if(Environment::plotPotential())    IO::writeDataInParallel(root_grid, Environment::timestep, "potential");
        if(Environment::plotRho())          IO::writeDataInParallel(root_grid, Environment::timestep, "rho");
        if(Environment::plotEfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "efield");
        if(Environment::plotDensity())      IO::writeDataInParallel(root_grid, Environment::timestep, "density");
        if(Environment::plotEnergy())       IO::plotEnergy(root_grid, Environment::timestep);
        if(Environment::plotEnergyDist())   IO::plotParticleEnergyDistribution(root_grid->getParticles());
        if(Environment::plotVelocityDist()) IO::plotParticleVelocityDistribution(root_grid->getParticles());

        // 磁場プロット
        if(Environment::solver_type == "EM" && Environment::plotBfield()) {
            IO::writeDataInParallel(root_grid, Environment::timestep, "bfield");
        }

        IO::plotObjectsData(root_grid);
        IO::plotValidParticleNumber(root_grid);
        time_counter->end();
    }

    if( Environment::isRootNode ) {
        cout << "===== End Main Loop =====" << endl;
    }

    time_counter->destroy();

    return 0;
}
#endif
