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
    std::string input_filename = "input.json";
    if (argc > 1) {
        input_filename = argv[1];
    }

    //! MPI Environmentを初期化
    MPIw::Environment mpiEnv(argc, argv);

    std::shared_ptr<RootGrid> root_grid = Initializer::initTDPIC(input_filename);
    auto time_counter = Utils::TimeCounter::getInstance();

    if( Environment::isRootNode ) {
        cout << "===== Begin Main Loop =====" << endl;
        time_counter->enableReport();
    }

    for(; Environment::timestep <= Environment::getEndTimestep(); ++Environment::timestep) {
        if( Environment::isRootNode ) {
            cout << "====== Iteration " << Environment::timestep << " =====" << endl;
        }
        root_grid->mainLoop();

        time_counter->begin("plotData");
        if (Environment::plotPotential())    IO::writeDataInParallel(root_grid, Environment::timestep, "potential");
        if (Environment::plotRho())          IO::writeDataInParallel(root_grid, Environment::timestep, "rho");
        if (Environment::plotEfield())       IO::writeDataInParallel(root_grid, Environment::timestep, "efield");
        if (Environment::plotDensity())      IO::writeDataInParallel(root_grid, Environment::timestep, "density");
        if (Environment::plotEnergy())       IO::plotEnergy(root_grid, Environment::timestep);
        if (Environment::plotEnergyDist())   IO::plotParticleEnergyDistribution(root_grid->getParticles());
        if (Environment::plotVelocityDist()) IO::plotParticleVelocityDistribution(root_grid->getParticles());

        // 電磁計算時の追加プロット
        if (Environment::isEMMode()) {
            if (Environment::plotBfield()) IO::writeDataInParallel(root_grid, Environment::timestep, "bfield");
            if (Environment::plotCurrent()) {
                root_grid->updateReferenceCurrent();
                IO::writeDataInParallel(root_grid, Environment::timestep, "current");
            }
        }

        IO::plotObjectsData(root_grid);
        IO::plotValidParticleNumber(root_grid);
        time_counter->end();
    }

    Environment::saveInfo();
    root_grid->saveResumeData();

    if( Environment::isRootNode ) {
        cout << "===== End Main Loop =====" << endl;
    }

    time_counter->destroy();

    return 0;
}
#endif
