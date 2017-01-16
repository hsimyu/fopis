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

    // load parameter from json
    std::string filename = "input.json";
    auto inputs = Utils::readJSONFile(filename);
    Environment* env = Initializer::loadEnvironment(inputs);

    // EnvironmentにMPI::Environment情報をセット
    Initializer::setMPIInfoToEnv(env);

    if( env->isRootNode ) {
        cout << "---    [ TDPIC ]      --" << endl;
        cout << env << endl;

        //! データ書き込み用ディレクトリを作成
        Utils::createDir("data");
    }

    // initialize normalizer
    Utils::Normalizer::x_unit = env->dx;
    Utils::Normalizer::t_unit = env->dt;
    Utils::Normalizer::e_unit = e;

    ParticleType* ptype;
    Grid* root_grid;

    if(env->jobtype == "new") {
        ptype = Initializer::loadParticleType(inputs, env);
        root_grid = Initializer::initializeGrid(env);

        // generate root field
        Initializer::initializeRootField(env, root_grid);
    } else {
        // load ptype?
        // load particle
        // load grid
        // load field
    }

    if( env->isRootNode ) {
        cout << "--  End Initializing  --" << endl;
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

        // IO::writeData(root_grid, 0);
#endif

        // Level 2まで
        // patchesは4個?
        root_grid->makeChild(2, 2, 2, 5, 5, 5);
        root_grid->makeChild(4, 8, 8, 9, 9, 9);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);
        root_grid->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);
        root_grid->getChildren()[0]->getChildren()[0]->makeChild(2, 2, 2, 3, 3, 3);
        cout << root_grid << endl;

        int* test = root_grid->getChildOfPatches();
        cout << root_grid->getNumOfPatches() << endl;
    }

    IO::writeDataInParallel(root_grid, 0, "potential");

    return 0;
}
#endif
