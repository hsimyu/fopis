/**
 * @mainpage
 * TDPIC: 三次元Electromagnetic Full Particle-in-cell code
 * 帯電解析用 & マルチグリッド解析用のPICコード
 **/

#include <tdpic.h>

using std::cout;
using std::endl;

#ifndef BUILD_TEST
int main(int argc, char* argv[]){
    cout << "---    [ TDPIC ]      --" << endl;
    cout << "-- Start Initializing --" << endl;

    std::string filename = "input.json";

    // load parameter from json
    auto inputs = Utils::readJSONFile(filename);
    Environment* env = Initializer::loadEnvironment(inputs);
    cout << env << endl;

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

    cout << "--  End Initializing  --" << endl;
    cout << "--  Begin A Loop  --" << endl;

    // particle -> space charge
    root_grid->updateRho(env);

    // space charge -> potential
    root_grid->getField()->solvePoisson(env);

    // potential -> efield
    root_grid->getField()->updateEfield(env);

    cout << "--  End A Loop  --" << endl;

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
    cout << "-- Bx --" << endl;
    IO::print3DArray( root_grid->getField()->getBx(), env->cell_x + 1, env->cell_y + 2, env->cell_z + 2);
    cout << "-- By --" << endl;
    IO::print3DArray( root_grid->getField()->getBy(), env->cell_x + 2, env->cell_y + 1, env->cell_z + 2);
    cout << "-- Bz --" << endl;
    IO::print3DArray( root_grid->getField()->getBz(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 1);
#endif

    return 0;
}
#endif
