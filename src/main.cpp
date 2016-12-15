/**
 * @mainpage
 * TDPIC: 三次元Electromagnetic Full Particle-in-cell code
 * 帯電解析用 & マルチグリッド解析用のPICコード
 **/

#include <tdpic.h>

using std::cout;
using std::endl;

int main(int argc, char* argv[]){
    cout << "---    [ TDPIC ]      --" << endl;
    cout << "-- Start Initializing --" << endl;

    std::string filename = "input.json";

    // load parameter from json
    auto inputs = Utils::readJSONFile(filename);
    Environment* env = Initializer::loadEnvironment(inputs);
    cout << env << endl;

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

    cout << "--  End A Loop  --" << endl;

#ifdef DEBUG
    IO::print3DArray( root_grid->getField()->getRho(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
    IO::print3DArray( root_grid->getField()->getPhi(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
    IO::outputParticlePositions( env, root_grid->particles );
#endif

    return 0;
}
