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

        // particle -> space charge
        root_grid->updateRho(env);
    } else {
        // load ptype?
        // load particle
        // load grid
        // load field
    }

#ifdef DEBUG
    IO::print3DArray( root_grid->getField()->getRho(), env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
    IO::outputParticlePositions( env, root_grid->particles );
#endif

    cout << "--  End Initializing  --" << endl;
    return 0;
}
