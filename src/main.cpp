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
    Field* field;
    Grid* root_grid;
    if(env->jobtype == "new") {
        ptype = Initializer::loadParticleType(inputs, env);
        field = Initializer::initializeField(env);
        root_grid = Initializer::initializeGrid(env, ptype);
    } else {
        // load ptype?
        // load particle
        // load grid
        // load field
    }

#ifdef DEBUG
    Utils::print3DArray( field->getRho() );
#endif

    cout << "--  End Initializing  --" << endl;
    return 0;
}
