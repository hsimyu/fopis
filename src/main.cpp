#include "tdpic.h"

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
    if(env->jobtype == "new") {
        ptype = Initializer::loadParticleType(inputs, env);
        field = Initializer::initializeField(env);
    } else {
        // load ptype?
        // load particle
        // load grid
        // load field
    }
    Utils::print3DArray( field->getRho() );

    // std::unique_ptr<Particle[]> particles(new Particle[ptype.getTotalNumber()]);
    // Initializer::setParticleInfo(&pinfo);

    cout << "--  End Initializing  --" << endl;
    return 0;
}
