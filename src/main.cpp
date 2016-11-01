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

    ParticleType* ptype = Initializer::loadParticleType(inputs, env);

    Field* field = Initializer::initializeField(env);
    Utils::print3DArray( field->getRho() );

    // std::unique_ptr<Particle[]> particles(new Particle[ptype.getTotalNumber()]);
    // Initializer::setParticleInfo(&pinfo);

    cout << "--  End Initializing  --" << endl;
    return 0;
}
