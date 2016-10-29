#include "tdpic.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){
    cout << "---    [ TDPIC ]      --" << endl;
    cout << "-- Start Initializing --" << endl;

    int nr = 1;
    double density = 1.0e6;
    const double dx = 0.20;

    // TODO: Input parameter from the file
    int nx = 8, ny = 8, nz = 1;
    //auto inputs = Utils::readInputFile("input.json");
    //cout << inputs.get<int>("Image.Width") << endl;

    int proc_x = 1, proc_y = 1, proc_z = 1;
    // TODO: get domain decomposition size
    // TODO: これは全体の長さなのであとで各領域ごとの長さにする
    int cell_x = nx/proc_x, cell_y = ny/proc_y, cell_z = nz/proc_z;

    FieldPointers field;
    Initializer::setFieldPointers(&field, cell_x, cell_y, cell_z);
    Utils::print3DArray( field.getRho() );

    std::string name_array[3] = {"electron", "proton", "beam"};
    ParticleType ptype[3];
    for(int i=0; i<3; i++){
        ptype[i].setId(i);
        ptype[i].setName(name_array[i]);
        ptype[i].calcTotalNumber(cell_x, cell_y, cell_z, nr);
        cout << ptype[i];
        Utils::printTotalMemory(ptype[i]);
    }

    // std::unique_ptr<Particle[]> particles(new Particle[ptype.getTotalNumber()]);
    // Initializer::setParticleInfo(&pinfo);

    cout << "--  End Initializing  --" << endl;
    return 0;
}
