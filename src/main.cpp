#include <iostream>
#include <math.h>
#include <boost/multi_array.hpp>
#include "tdpic.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){
    int nr = 20;
    double density = 1.0e6;
    const double dx = 0.20;
    int num_particle = Initializer::getSizeOfSuperParticle(nr, density, dx);
    cout << "num: " << num_particle << endl;

    // TODO: Input parameter from the file
    int nx = 8, ny = 8, nz = 1;

    // TODO: get domain decomposition size
    int proc_x = 1, proc_y = 1, proc_z = 1;
    int cell_x = nx/proc_x, cell_y = ny/proc_y, cell_z = nz/proc_z;

    // TODO: これは全体の長さなのであとで各領域ごとの長さにする
    FieldPointers field;
    Initializer::setFieldPointers(&field, cell_x, cell_y, cell_z);
    Utils::print3DArray( field.getRho() );

    // ParticleInfo pinfo;
    // Initializer::setParticleInfo(&pinfo);

    return 0;
}
