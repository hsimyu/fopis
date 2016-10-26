#include <iostream>
#include <math.h>
#include <boost/multi_array.hpp>
#include "2dpic.h"

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
    // threeD_array ex(boost::extents[cell_x][cell_y][cell_z]);
    // threeD_array phi(boost::extents[cell_x][cell_y][cell_z]);
    // threeD_array phi(boost::extents[cell_x][cell_y][cell_z]);

    // // global variables
    // MatrixXd phi(cell_x, cell_y, cell_z), rho(cell_x, cell_y, cell_z);
    // MatrixXd ex(cell_x - 1, cell_y, cell_z), ey(cell_x, cell_y - 1, cell_z), ez(cell_x, cell_y, cell_z - 1);
    // MatrixXd bx(cell_x, cell_y - 1, cell_z - 1), by(cell_x - 1, cell_y, cell_z - 1), bz(cell_x - 1, cell_y - 1, cell_z);
    // phi = MatrixXd::Zero(cell_x, cell_y, cell_z);
    // rho = MatrixXd::Zero(cell_x, cell_y, cell_z);

    // ex = MatrixXd::Zero(cell_x - 1, cell_y, cell_z);
    // ey = MatrixXd::Zero(cell_x, cell_y - 1, cell_z);
    // ez = MatrixXd::Zero(cell_x, cell_y, cell_z - 1);

    // bx = MatrixXd::Zero(cell_x, cell_y - 1, cell_z - 1);
    // by = MatrixXd::Zero(cell_x - 1, cell_y, cell_z - 1);
    // bz = MatrixXd::Zero(cell_x - 1, cell_y - 1, cell_z);
    // cout << "-- phi --" << endl << phi << endl;
    // cout << "-- rho --" << endl << rho << endl;
    //(*field.getRho)[0][0][0] = 1.0;
    Utils::print3DArray( field.getPhi() );
    //Utils::print3DArray( field.getRho() );
    // Utils::print3DArray(rho);

    return 0;
}
