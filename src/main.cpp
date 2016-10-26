#include <iostream>
#include <math.h>
#include "2dpic.h"
#include <Eigen/Core>

using std::cout;
using std::endl;
using Eigen::MatrixXd;

int main(int argc, char* argv[]){
    int nr = 20;
    double density = 1.0e6;
    double dx = 0.20;
    int num_particle = Initializer::getSizeOfSuperParticle(nr, density, dx);
    cout << "num: " << num_particle << endl;

    // Input parameter from the file
    int nx = 8, ny = 8;

    // Initialize particle and field
    //double phi[nx][ny];
    // auto phi = new double(nx, ny);
    MatrixXd phi(nx, ny), rho(nx, ny);
    cout << "-- phi --" << endl << phi << endl;
    cout << "-- rho --" << endl << rho << endl;

    // Utils::print2DArray<double>(phi, nx, ny);
    return 0;
}
