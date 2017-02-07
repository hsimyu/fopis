#include "global.hpp"
#include "environment.hpp"
#include "mpiw.hpp"

// static variables
int Environment::max_particle_num;
int Environment::num_of_particle_types;
double Environment::dx;
double Environment::dt;
int Environment::nx, Environment::ny, Environment::nz;
int Environment::proc_x, Environment::proc_y, Environment::proc_z;
int Environment::cell_x, Environment::cell_y, Environment::cell_z;
int Environment::max_iteration;
bool Environment::isRootNode;
bool Environment::onLowXedge, Environment::onHighXedge;
bool Environment::onLowYedge, Environment::onHighYedge;
bool Environment::onLowZedge, Environment::onHighZedge;
std::string Environment::jobtype;
std::string Environment::solver_type;
std::string Environment::boundary;
std::string Environment::dimension;
ParticleType* Environment::ptype;

void Environment::printInfo(void){
    cout << "[Environment]" << endl;
    cout << "    jobtype: " << jobtype << endl;
    cout << "   max_pnum: " << max_particle_num << endl;
    cout << "  iteration: " << max_iteration << endl;
    cout << "         dx: " << (format("%8.2f") % dx).str() << "   m" << endl;
    cout << "         dt: " << (format("%6.2e") % dt).str() << " sec" << endl;
    cout << " nx, ny, nz: " << format("%1%x%2%x%3%") % nx % ny % nz << " grids [total]" << endl;
    cout << "    process: " << format("%1%x%2%x%3%") % proc_x % proc_y % proc_z << " = " << (proc_x * proc_y * proc_z) << " procs" << endl;
    cout << "       cell: " << format("%1%x%2%x%3%") % cell_x % cell_y % cell_z << " grids [/proc] " << endl;
    cout << "    cell(+): " << format("%1%x%2%x%3%") % (cell_x + 2) % (cell_y + 2) % (cell_z + 2) << " grids [/proc] (with glue cells) " << endl;
}
