#include <tdpic.h>

// static variables
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
    cout << "[Environment]" << std::endl;
    cout << "    jobtype: " << jobtype << std::endl;
    cout << "  iteration: " << max_iteration << std::endl;
    cout << "         dx: " << (format("%8.2f") % dx).str() << "   m" << std::endl;
    cout << "         dt: " << (format("%6.2e") % dt).str() << " sec" << std::endl;
    cout << " nx, ny, nz: " << format("%1%x%2%x%3%") % nx % ny % nz << " grids [total]" << std::endl;
    cout << "    process: " << format("%1%x%2%x%3%") % proc_x % proc_y % proc_z << " = " << (proc_x * proc_y * proc_z) << " procs" << std::endl;
    cout << "       cell: " << format("%1%x%2%x%3%") % cell_x % cell_y % cell_z << " grids [/proc] " << std::endl;
    cout << "    cell(+): " << format("%1%x%2%x%3%") % (cell_x + 2) % (cell_y + 2) % (cell_z + 2) << " grids [/proc] (with glue cells) " << std::endl;
}

std::string Environment::rankStr(void) {
    return MPIw::Environment::rankStr();
}
