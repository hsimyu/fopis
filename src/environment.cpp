#include "global.hpp"
#include "environment.hpp"
#include "mpiw.hpp"

// static variables
int Environment::max_particle_num;
int Environment::num_of_particle_types;
int Environment::timestep;
double Environment::dx;
double Environment::dt;
int Environment::nx, Environment::ny, Environment::nz;
int Environment::proc_x, Environment::proc_y, Environment::proc_z;
int Environment::cell_x, Environment::cell_y, Environment::cell_z;
int Environment::max_iteration;
int Environment::plot_energy_dist_width, Environment::plot_velocity_dist_width;
int Environment::plot_potential_width, Environment::plot_rho_width;
int Environment::plot_efield_width, Environment::plot_bfield_width;
int Environment::plot_energy_width, Environment::plot_particle_width;
int Environment::plot_density_width;
bool Environment::isRootNode;
bool Environment::onLowXedge, Environment::onHighXedge;
bool Environment::onLowYedge, Environment::onHighYedge;
bool Environment::onLowZedge, Environment::onHighZedge;
std::string Environment::jobtype;
std::string Environment::solver_type;
std::string Environment::boundary;
std::string Environment::dimension;
ParticleType* Environment::ptype;

/*
* @params
* const std::string axis: "x", "y", or "z"
* const std::string low_or_up: "l" or "u"
*/
std::string Environment::getBoundaryCondition(const std::string axis, const std::string low_or_up) {
    int axisIndex;

    if (axis == "x") {
        axisIndex = 0;
    } else if (axis == "y") {
        axisIndex = 2;
    } else if (axis == "z") {
        axisIndex = 4;
    } else {
        throw std::invalid_argument( (format("Unknown axis type was passed. axis = %s, low_or_up = %s") % axis % low_or_up).str() );
    }

    std::string res;

    if (low_or_up == "l") {
        res = boundary.substr(axisIndex, 1);
    } else if (low_or_up == "u") {
        res = boundary.substr(axisIndex + 1, 1);
    } else {
        throw std::invalid_argument( (format("Unknown low_or_up type was passed. axis = %s, low_or_up = %s") % axis % low_or_up).str() );
    }

    return res;
}

void Environment::printInfo(void){
    cout << "[Environment]" << endl;
    cout << "      jobtype: " << jobtype << endl;
    cout << "     max_pnum: " << max_particle_num << endl;
    cout << "    iteration: " << max_iteration << endl;
    cout << "boundary cond: " << boundary << endl;
    cout << "           dx: " << (format("%8.2f") % dx).str() << "   m" << endl;
    cout << "           dt: " << (format("%6.2e") % dt).str() << " sec" << endl;
    cout << "   nx, ny, nz: " << format("%1%x%2%x%3%") % nx % ny % nz << " grids [total]" << endl;
    cout << "      process: " << format("%1%x%2%x%3%") % proc_x % proc_y % proc_z << " = " << (proc_x * proc_y * proc_z) << " procs" << endl;
    cout << "         cell: " << format("%1%x%2%x%3%") % cell_x % cell_y % cell_z << " grids [/proc] " << endl;
    cout << "      cell(+): " << format("%1%x%2%x%3%") % (cell_x + 2) % (cell_y + 2) % (cell_z + 2) << " grids [/proc] (with glue cells) " << endl << endl;

    cout << "[IO width]" << endl;
    cout << "  energy_dist: " << plot_energy_dist_width << endl;
    cout << "velocity_dist: " << plot_velocity_dist_width << endl;
    cout << "    potential: " << plot_potential_width << endl;
    cout << "          rho: " << plot_rho_width << endl;
    cout << "       efield: " << plot_efield_width << endl;
    cout << "       bfield: " << plot_bfield_width << endl;
    cout << "     particle: " << plot_particle_width << endl;
    cout << "       energy: " << plot_energy_width << endl;
}
