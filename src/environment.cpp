#include "global.hpp"
#include "environment.hpp"
#include "particle.hpp"
#include "utils.hpp"
#include "normalizer.hpp"
#include "mpiw.hpp"

// static variables
int Environment::max_particle_num;
int Environment::num_of_particle_types;
int Environment::timestep;
int Environment::max_timestep;
double Environment::dx;
double Environment::dt;
int Environment::nx, Environment::ny, Environment::nz;
int Environment::proc_x, Environment::proc_y, Environment::proc_z;
int Environment::cell_x, Environment::cell_y, Environment::cell_z;
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
std::vector<ObjectInfo_t> Environment::objects_info;
std::vector<ParticleType> Environment::ptype;

/*
* @params
* const std::string axis: "x", "y", or "z"
* const AXIS_SIDE low_or_up: AXIS_SIDE::low or AXIS_SIDE::up
*/
std::string Environment::getBoundaryCondition(const AXIS axis, const AXIS_SIDE low_or_up) {
    int axisIndex = Utils::getAxisIndex(axis);

    // 列挙型AXIS_SIDEから対応するindexを取る
    int lowOrUpIndex = Utils::getLowOrUpIndex(low_or_up);

    // "DDPPDD"のような文字列が入力されるので
    // 軸 + 上下からindexを判別して文字列を抜き出す
    return boundary.substr(2 * axisIndex + lowOrUpIndex, 1);
}

bool Environment::isOnEdge(const AXIS axis, const AXIS_SIDE low_or_up) {
    switch(axis) {
        case AXIS::x:
            if (low_or_up == AXIS_SIDE::low) {
                return onLowXedge;
            } else {
                return onHighXedge;
            }
            break;
        case AXIS::y:
            if (low_or_up == AXIS_SIDE::low) {
                return onLowYedge;
            } else {
                return onHighYedge;
            }
            break;
        case AXIS::z:
            if (low_or_up == AXIS_SIDE::low) {
                return onLowZedge;
            } else {
                return onHighZedge;
            }
            break;
        default:
            break;
    }

    return false;
}

//! 現在の領域の上限/下限が、
//! 計算空間全体の境界でない or 計算空間全体の境界であるが周期境界である
//! かどうかをチェックする。
//! どちらかが満たされていれば、その方向の端の要素は境界とみなす必要がない (Iteration時の判定に用いる)
bool Environment::isNotBoundary(const AXIS axis, const AXIS_SIDE low_or_up) {
    return ( (!Environment::isOnEdge(axis, low_or_up)) || (Environment::getBoundaryCondition(axis, low_or_up) == "P" ) );
}

bool Environment::isBoundary(const AXIS axis, const AXIS_SIDE low_or_up) {
    return !isNotBoundary(axis, low_or_up);
}


void Environment::printInfo(void){
    cout << "[Environment]" << endl;
    cout << "      jobtype: " << jobtype << endl;
    cout << "     max pnum: " << max_particle_num << endl;
    cout << " max timestep: " << max_timestep << endl;
    cout << "boundary cond: " << boundary << endl;
    cout << "           dx: " << (format("%8.2f") % dx).str() << "   m" << endl;
    cout << "           dt: " << (format("%6.2e") % dt).str() << " sec" << endl;
    cout << "   nx, ny, nz: " << format("%1%x%2%x%3%") % nx % ny % nz << " grids [total]" << endl;
    cout << "      process: " << format("%1%x%2%x%3%") % proc_x % proc_y % proc_z << " = " << (proc_x * proc_y * proc_z) << " procs" << endl;
    cout << "         cell: " << format("%1%x%2%x%3%") % cell_x % cell_y % cell_z << " grids [/proc] " << endl;
    cout << "      cell(+): " << format("%1%x%2%x%3%") % (cell_x + 2) % (cell_y + 2) % (cell_z + 2) << " grids [/proc] (with glue cells) " << endl << endl;

    cout << "    [IO width]" << endl;
    cout << "      energy_dist: " << plot_energy_dist_width << endl;
    cout << "    velocity_dist: " << plot_velocity_dist_width << endl;
    cout << "        potential: " << plot_potential_width << endl;
    cout << "              rho: " << plot_rho_width << endl;
    cout << "           efield: " << plot_efield_width << endl;
    cout << "           bfield: " << plot_bfield_width << endl;
    cout << "         particle: " << plot_particle_width << endl;
    cout << "           energy: " << plot_energy_width << endl << endl;

    cout << "    [Objects]" << endl;
    for(const auto& object_info : objects_info) {
        cout << "        [" << object_info.name << "]" << endl;
        cout << "                file_name: "<< object_info.file_name << endl;
        cout << "             surface_type: "<< object_info.surface_type << endl;
        cout << "            history_width: "<< object_info.history_width << endl << endl;
    }
}

void Environment::checkCFLCondition(void) {
    cout << "[CFL Validation]" << endl;

    const auto courant = Normalizer::c;
    cout << "    Courant number: " << courant << endl;

    if ( courant > 1.0 ) {
        cout << "    [ERROR] Courant number exceeds 1.0. " << endl;
        MPIw::Environment::abort(1);
    }
}

void Environment::checkPlasmaInfo(void) {
    cout << "[Plasma Info]" << endl;

    double total_debye = 0.0;

    for (int pid = 0; pid < num_of_particle_types; pid++) {
        ParticleType& pt = ptype[pid];
        cout << pt;

        // プラズマ特徴量のチェック
        double debye = pt.calcDebyeLength();
        
        cout << "    [Debye Length]:" << endl;
        cout << "        Debye Length: " << debye << " m" << endl;
        cout << "        Debye / dx = "  << debye / dx << " > 1.0?: " << ((debye / dx > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
        cout << "        (nx * dx) / Debye = " << (nx * dx) / debye << " > 1.0?: "
            << ((((nx     * dx) / debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
        cout << "        (ny * dx) / Debye = " << (ny * dx) / debye << " > 1.0?: "
            << ((((ny     * dx) / debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
        cout << "        (nz * dx) / Debye = " << (nz * dx) / debye << " > 1.0?: "
            << ((((nz * dx) / debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl << endl;

        // 1/ld_total^2 = \sigma 1/ld_i^2
        if (pt.getType() == "ambient") total_debye += 1.0/pow(debye, 2);

        cout << "    [Thermal Velocity]" << endl;
        cout << "        Vth = " << pt.calcThermalVelocity() << " < 0.4?: "
            << (pt.calcThermalVelocity() < 1.0 ? "OK" : "*NOT SATISFIED*") << endl << endl;

        cout << "    [Plasma Frequency]" << endl;
        double omega_p = pt.calcPlasmaFrequency();
        cout << "        omega_p = " << omega_p << endl;
        cout << "        1 / omega_p = " << 1.0 / omega_p << endl;
        cout << "        1 / (omega_p * dt) = " << 1.0 / (omega_p * dt) << " > 5.0?: "
            << ( (1.0 / (omega_p * dt)) > 5.0 ? "OK" : "*NOT SATISFIED*") << endl;
        cout << "        1 / (omega_p * dt) < total_timestep?: "
            << ( (1.0 / (omega_p * dt)) < max_timestep ? "OK" : "*NOT SATISFIED*") << endl;
        cout << endl;
    }

    // 1/ld_total^2 -> ld_total
    total_debye = sqrt(1.0/total_debye);

    cout << "[Total Debye Length]:" << endl;
    cout << "Total Debye Length: " << total_debye << " m" <<  endl;
    cout << "    Total Debye / dx = "  << total_debye / dx << " > 1.0?: " << ((total_debye / dx > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
    cout << "    (nx * dx) / Total Debye = " << (nx * dx) / total_debye << " > 1.0?: "
            << ((((nx * dx) / total_debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
    cout << "    (ny * dx) / Total Debye = " << (ny * dx) / total_debye << " > 1.0?: "
            << ((((ny * dx) / total_debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl;
    cout << "    (nz * dx) / Total Debye = " << (nz * dx) / total_debye << " > 1.0?: "
            << ((((nz * dx) / total_debye) > 1.0) ? "OK" : "*NOT SATISFIED*") << endl << endl;
}
