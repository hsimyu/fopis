#include "tdpic.h"

std::ostream& operator<<(std::ostream& ost, const Environment* env){
    ost << "[Environment]" << std::endl;
    ost << "         nx: " << env->nx << std::endl;
    ost << "         ny: " << env->ny << std::endl;
    ost << "         nz: " << env->nz << std::endl;
    ost << "     proc_x: " << env->proc_x << std::endl;
    ost << "     proc_y: " << env->proc_y << std::endl;
    ost << "     proc_z: " << env->proc_z << std::endl;
    ost << "     cell_x: " << env->cell_x - 2 << " ("<< env->cell_x << ")" << std::endl;
    ost << "     cell_y: " << env->cell_y - 2 << " ("<< env->cell_y << ")" << std::endl;
    ost << "     cell_z: " << env->cell_z - 2 << " ("<< env->cell_z << ")" << std::endl;
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Environment& env){
    ost << "[Environment]" << std::endl;
    ost << "         nx: " << env.nx << std::endl;
    ost << "         ny: " << env.ny << std::endl;
    ost << "         nz: " << env.nz << std::endl;
    ost << "     proc_x: " << env.proc_x << std::endl;
    ost << "     proc_y: " << env.proc_y << std::endl;
    ost << "     proc_z: " << env.proc_z << std::endl;
    ost << "     cell_x: " << env.cell_x << std::endl;
    ost << "     cell_y: " << env.cell_y << std::endl;
    ost << "     cell_z: " << env.cell_z << std::endl;
    return ost;
}
