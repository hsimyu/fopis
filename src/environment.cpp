#include <tdpic.h>
#include <boost/format.hpp>
using boost::format;

std::ostream& operator<<(std::ostream& ost, const Environment* env){
    ost << "[Environment]" << std::endl;
    ost << "    jobtype: " << env->jobtype << std::endl;
    ost << "  iteration: " << env->max_iteration << std::endl;
    ost << "         dx: " << (format("%8.2f") % env->dx).str() << "   m" << std::endl;
    ost << "         dt: " << (format("%6.2e") % env->dt).str() << " sec" << std::endl;
    ost << " nx, ny, nz: " << format("%1%x%2%x%3%") % env->nx % env->ny % env->nz << " grids [total]" << std::endl;
    ost << "    process: " << format("%1%x%2%x%3%") % env->proc_x % env->proc_y % env->proc_z << " = " << (env->proc_x * env->proc_y * env->proc_z) << " procs" << std::endl;
    ost << "       cell: " << format("%1%x%2%x%3%") % (env->cell_x - 2) % (env->cell_y - 2) % (env->cell_z - 2) << " grids [/proc] " << std::endl;
    ost << "    cell(+): " << format("%1%x%2%x%3%") % env->cell_x % env->cell_y % env->cell_z << " grids [/proc] (with glue cells) " << std::endl;
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Environment& env){
    ost << "[Environment]" << std::endl;
    ost << "    jobtype: " << env.jobtype << std::endl;
    ost << "  iteration: " << env.max_iteration << std::endl;
    ost << "         dx: " << (format("%8.2f") % env.dx).str() << "   m" << std::endl;
    ost << "         dt: " << (format("%6.2e") % env.dt).str() << " sec" << std::endl;
    ost << " nx, ny, nz: " << format("%1%x%2%x%3%") % env.nx % env.ny % env.nz << " grids [total]" << std::endl;
    ost << "    process: " << format("%1%x%2%x%3%") % env.proc_x % env.proc_y % env.proc_z << " = " << (env.proc_x * env.proc_y * env.proc_z) << " procs" << std::endl;
    ost << "       cell: " << format("%1%x%2%x%3%") % (env.cell_x - 2) % (env.cell_y - 2) % (env.cell_z - 2) << " grids [/proc] " << std::endl;
    ost << "    cell(+): " << format("%1%x%2%x%3%") % env.cell_x % env.cell_y % env.cell_z << " grids [/proc] (with glue cells) " << std::endl;
    return ost;
}
