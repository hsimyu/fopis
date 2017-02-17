#ifndef __TDPIC_ENVIRONMENT_H_INCLUDED__
#define __TDPIC_ENVIRONMENT_H_INCLUDED__
#include "mpiw.hpp"
#include <string>

class ParticleType;

struct Environment {
    private:
        static bool isPlot(const std::string type) {
            if(type == "potential"     && plot_potential_width     != 0) return (timestep % plot_potential_width == 0);
            if(type == "rho"           && plot_rho_width           != 0) return (timestep % plot_rho_width == 0);
            if(type == "efield"        && plot_efield_width        != 0) return (timestep % plot_efield_width == 0);
            if(type == "bfield"        && plot_bfield_width        != 0) return (timestep % plot_bfield_width == 0);
            if(type == "density"       && plot_density_width       != 0) return (timestep % plot_density_width == 0);
            if(type == "particle"      && plot_particle_width      != 0) return (timestep % plot_particle_width == 0);
            if(type == "energy"        && plot_energy_width        != 0) return (timestep % plot_energy_width == 0);
            if(type == "energy_dist"   && plot_energy_dist_width   != 0) return (timestep % plot_energy_dist_width == 0);
            if(type == "velocity_dist" && plot_velocity_dist_width != 0) return (timestep % plot_velocity_dist_width == 0);
            return false;
        }

    public:
        static int max_particle_num;
        static int num_of_particle_types;
        static int timestep;
        static double dx;
        static double dt;
        static int nx, ny, nz;
        static int proc_x, proc_y, proc_z;
        static int cell_x, cell_y, cell_z;
        static int max_iteration;
        static int plot_energy_dist_width, plot_velocity_dist_width;
        static int plot_potential_width, plot_rho_width;
        static int plot_efield_width, plot_bfield_width;
        static int plot_particle_width, plot_energy_width;
        static int plot_density_width;
        static bool isRootNode;

        static bool onLowXedge, onHighXedge;
        static bool onLowYedge, onHighYedge;
        static bool onLowZedge, onHighZedge;

        static std::string jobtype;
        static std::string solver_type;
        static std::string boundary;
        static std::string dimension;

        static ParticleType* ptype;

        //! 中身はMPI::Environment::getRankStr()と同様
        static std::string rankStr(void) {
            return MPIw::Environment::rankStr();
        }

        static void printInfo(void);

        // plot timing
        static bool plotPotential(void)    { return isPlot("potential"); }
        static bool plotRho(void)          { return isPlot("rho"); }
        static bool plotEfield(void)       { return isPlot("efield"); }
        static bool plotBfield(void)       { return isPlot("bfield"); }
        static bool plotDensity(void)      { return isPlot("density"); }
        static bool plotParticle(void)     { return isPlot("particle"); }
        static bool plotEnergy(void)       { return isPlot("energy"); }
        static bool plotEnergyDist(void)   { return isPlot("energy_dist"); }
        static bool plotVelocityDist(void) { return isPlot("velocity_dist"); }
};
#endif
