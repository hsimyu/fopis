#ifndef __TDPIC_ENVIRONMENT_H_INCLUDED__
#define __TDPIC_ENVIRONMENT_H_INCLUDED__
#include "mpiw.hpp"
#include "particle_type.hpp"

struct Environment {
    private:
        static bool isPlotTimestep(const std::string type) {
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

        //! 内部的には各粒子の種類毎にリストを持つ
        static std::vector<std::shared_ptr<AmbientParticleType>> ambient_particles;
        static std::vector<std::shared_ptr<BeamParticleType>> beam_particles;

    public:
        static int max_particle_num;
        static int num_of_particle_types;
        static int timestep;
        static int max_timestep;
        static double dx;
        static double dt;
        static int nx, ny, nz;
        static int proc_x, proc_y, proc_z;
        static int cell_x, cell_y, cell_z;
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

        //! 物体情報
        static std::vector<ObjectInfo_t> objects_info;

        //! 中身はMPI::Environment::getRankStr()と同様
        static std::string rankStr(void) {
            return MPIw::Environment::rankStr();
        }

        static int getSpecifiedRankFromXYZRanks(const int xrank, const int yrank, const int zrank) {
            return (xrank + proc_x * yrank + proc_x * proc_y * zrank);
        }

        static double getDataTime() { return static_cast<double>(timestep) * dt; }

        static void addAmbientParticleType(const std::shared_ptr<AmbientParticleType>& ptr) {
            ambient_particles.push_back(ptr);
            ++num_of_particle_types;
        }
        static int getNumOfAmbientParticles() {return ambient_particles.size();}

        static void addBeamParticleType(const std::shared_ptr<BeamParticleType>& ptr) {
            beam_particles.push_back(ptr);
            ++num_of_particle_types;
        }
        static int getNumOfBeamParticles() {return beam_particles.size();}

        static std::shared_ptr<ParticleType> getParticleType(const int pid) {
            if (pid < ambient_particles.size()) {
                return std::static_pointer_cast<ParticleType>(ambient_particles[pid]);
            } else {
                return std::static_pointer_cast<ParticleType>(beam_particles[pid - ambient_particles.size()]);
            }
        }

        static std::shared_ptr<AmbientParticleType> getAmbientParticleType(const int pid) {
            if (pid < ambient_particles.size()) {
                return ambient_particles[pid];
            } else {
                throw std::invalid_argument("[ERROR] The particle id passed to getAmbientParticleType() exceeds the ambient particle count.");
            }
        }

        static std::shared_ptr<BeamParticleType> getBeamParticleType(const int pid) {
            const auto ambient_size = ambient_particles.size();
            if (pid >= ambient_size && pid - ambient_size < beam_particles.size()) {
                return beam_particles[pid - ambient_size];
            } else {
                throw std::invalid_argument("[ERROR] The particle id passed to getBeamParticleType() exceeds the beam particle count.");
            }
        }

        //! 現在のプロセスが担当する領域の実座標（glueセルなし）を返す
        static int getAssignedXBegin(void) { return ::MPIw::Environment::xrank * cell_x; }
        static int getAssignedXEnd(void) { return (::MPIw::Environment::xrank + 1) * cell_x - 1; }
        static int getAssignedYBegin(void) { return ::MPIw::Environment::yrank * cell_y; }
        static int getAssignedYEnd(void) { return (::MPIw::Environment::yrank + 1) * cell_y - 1; }
        static int getAssignedZBegin(void) { return ::MPIw::Environment::zrank * cell_z; }
        static int getAssignedZEnd(void) { return (::MPIw::Environment::zrank + 1) * cell_z - 1; }

        //! 各方向の境界条件を取得する
        static std::string getBoundaryCondition(const AXIS, const AXIS_SIDE);

        //! 対応するEdgeに乗っているかどうかを返す
        static bool isOnEdge(const AXIS, const AXIS_SIDE);

        //! 現在の領域の上限/下限が、計算空間全体の境界でないか、計算空間全体の境界だが周期境界かどうかをチェックする
        static bool isNotBoundary(const AXIS, const AXIS_SIDE);
        static bool isBoundary(const AXIS, const AXIS_SIDE);

        static void printInfo(void);
        static void checkPlasmaInfo(void);
        static void checkCFLCondition(void);

        // plot timing
        static bool plotPotential(void)    { return isPlotTimestep("potential"); }
        static bool plotRho(void)          { return isPlotTimestep("rho"); }
        static bool plotEfield(void)       { return isPlotTimestep("efield"); }
        static bool plotBfield(void)       { return isPlotTimestep("bfield"); }
        static bool plotDensity(void)      { return isPlotTimestep("density"); }
        static bool plotParticle(void)     { return isPlotTimestep("particle"); }
        static bool plotEnergy(void)       { return isPlotTimestep("energy"); }
        static bool plotEnergyDist(void)   { return isPlotTimestep("energy_dist"); }
        static bool plotVelocityDist(void) { return isPlotTimestep("velocity_dist"); }
};
#endif