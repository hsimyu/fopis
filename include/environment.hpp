#ifndef __TDPIC_ENVIRONMENT_H_INCLUDED__
#define __TDPIC_ENVIRONMENT_H_INCLUDED__
#include "mpiw.hpp"
#include "particle_type.hpp"
#include <array>

struct Environment {
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

        using ParticleTypePtr = std::shared_ptr<ParticleType>;
        using EmissionParticleTypePtr = std::shared_ptr<EmissionParticleType>;
        using AmbientParticlePtr = std::shared_ptr<AmbientParticleType>;
        using BeamParticlePtr = std::shared_ptr<BeamParticleType>;

        using AmbientParticleList = std::vector<AmbientParticlePtr>;
        using BeamParticleList = std::vector<BeamParticlePtr>;

        //! イテレータ経由で使う
        static AmbientParticleList getAmbientParticleTypes() { return ambient_particles; }
        static BeamParticleList getBeamParticleTypes() { return beam_particles; }

        //! インデックスを探して使う
        static ParticleTypePtr getParticleType(const int pid);
        static EmissionParticleTypePtr getEmissionParticleType(const int pid);
        static AmbientParticlePtr getAmbientParticleType(const int pid);
        static BeamParticlePtr getBeamParticleType(const int pid);

        //! 現在のプロセスが担当する領域の実座標（glueセルなし）を返す
        static int getAssignedXBegin(void) { return ::MPIw::Environment::xrank * cell_x; }
        static int getAssignedXEnd(void) { return (::MPIw::Environment::xrank + 1) * cell_x - 1; }
        static int getAssignedYBegin(void) { return ::MPIw::Environment::yrank * cell_y; }
        static int getAssignedYEnd(void) { return (::MPIw::Environment::yrank + 1) * cell_y - 1; }
        static int getAssignedZBegin(void) { return ::MPIw::Environment::zrank * cell_z; }
        static int getAssignedZEnd(void) { return (::MPIw::Environment::zrank + 1) * cell_z - 1; }

        static std::array<int, 3> getRelativePositionOnRootWithGlue(const int i, const int j, const int k) {
            //! Glueセルありで返す
            return std::array<int, 3>{{
                i - Environment::getAssignedXBegin() + 1,
                j - Environment::getAssignedYBegin() + 1,
                k - Environment::getAssignedZBegin() + 1
            }};
        }

        static std::array<double, 3> getRelativePositionOnRootWithoutGlue(const double x, const double y, const double z) {
            //! Glueセルなしで返す
            return std::array<double, 3>{{
                x - Environment::getAssignedXBegin(),
                y - Environment::getAssignedYBegin(),
                z - Environment::getAssignedZBegin()
            }};
        }

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
        static AmbientParticleList ambient_particles;
        static BeamParticleList beam_particles;
};
#endif
