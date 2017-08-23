#include "global.hpp"
#include "initialize.hpp"
#include "environment.hpp"
#include "utils.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "mpiw.hpp"
#include "normalizer.hpp"
#include <tdpic_configure.h>

namespace Initializer {
    void initTDPIC(Grid*& root_grid){
#ifndef BUILD_TEST
        // load parameter from json
        std::string filename = "input.json";
        auto inputs = Utils::readJSONFile(filename);
        Initializer::loadEnvironment(inputs);

        // EnvironmentにMPIw::Environment情報をセット
        Initializer::setMPIInfoToEnv();

        // 粒子情報セット
        Initializer::loadParticleType(inputs);
#else
        Initializer::setTestEnvirontment();
        Initializer::setMPIInfoToEnv();
#endif

        //! initialize normalizer
        //! normalizerのセットはGridの生成より先
        Normalizer::setLengthUnit(Environment::dx);
        Normalizer::setTimeUnit(Environment::dt);
        Normalizer::setMassUnit(me);
        Normalizer::setChargeUnit(e);

        if( Environment::isRootNode ) {
            cout << "---    [ TDPIC v" << TDPIC_VERSION << " ]      --" << endl;
            cout << "    Built on: " << TDPIC_DATE << endl;
            cout << "Git Revision: " << TDPIC_REVISION << endl << endl;

            Environment::printInfo();
            Environment::checkPlasmaInfo();
            Environment::checkCFLCondition();
            //! データ書き込み用ディレクトリを作成
            Utils::createDir("data");
        }

        if(Environment::jobtype == "new") {
            root_grid = Initializer::initializeGrid();
        } else {
            // load from files
        }

        if( Environment::isRootNode ) {
            cout << "--  End Initializing  --" << endl;
        }
    }

    void setTestEnvirontment() {
        // test input
        Environment::jobtype = "new";
        Environment::nx = 16;
        Environment::ny = 16;
        Environment::nz = 16;
        Environment::proc_x = 1;
        Environment::proc_y = 1;
        Environment::proc_z = 1;
        Environment::solver_type = "EM";
        Environment::boundary = "DDDDDD";
        Environment::dimension = "3D";
        Environment::dx = 0.1;
        Environment::dt = 1e-8;
        Environment::cell_x = Environment::nx/Environment::proc_x;
        Environment::cell_y = Environment::ny/Environment::proc_y;
        Environment::cell_z = Environment::nz/Environment::proc_z;

        Environment::max_particle_num = MAX_PARTICLE_NUM;
        Environment::num_of_particle_types = 2;
        Environment::ptype = new ParticleType[Environment::num_of_particle_types];

        Environment::ptype[0].setId(0);
        Environment::ptype[0].setName("Electron");
        Environment::ptype[0].setType("ambient");
        Environment::ptype[0].setMass(1.0);
        Environment::ptype[0].setCharge(-1.0);
        Environment::ptype[0].setTemperature(1.0);
        Environment::ptype[0].setDensity(1.0e6);
        Environment::ptype[0].setPcell(20);
        Environment::ptype[0].calcTotalNumber();
        Environment::ptype[0].calcSize();

        Environment::ptype[1].setId(1);
        Environment::ptype[1].setName("Proton");
        Environment::ptype[1].setType("ambient");
        Environment::ptype[1].setMass(1836.0 * 1.0);
        Environment::ptype[1].setCharge(1.0);
        Environment::ptype[1].setTemperature(1.0);
        Environment::ptype[1].setDensity(1.0e6);
        Environment::ptype[1].setPcell(20);
        Environment::ptype[1].calcTotalNumber();
        Environment::ptype[1].calcSize();
    }

    void setMPIInfoToEnv() {
        // for visiblity
        const int npux = Environment::proc_x;
        const int npuy = Environment::proc_y;
        const int npuz = Environment::proc_z;

        if(MPIw::Environment::numprocs != npux * npuy * npuz) {
            if(MPIw::Environment::rank == 0) {
                cout << format("[ERROR] Allocated Process Number [%d] is different from [%d] inputted from json.") % MPIw::Environment::numprocs % (npux * npuy * npuz) << endl;
            }
            MPIw::Environment::abort(1);
        }

        // Processの積み方は
        // x->y->z

        for(int k = 0; k < npuz; ++k){
            for(int j = 0; j < npuy; ++j){
                for(int i = 0; i < npux; ++i){
                    if( (i + npux * j + npux * npuy * k) == MPIw::Environment::rank ){
                        MPIw::Environment::xrank = i;
                        MPIw::Environment::yrank = j;
                        MPIw::Environment::zrank = k;

                        // set adjacent ranks
                        MPIw::Environment::adj[0] = ( ((npux + i-1) % npux) + npux*j + npux*npuy*k);
                        MPIw::Environment::adj[1] = ( ((i+1) % npux) + npux*j + npux*npuy*k);

                        MPIw::Environment::adj[2] = (i + npux*((npuy+j-1) % npuy) + npux*npuy*k);
                        MPIw::Environment::adj[3] = (i + npux*((j+1) % npuy) + npux*npuy*k);

                        MPIw::Environment::adj[4] = (i + npux*j + npux*npuy*((npuz+k-1) % npuz));
                        MPIw::Environment::adj[5] = (i + npux*j + npux*npuy*((k+1) % npuz));
                    }
                }
            }
        }

        Environment::onLowXedge = (MPIw::Environment::xrank == 0);
        Environment::onHighXedge = (MPIw::Environment::xrank == Environment::proc_x - 1);
        Environment::onLowYedge = (MPIw::Environment::yrank == 0);
        Environment::onHighYedge = (MPIw::Environment::yrank == Environment::proc_y - 1);
        Environment::onLowZedge = (MPIw::Environment::zrank == 0);
        Environment::onHighZedge = (MPIw::Environment::zrank == Environment::proc_z - 1);

        // 0のノードをルートとして扱う
        Environment::isRootNode = (MPIw::Environment::rank == 0);

        if(Environment::isRootNode) cout << format("    [MPIINFO] allocated processes: %d") % MPIw::Environment::numprocs << endl;
    }

    void loadEnvironment(picojson::object& inputs){
        auto env_inputs = inputs["Environment"].get<picojson::object>();

        //! 読み込まなくてよい部分
        Environment::max_particle_num = MAX_PARTICLE_NUM;
        Environment::timestep = 1;

        for (auto it = env_inputs.begin(); it != env_inputs.end(); ++it) {
            // string で switch したい...
            if(it->first == "nx"){
                Environment::nx = static_cast<int>(it->second.get<double>());
            } else if(it->first == "ny"){
                Environment::ny = static_cast<int>(it->second.get<double>());
            } else if(it->first == "nz"){
                Environment::nz = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_x"){
                Environment::proc_x = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_y"){
                Environment::proc_y = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_z"){
                Environment::proc_z = static_cast<int>(it->second.get<double>());
            } else if(it->first == "dt"){
                Environment::dt = it->second.get<double>();
            } else if(it->first == "dx"){
                Environment::dx = it->second.get<double>();
            } else if(it->first == "max_timestep"){
                Environment::max_timestep = static_cast<int>(it->second.get<double>());
            } else if(it->first == "job_type"){
                Environment::jobtype = it->second.to_str();
            } else if(it->first == "solver_type"){
                Environment::solver_type = it->second.to_str();
            } else if(it->first == "boundary"){
                Environment::boundary = it->second.to_str();
            } else if(it->first == "dimension"){
                Environment::dimension = it->second.to_str();
            } else {
                std::cout <<"Unsupportted Key [" << it->first << "] is in json." << std::endl;
            }
        }

        // 1プロセスあたりのグリッド数
        // これに2を加えた数がのりしろ分になる
        Environment::cell_x = Environment::nx/Environment::proc_x;
        Environment::cell_y = Environment::ny/Environment::proc_y;
        Environment::cell_z = Environment::nz/Environment::proc_z;

        auto io_inputs = inputs["IO"].get<picojson::object>();
        for(auto it = io_inputs.begin(); it != io_inputs.end(); ++it){
            if(it->first == "plot_energy_dist_width"){
                Environment::plot_energy_dist_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_velocity_dist_width"){
                Environment::plot_velocity_dist_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_potential_width"){
                Environment::plot_potential_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_rho_width"){
                Environment::plot_rho_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_efield_width"){
                Environment::plot_efield_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_bfield_width"){
                Environment::plot_bfield_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_density_width"){
                Environment::plot_density_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_particle_width"){
                Environment::plot_particle_width = static_cast<int>(it->second.get<double>());
            } else if(it->first == "plot_energy_width"){
                Environment::plot_energy_width = static_cast<int>(it->second.get<double>());
            } else {
                std::cout <<"Unsupportted Key [" << it->first << "] is in json." << std::endl;
            }
        }

        //! 物体情報
        auto object_inputs = inputs["Object"].get<picojson::object>();
        for(auto it = object_inputs.begin(); it != object_inputs.end(); ++it){
            const auto obj_name = it->first;
            auto obj_info = it->second.get<picojson::object>();

            Environment::ObjectInfo_t obj;
            obj.name = obj_name;
            obj.file_name = obj_info["file_name"].to_str();
            obj.surface_type = obj_info["surface_type"].to_str();
            obj.history_width = static_cast<unsigned int>(obj_info["history_width"].get<double>());
            
            Environment::objects_info.push_back( std::move(obj) );
        }
    }

    void loadParticleType(picojson::object& inputs){
        auto plasma_inputs = inputs["Plasma"].get<picojson::object>();
        int particle_types = 0;
        for(auto i = plasma_inputs.begin(); i != plasma_inputs.end(); ++i){
            ++particle_types;
        }
        Environment::num_of_particle_types = particle_types;
        ParticleType* ptype = new ParticleType[particle_types];
        int ii = 0;
        for(auto it = plasma_inputs.begin(); it != plasma_inputs.end(); ++it){
            std::string name = it->first;
            auto plasma = it->second.get<picojson::object>();
            ptype[ii].setId(ii);
            ptype[ii].setName(name);
            ptype[ii].setType( plasma["type"].to_str() );
            ptype[ii].setMass( plasma["mass"].get<double>() );
            ptype[ii].setCharge( plasma["charge"].get<double>() );
            ptype[ii].setTemperature( plasma["temperature"].get<double>() );
            ptype[ii].setDensity( plasma["density"].get<double>() );
            ptype[ii].setPcell( static_cast<int>((plasma["particle_per_cell"].get<double>() )) );
            ptype[ii].calcTotalNumber();
            ptype[ii].calcSize();

            ++ii;
        }

        //! @note: staticなポインタって持っても大丈夫？
        Environment::ptype = ptype;
    }

    Grid* initializeGrid(void){
        return new Grid();
    }
}

