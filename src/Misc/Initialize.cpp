#include "global.hpp"
#include "initialize.hpp"
#include "environment.hpp"
#include "utils.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "mpiw.hpp"
#include "normalizer.hpp"
#include <tdpic_configure.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Initializer {
    std::shared_ptr<RootGrid> initTDPIC(const std::string& input_filename) {
        // load parameter from json
        auto inputs = Utils::readJSONFile(input_filename);
        Initializer::loadEnvironment(inputs);

        // EnvironmentにMPIw::Environment情報をセット
        Initializer::setMPIInfoToEnv();

        // EnviromentにOpenMPスレッド数をセット
        #ifdef _OPENMP
        Environment::num_threads = omp_get_max_threads();
        #else
        Environment::num_threads = 1;
        #endif

        //! initialize normalizer
        //! normalizerのセットはGridの生成より先
        Normalizer::setLengthUnit(Environment::dx);
        Normalizer::setTimeUnit(Environment::dt);
        Normalizer::setMassUnit(me);
        Normalizer::setChargeUnit(e);

        //! Static Field読み込み
        loadStaticField(inputs);

        // 粒子情報セット
        Initializer::loadParticleType(inputs);

        //! 継続計算の場合はTimestepをロード
        if (Environment::jobtype == "load") {
            Environment::loadInfo();
        }

        if (Environment::isRootNode) {
            cout << "---    [ TDPIC " << TDPIC_VERSION << " ]     --" << "\n";
            cout << "  Built date: " << TDPIC_DATE << "\n";
            cout << "  Git Revision: " << TDPIC_REVISION << "\n";
            cout << "  MPI processes: " << format("%d") % MPIw::Environment::numprocs << "\n";
            cout << "  OpenMP processes: " << format("%d") % Environment::num_threads << "\n" << endl;

            Environment::printInfo();
            Environment::checkPlasmaInfo();

            //! EMの場合はFDTD用の安定性条件をチェック
            if (Environment::isEMMode()) {
                Environment::checkCFLCondition();
            }

            //! データ書き込み用ディレクトリを作成
            Utils::createDir("data");
            Utils::createDir("data/raw_data");

            //! 継続計算データ用ディレクトリを作成
            Utils::createDir("resume");
            Utils::createDir("resume/objects");
        }

        std::shared_ptr<RootGrid> root_grid;

        if (Environment::jobtype == "new") {
            root_grid = std::make_shared<RootGrid>();
            if (Environment::isRootNode) {
                cout << "--  New Data Initializing --" << endl;
            }
            root_grid->initialize();
        } else if (Environment::jobtype == "load") {
            root_grid = std::make_shared<RootGrid>();

            if (Environment::isRootNode) {
                cout << "--  Resume Data (t = " << Environment::timestep << ") Loading --" << endl;
            }

            root_grid->loadResumeData();

            if (Environment::isRootNode) {
                cout << "--  Resume Data Initializing --" << endl;
            }

            root_grid->initializeForLoad();
        } else {
            if (Environment::isRootNode) {
                cout << "[ERROR] Unknown jobtype, neither ``new'' or ``load'', was inputted." << endl;
            }
            MPIw::Environment::abort(1);
        }

        if (Environment::isRootNode) {
            cout << "--  End Initializing  --" << endl;
        }

        return root_grid;
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

        Environment::setOnEdge(AXIS::x, AXIS_SIDE::low, (MPIw::Environment::xrank == 0));
        Environment::setOnEdge(AXIS::x, AXIS_SIDE::up, (MPIw::Environment::xrank == Environment::proc_x - 1));
        Environment::setOnEdge(AXIS::y, AXIS_SIDE::low, (MPIw::Environment::yrank == 0));
        Environment::setOnEdge(AXIS::y, AXIS_SIDE::up, (MPIw::Environment::yrank == Environment::proc_y - 1));
        Environment::setOnEdge(AXIS::z, AXIS_SIDE::low, (MPIw::Environment::zrank == 0));
        Environment::setOnEdge(AXIS::z, AXIS_SIDE::up, (MPIw::Environment::zrank == Environment::proc_z - 1));

        // 0のノードをルートとして扱う
        Environment::isRootNode = (MPIw::Environment::rank == 0);
    }

    void loadEnvironment(picojson::object& inputs){
        //! 読み込まなくてよい部分
        Environment::initial_timestep = 0;
        Environment::timestep = 1;

        //! Environment変数読み込み
        {
            auto env_inputs = inputs["Environment"].get<picojson::object>();
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
        }

        //! オプション
        {
            auto options_inputs = inputs["Options"].get<picojson::object>();
            for(auto it = options_inputs.begin(); it != options_inputs.end(); ++it){
                if (it->first == "use_existing_capacity_matrix") {
                    Environment::getOptions().setUseExistingCapacityMatrix(it->second.get<bool>());
                } else if (it->first == "maximum_poisson_post_loop") {
                    Environment::getOptions().setMaximumPoissonPostLoop( static_cast<unsigned int>(it->second.get<double>()) );
                } else if (it->first == "maximum_poisson_pre_loop") {
                    Environment::getOptions().setMaximumPoissonPreLoop( static_cast<unsigned int>(it->second.get<double>()) );
                }
            }
        }

        // 1プロセスあたりのグリッド数
        // これに2を加えた数がのりしろ分になる
        Environment::cell_x = Environment::nx/Environment::proc_x;
        Environment::cell_y = Environment::ny/Environment::proc_y;
        Environment::cell_z = Environment::nz/Environment::proc_z;

        //! IO関連
        {
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
                } else if(it->first == "plot_current_width"){
                    Environment::plot_current_width = static_cast<int>(it->second.get<double>());
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
        }

        //! 物体情報
        {
            auto object_inputs = inputs["Object"].get<picojson::object>();
            for(auto it = object_inputs.begin(); it != object_inputs.end(); ++it){
                const auto obj_name = it->first;
                auto obj_info = it->second.get<picojson::object>();

                ObjectInfo_t obj;
                obj.name = obj_name;
                for(auto it = obj_info.begin(); it != obj_info.end(); ++it) {

                    if(it->first == "file_name") {
                        obj.file_name = it->second.to_str();
                    } else if (it->first == "surface_type") {
                        obj.surface_type = it->second.to_str();
                    } else if (it->first == "plot_potential_mapping_width") {
                        obj.plot_potential_mapping_width = static_cast<unsigned int>(it->second.get<double>());
                    } else if (it->first == "is_potential_fixed") {
                        obj.is_potential_fixed = it->second.get<bool>();
                    } else if (it->first == "fixed_potential") {
                        obj.fixed_potential = it->second.get<double>();
                    } else if (it->first == "initial_potential_offset") {
                        obj.initial_potential_offset = it->second.get<double>();
                    } else if (it->first == "emit_particles") {
                        auto particle_names = it->second.get<picojson::object>();

                        for(auto pit = particle_names.begin(); pit != particle_names.end(); ++pit) {
                            const std::string pname = pit->first;

                            ParticleEmissionInfo pinfo;

                            auto inner = pit->second.get<picojson::object>();
                            for(auto inner_pit = inner.begin(); inner_pit != inner.end(); ++inner_pit) {
                                if (inner_pit->first == "emission_position") {
                                    pinfo.emission_position = Utils::convertPicoJSONArrayToVectorDouble( inner_pit->second.get<picojson::array>() );
                                } else if (inner_pit->first == "emission_vector") {
                                    pinfo.emission_vector = Utils::convertPicoJSONArrayToVectorDouble( inner_pit->second.get<picojson::array>() );
                                }
                            }
                            obj.emit_particle_info[pname] = pinfo;
                        }
                    } else if (it->first == "materials") {
                        auto material_names = it->second.get<picojson::object>();

                        for(auto mit = material_names.begin(); mit != material_names.end(); ++mit) {
                            obj.materials[ std::stoi(mit->first) ] = mit->second.to_str();
                        }
                    }
                }
                Environment::objects_info.push_back( std::move(obj) );
            }
        }
    }

    void loadStaticField(picojson::object& inputs) {
        if (inputs.count("Field") > 0) {
            auto field_inputs = inputs["Field"].get<picojson::object>();
            auto& static_field = Environment::getStaticField();

            for(auto it = field_inputs.begin(); it != field_inputs.end(); ++it){
                if (it->first == "static_bfield") {
                    auto static_bfield = Utils::convertPicoJSONArrayToVectorDouble( it->second.get<picojson::array>() );
                    const auto bfield_norm = Normalizer::normalizeBfield(1.0);
                    constexpr double bfield_input_unit = 1e-9; // nT
                    static_bfield[0] *= bfield_norm * bfield_input_unit;
                    static_bfield[1] *= bfield_norm * bfield_input_unit;
                    static_bfield[2] *= bfield_norm * bfield_input_unit;
                    static_field.setStaticBfield(static_bfield);
                } else if (it->first == "shine_vector") {
                    static_field.setShineVector(
                        Utils::convertPicoJSONArrayToVectorDouble( it->second.get<picojson::array>() )
                    );
                }
            }
        }
    }

    void loadParticleType(picojson::object& inputs) {
        auto plasma_inputs = inputs["Plasma"].get<picojson::object>();
        std::vector<ParticleType*> ptype;

        int id = 0;
        for(auto it = plasma_inputs.begin(); it != plasma_inputs.end(); ++it) {
            std::string name = it->first;
            auto plasma = it->second.get<picojson::object>();
            const auto& type = plasma["type"].to_str();

            if (type == "ambient") {
                auto ambient = std::make_shared<AmbientParticleType>();
                ambient->setId(id);
                ambient->setName(name);
                ambient->setType(type);
                ambient->setMass(plasma["mass"].get<double>());
                ambient->setCharge(plasma["charge"].get<double>());
                ambient->setTemperature(plasma["temperature"].get<double>());
                ambient->setDensity(plasma["density"].get<double>());

                if (plasma.count("drift_velocity") > 0) {
                    auto drift_velocity = Utils::convertPicoJSONArrayToVectorDouble( plasma["drift_velocity"].get<picojson::array>());

                    //! ドリフト速度を正規化 / 1e3かけてkm/s単位からm/s単位に変換
                    drift_velocity[0] *= 1e3 * Normalizer::normalizeVelocity(1.0);
                    drift_velocity[1] *= 1e3 * Normalizer::normalizeVelocity(1.0);
                    drift_velocity[2] *= 1e3 * Normalizer::normalizeVelocity(1.0);

                    ambient->setDriftVelocity(drift_velocity);
                }

                if (plasma.count("injection_axis") > 0) {
                    ambient->setInjectionAxis(
                        Utils::convertPicoJSONArrayToVectorBool( plasma["injection_axis"].get<picojson::array>())
                    );
                }

                ambient->setPcell(static_cast<int>((plasma["particle_per_cell"].get<double>())));
                ambient->updateSize();
                Environment::addAmbientParticleType(ambient);
            } else if (type == "beam") {
                auto beam = std::make_shared<BeamParticleType>();
                beam->setId(id);
                beam->setName(name);
                beam->setType(type);
                beam->setMass(plasma["mass"].get<double>());
                beam->setCharge(plasma["charge"].get<double>());
                beam->setTemperature(plasma["temperature"].get<double>());
                beam->setDensity(plasma["density"].get<double>());
                beam->setPcell(static_cast<int>((plasma["particle_per_cell"].get<double>())));
                beam->updateSize();

                beam->setAcceleratingPotential( Normalizer::normalizePotential(plasma["accel_potential"].get<double>()) );
                beam->setBeamCurrent( Normalizer::normalizeCurrent(plasma["beam_current"].get<double>()) );
                beam->setBeamDivergence( plasma["beam_divergence"].get<double>() );
                beam->setEmissionRadius( plasma["emission_radius"].get<double>() );
                beam->setEmissionType( plasma["emission_type"].to_str() );
                Environment::addBeamParticleType(beam);
            } else if (type == "photoelectron") {
                auto photo = std::make_shared<PhotoElectronParticleType>();
                photo->setId(id);
                photo->setName(name);
                photo->setType(type);
                photo->setMass(plasma["mass"].get<double>());
                photo->setCharge(plasma["charge"].get<double>());
                photo->setTemperature(plasma["temperature"].get<double>());
                photo->setDensity(plasma["density"].get<double>());
                photo->setPcell(static_cast<int>((plasma["particle_per_cell"].get<double>())));
                photo->updateSize();

                photo->setCurrentDensity(Normalizer::normalizeCurrentDensity(plasma["current_density"].get<double>()));
                Environment::addPhotoElectronParticleType(photo);
            } else if (type == "secondary") {
                auto see = std::make_shared<SecondaryParticleType>();
                see->setId(id);
                see->setName(name);
                see->setType(type);
                see->setMass(plasma["mass"].get<double>());
                see->setCharge(plasma["charge"].get<double>());
                see->setTemperature(plasma["temperature"].get<double>());
                see->setDensity(plasma["density"].get<double>());
                see->setPcell(static_cast<int>((plasma["particle_per_cell"].get<double>())));
                see->updateSize();

                Environment::addSecondaryParticleType(see);
            }

            ++id;
        }

        //! 粒子種の並び順をリセットする
        Environment::resetParticleTypeOrder();
    }
}

