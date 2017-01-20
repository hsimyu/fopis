#include <math.h>
#include <tdpic.h>

namespace Initializer {
    void initTDPIC(Environment*& env, ParticleType*& ptype, Grid*& root_grid){
#ifndef BUILD_TEST
        // load parameter from json
        std::string filename = "input.json";
        auto inputs = Utils::readJSONFile(filename);
        env = Initializer::loadEnvironment(inputs);
        ptype = Initializer::loadParticleType(inputs, env);
#else
        Initializer::setTestEnvirontment(env, ptype);
#endif

        // EnvironmentにMPI::Environment情報をセット
        Initializer::setMPIInfoToEnv(env);

        if( env->isRootNode ) {
            cout << "---    [ TDPIC ]      --" << endl;
            cout << env << endl;
            //! データ書き込み用ディレクトリを作成
            Utils::createDir("data");
        }

        if(env->jobtype == "new") {
            root_grid = Initializer::initializeGrid(env);
        } else {
            // load ptype?
            // load particle
            // load grid
            // load field
        }

        // initialize normalizer
        Utils::Normalizer::x_unit = env->dx;
        Utils::Normalizer::t_unit = env->dt;
        Utils::Normalizer::e_unit = e;

        if( env->isRootNode ) {
           cout << "--  End Initializing  --" << endl;
        }
    }

    void setTestEnvirontment(Environment*& env, ParticleType*& ptype) {
        env = new Environment;

        // test input
        env->jobtype = "new";
        env->nx = 16;
        env->ny = 16;
        env->nz = 16;
        env->proc_x = 1;
        env->proc_y = 1;
        env->proc_z = 1;
        env->solver_type = "EM";
        env->boundary = "DDDDDD";
        env->dimension = "3D";
        env->dx = 0.1;
        env->dt = 1e-8;
        env->cell_x = env->nx/env->proc_x;
        env->cell_y = env->ny/env->proc_y;
        env->cell_z = env->nz/env->proc_z;

        env->num_of_particle_types = 2;
        ptype = new ParticleType[env->num_of_particle_types];

        ptype[0].setId(0);
        ptype[0].setName("Electron");
        ptype[0].setType("ambient");
        ptype[0].setMass(1.0);
        ptype[0].setCharge(-1.0);
        ptype[0].setTemperature(1.0);
        ptype[0].setDensity(1.0e6);
        ptype[0].setPcell(20);
        ptype[0].calcTotalNumber(env);
        ptype[0].calcSize(env);

        ptype[1].setId(1);
        ptype[1].setName("Proton");
        ptype[1].setType("ambient");
        ptype[1].setMass(1836.0 * 1.0);
        ptype[1].setCharge(1.0);
        ptype[1].setTemperature(1.0);
        ptype[1].setDensity(1.0e6);
        ptype[1].setPcell(20);
        ptype[1].calcTotalNumber(env);
        ptype[1].calcSize(env);
    }

    double getSizeOfSuperParticle(int nr, double density, const double dx){
        return pow(dx, 2.0) * density / nr;
    }

    void setMPIInfoToEnv(Environment* env) {
        int rank = MPI::Environment::rank;
        int numprocs = MPI::Environment::numprocs;

        if(numprocs != env->proc_x * env->proc_y * env->proc_z) {
            if(rank == 0) {
                cout << format("[ERROR] Allocated Process Number [%d] is different from [%d] inputted from json.") % numprocs % (env->proc_x * env->proc_y * env->proc_z) << endl;
            }
            MPI::Environment::exitWithFinalize(0);
        }

        env->numprocs = numprocs;
        env->rank = rank;

        int xrank = -1, yrank = -1, zrank = -1;

        // Processの積み方は
        // x->y->z
        for(int k = 0; k < env->proc_z; ++k){
            for(int j = 0; j < env->proc_y; ++j){
                for(int i = 0; i < env->proc_x; ++i){
                    if( (i + env->proc_x * j + env->proc_x * env->proc_y * k) == rank ){
                        xrank = i;
                        yrank = j;
                        zrank = k;
                    }
                }
            }
        }

        env->xrank = xrank;
        env->yrank = yrank;
        env->zrank = zrank;

        env->onLowXedge = (xrank == 0); env->onHighXedge = (xrank == env->proc_x - 1);
        env->onLowYedge = (yrank == 0); env->onHighYedge = (yrank == env->proc_y - 1);
        env->onLowZedge = (zrank == 0); env->onHighZedge = (zrank == env->proc_z - 1);

        // 0のノードをルートとして扱う
        env->isRootNode = (rank == 0);

        if(env->isRootNode) cout << format("[MPIINFO] allocated processes: %d") % numprocs << endl;
    }

    Environment* loadEnvironment(picojson::object& inputs){
        Environment* env = new Environment;
        auto env_inputs = inputs["Environment"].get<picojson::object>();

        for(auto it = env_inputs.begin();
                 it != env_inputs.end();
                 ++it){
            // string で switch したい...
            if(it->first == "nx"){
                env->nx = static_cast<int>(it->second.get<double>());
            } else if(it->first == "ny"){
                env->ny = static_cast<int>(it->second.get<double>());
            } else if(it->first == "nz"){
                env->nz = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_x"){
                env->proc_x = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_y"){
                env->proc_y = static_cast<int>(it->second.get<double>());
            } else if(it->first == "proc_z"){
                env->proc_z = static_cast<int>(it->second.get<double>());
            } else if(it->first == "dt"){
                env->dt = it->second.get<double>();
            } else if(it->first == "dx"){
                env->dx = it->second.get<double>();
            } else if(it->first == "max_iteration"){
                env->max_iteration = static_cast<int>(it->second.get<double>());
            } else if(it->first == "job_type"){
                env->jobtype = it->second.to_str();
            } else if(it->first == "solver_type"){
                env->solver_type = it->second.to_str();
            } else if(it->first == "boundary"){
                env->boundary = it->second.to_str();
            } else if(it->first == "dimension"){
                env->dimension = it->second.to_str();
            } else {
                std::cout <<"Unsupportted Key [" << it->first << "] is in json." << std::endl;
            }
        }

        // 1プロセスあたりのグリッド数
        // これに2を加えた数がのりしろ分になる
        env->cell_x = env->nx/env->proc_x;
        env->cell_y = env->ny/env->proc_y;
        env->cell_z = env->nz/env->proc_z;

        return env;
    }

    ParticleType* loadParticleType(picojson::object& inputs, Environment* env){
        auto plasma_inputs = inputs["Plasma"].get<picojson::object>();
        int particle_types = 0;
        for(auto i = plasma_inputs.begin(); i != plasma_inputs.end(); ++i){
            ++particle_types;
        }
        env->num_of_particle_types = particle_types;
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

            ptype[ii].calcTotalNumber(env);
            ptype[ii].calcSize(env);

            // print particle info
            if(env->isRootNode) {
                std::cout << ptype[ii];
                Utils::printTotalMemory(ptype[ii]);
            }

            ++ii;
        }

        // add pointer to ptype
        env->ptype = ptype;

        return ptype;
    }

    Grid* initializeGrid(const Environment* env){
        return new Grid(env);
    }
}

