#include <math.h>
#include <tdpic.h>

namespace Initializer {
    double getSizeOfSuperParticle(int nr, double density, const double dx){
        return pow(dx, 2.0) * density / nr;
    }

    void initializeRootField(const Environment* env, Grid* grid){
        Field* field = new Field;

        const int cx = env->cell_x + 2;
        const int cy = env->cell_y + 2;
        const int cz = env->cell_z + 2;

        threeDArray phi = Utils::create3DArray(cx, cy, cz);
        threeDArray rho = Utils::create3DArray(cx, cy, cz);
        field->setPhi(phi);
        field->setRho(rho);

        threeDArray ex = Utils::create3DArray(cx-1, cy, cz);
        threeDArray ey = Utils::create3DArray(cx, cy-1, cz);
        threeDArray ez = Utils::create3DArray(cx, cy, cz-1);

        field->setEx(ex);
        field->setEy(ey);
        field->setEz(ez);

        threeDArray bx = Utils::create3DArray(cx, cy-1, cz-1);
        threeDArray by = Utils::create3DArray(cx-1, cy, cz-1);
        threeDArray bz = Utils::create3DArray(cx-1, cy-1, cz);

        field->setBx(bx);
        field->setBy(by);
        field->setBz(bz);

        grid->setField(field);
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
            std::cout << ptype[ii];
            Utils::printTotalMemory(ptype[ii]);
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

