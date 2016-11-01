#include <math.h>
#include "tdpic.h"

namespace Initializer {
    double getSizeOfSuperParticle(int nr, double density, const double dx){
        return pow(dx, 2.0) * density / nr;
    }

    Field* initializeField(const Environment* env){
        Field* field = new Field;

        const int cx = env->cell_x;
        const int cy = env->cell_y;
        const int cz = env->cell_z;
        auto extents = boost::extents[cx][cy][cz];
        threeD_array* phi = new threeD_array(extents);
        threeD_array* rho = new threeD_array(extents);

        field->setPhi(phi);
        field->setRho(rho);

        threeD_array* ex = new threeD_array(boost::extents[cx-1][cy][cz]);
        threeD_array* ey = new threeD_array(boost::extents[cx][cy-1][cz]);
        threeD_array* ez = new threeD_array(boost::extents[cx][cy][cz-1]);

        field->setEx(ex);
        field->setEy(ey);
        field->setEz(ez);

        threeD_array* bx = new threeD_array(boost::extents[cx][cy-1][cz-1]);
        threeD_array* by = new threeD_array(boost::extents[cx-1][cy][cz-1]);
        threeD_array* bz = new threeD_array(boost::extents[cx-1][cy-1][cz]);

        field->setBx(bx);
        field->setBy(by);
        field->setBz(bz);

        return field;
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

        // のりしろ1セルを両側においておく
        env->cell_x = env->nx/env->proc_x + 2;
        env->cell_y = env->ny/env->proc_y + 2;
        env->cell_z = env->nz/env->proc_z + 2;

        return env;
    }

    ParticleType* loadParticleType(picojson::object& inputs, Environment* env){
        auto plasma_inputs = inputs["Plasma"].get<picojson::object>();
        int particle_types = 0;
        for(auto i = plasma_inputs.begin(); i != plasma_inputs.end(); ++i){
            ++particle_types;
        }
        env->particle_types = particle_types;
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

        return ptype;
    }
}

