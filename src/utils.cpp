#include "tdpic.h"

using std::cout;
using std::endl;

namespace Utils {
    static std::string computeMemory(double mem){
        std::string suffix = "B";
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "kB";
        }
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "MB";
        }
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "GB";
        }

        return std::to_string(mem) + suffix;
    }

    void print3DArray(threeD_array* data){
        double* it = data->origin();
        for ( int i = 0 ; i < data->shape()[0] ; ++i ) {
            cout << "[" << i << "]" << endl;
            for ( int j = 0 ; j < data->shape()[1] ; ++j ) {
                for ( int k = 0 ; k < data->shape()[2] ; ++k, ++it ) {
                    cout << *it << " ";
                }
            }
            cout << endl;
        }
    }

    void printTotalMemory(const ParticleType& pinfo){
        cout << "Memory Info:" << endl;
        printParticleMemory(pinfo);
    }

    void printParticleMemory(const ParticleType& pinfo){
        static const double memory_per_particle = 50.0;

        double pmem = static_cast<double>(pinfo.getTotalNumber()) * memory_per_particle;

        cout << "    Particle[" << pinfo.getTotalNumber() << "]:";
        cout << computeMemory(pmem) << endl;
    }

    boost::property_tree::ptree readInputFile(const std::string& filename){
        boost::property_tree::ptree t;
        boost::property_tree::json_parser::read_json(filename, t);
        return t;
    }
}
