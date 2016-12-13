#include <tdpic.h>
#include <fstream>
#include <boost/format.hpp>

using std::cout;
using std::endl;
using boost::format;

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

        return (format("%6.2f") % mem).str() + suffix;
    }

    void print3DArray(threeD_array* data){
        double* it = data->origin();
        for ( int i = 1 ; i < data->shape()[0] - 1 ; ++i ) {
            cout << "[x:" << i << "] " << endl;
            for ( int j = 1 ; j < data->shape()[1] - 1; ++j ) {

                    if(j == 1) {
                        cout << "     [y/z]";
                        for ( int k = 1 ; k < data->shape()[2] - 1; ++k ) {
                            cout << "[" << k << "]";
                        }
                        cout << endl;
                    }
                    cout << "     [" << j << "]  ";
                for ( int k = 1 ; k < data->shape()[2] - 1; ++k, ++it ) {
                    cout << " " << (*data)[i][j][k] << " ";
                }
                cout << endl;
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

    std::string readFile(const std::string& filename){
        std::ifstream ifs;
        std::string res;
        std::string ifs_buffer;

        ifs.open(filename, std::ios::in);

        while(!ifs.eof()){
            std::getline(ifs, ifs_buffer);
            res += ifs_buffer;
        }

        return res;

    }

    picojson::value::object readJSONFile(const std::string& filename){
        std::string json = readFile(filename);

        picojson::value v;
        std::string error = picojson::parse(v, json);
        picojson::value::object& o = v.get<picojson::object>();

        return o;
    }
}
