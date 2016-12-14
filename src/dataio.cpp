#include <tdpic.h>
#include <boost/format.hpp>
#include <fstream>

using std::cout;
using std::cin;
using std::endl;
using boost::format;

namespace IO {
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

    void outputParticlePositions(const Environment* env, const ParticleType* ptype, const ParticleArray& parray, std::string filename){
        std::ofstream ofs(filename);

        for(int id = 0; id < env->particle_types; ++id){

            ofs << "## " << ptype[id].getName() << endl;

            for(int i = 0; i < parray[id].size(); ++i){
                ofs << format("%8.3f %8.3f %8.3f") % parray[id][i].getX() % parray[id][i].getY() % parray[id][i].getZ();
                ofs << format("%8.3f %8.3f %8.3f") % parray[id][i].getVX() % parray[id][i].getVY() % parray[id][i].getVZ();
                ofs << endl;
            }

            ofs << endl << endl;
        }
    }
}
