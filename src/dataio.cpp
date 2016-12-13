#include <tdpic.h>
#include <boost/format.hpp>

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

    void outputParticlePositions(const Environment* env, const ParticleArray& parray){
        cout << "-- PARTICLE POSITIONS --" << endl;
        for(int id = 0; id < env->particle_types; ++id){
            cout << "-- ID: " << format("%d") % id << " -- " << endl;
        }
    }
}
